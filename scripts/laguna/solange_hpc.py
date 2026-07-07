#!/usr/bin/env python3
"""
solange_hpc.py — SOLANGE HPC runner for the Laguna cluster (SLURM + GPU).

Scales the Phase-3A active space beyond the 4-qubit laptop proxy:
  RHF -> CASSCF(ncas, nelecas)/<basis>  ->  Jordan-Wigner Hamiltonian
       ->  exact active-space ground state (reference)
       -> (optional) statevector VQE on GPU (PennyLane lightning.gpu / default.qubit)

It emits TWO artifacts, both in the exact schema SOLANGE already consumes:

  1) <out>/jw_<KEY>.json          -> a jw_hamiltonians.json entry (Hamiltonian).
                                     Merge into backend/jw_hamiltonians.json + root.
  2) <out>/provenance_<KEY>.json  -> a signed P1-P9 record. The P8 seal is computed
                                     with the SAME sha256 algorithm as the SOLANGE
                                     backend (build_p8_seal), so the platform can
                                     re-verify integrity WITHOUT a GPU.

Design note (honesty): a minimal STO-3G model compound does not contain enough
orbitals for a large active space. Growing the active space means growing BOTH the
active space (--ncas/--nelecas) AND the basis (--basis: sto-3g -> 6-31g -> cc-pvdz).
This is a scientific choice, exposed as flags — not a silent parameter bump.

Compound geometries are imported from generate_expansion_jw.py (same repo), so this
runner stays consistent with the 31 Hamiltonians already in the platform.

USAGE (validate first with a tiny dry-run that reproduces the known 4q result):
  python solange_hpc.py --compound acetamide --basis sto-3g --ncas 2 --nelecas 2 \
      --key ARID2_LOF --side native --residue "Gln1118 dry-run" --out ./out
  # then scale up, e.g.:
  python solange_hpc.py --compound acetamide --basis 6-31g --ncas 8 --nelecas 8 \
      --key ARID2_LOF --side native --residue "Gln1118 amide" --vqe --out ./out
"""

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Import the existing, validated geometry table so we do not duplicate coordinates.
_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[2]          # scripts/laguna/solange_hpc.py -> repo root
sys.path.insert(0, str(_REPO_ROOT))
try:
    from generate_expansion_jw import GEOM, geom_to_string
except Exception as e:  # pragma: no cover
    print(f"ERROR: could not import GEOM from generate_expansion_jw.py "
          f"(expected at {_REPO_ROOT}). {e}", file=sys.stderr)
    raise


# ── Hardware autodetection ───────────────────────────────────────────────────

def detect_gpu(safety_gb: float = 8.0):
    """Return (gpu_name, vram_mb, max_qubits) or (None, None, None) if no GPU.

    Statevector needs 2**n complex128 = 2**n * 16 bytes, plus work buffers.
    We budget ~3x the statevector and reserve `safety_gb` for the framework.
    """
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name,memory.total",
             "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL, timeout=30).decode().strip()
    except Exception:
        return None, None, None
    # take the largest GPU visible
    best_name, best_mb = None, 0
    for line in out.splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 2:
            continue
        name, mb = parts[0], int(float(parts[1]))
        if mb > best_mb:
            best_name, best_mb = name, mb
    if best_mb == 0:
        return None, None, None
    usable_bytes = best_mb * 1024**2 - safety_gb * 1e9
    if usable_bytes <= 0:
        return best_name, best_mb, 0
    max_qubits = int(math.floor(math.log2(usable_bytes / (16 * 3))))
    return best_name, best_mb, max_qubits


# ── PySCF: RHF -> CASSCF(ncas, nelecas)/basis ────────────────────────────────

def run_casscf(compound, basis, ncas, nelecas, verbose=0):
    from pyscf import gto, scf, mcscf, ao2mo

    if compound not in GEOM:
        raise KeyError(f"Unknown compound '{compound}'. "
                       f"Available: {', '.join(sorted(GEOM))}")

    mol = gto.Mole()
    mol.atom    = geom_to_string(GEOM[compound])
    mol.basis   = basis
    mol.charge  = 0
    mol.spin    = 0
    mol.verbose = verbose
    mol.build()

    n_ao = mol.nao_nr()
    if ncas > n_ao:
        raise ValueError(
            f"ncas={ncas} exceeds available orbitals ({n_ao}) for {compound}/{basis}. "
            f"Use a larger basis (e.g. 6-31g, cc-pvdz) or a smaller active space.")

    mf = scf.RHF(mol)
    mf.max_cycle = 400
    mf.conv_tol  = 1e-10
    e_rhf = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"RHF did not converge for {compound}/{basis}")

    mc = mcscf.CASSCF(mf, ncas=ncas, nelecas=nelecas)
    mc.conv_tol        = 1e-9
    mc.conv_tol_grad   = 1e-5
    mc.max_cycle_macro = 200
    e_casscf = mc.kernel()[0]
    if not mc.converged:
        print(f"  WARNING: CASSCF did not fully converge for {compound} "
              f"CAS({nelecas},{ncas})/{basis} (JW terms still valid)", file=sys.stderr)

    if not mc.converged:
        # A non-converged CASSCF yields integrals inconsistent with e_casscf —
        # exactly the ARID2 failure. Retry with the second-order (Newton) solver.
        print(f"  CASSCF slow to converge — retrying with second-order (Newton) solver...",
              file=sys.stderr)
        mc = mcscf.CASSCF(mf, ncas=ncas, nelecas=nelecas).newton()
        mc.max_cycle_macro = 300
        e_casscf = mc.kernel()[0]
    if not mc.converged:
        raise RuntimeError(
            f"CASSCF failed to converge for {compound}/{basis} CAS({nelecas},{ncas}) "
            f"even with the Newton solver — refusing to emit an inconsistent Hamiltonian.")

    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), mc.ncas)

    # Consistency gate: active-space FCI of (h1e,h2e) must equal e_casscf - ecore.
    from pyscf import fci
    na = nelecas // 2
    e_fci = fci.direct_spin1.FCI().kernel(h1e, h2e, ncas, (na, nelecas - na), ecore=0.0)[0]
    if abs((ecore + e_fci) - e_casscf) > 1e-3:
        raise RuntimeError(
            f"{compound}: active-space FCI ({ecore + e_fci:.6f} Ha) != e_casscf "
            f"({e_casscf:.6f} Ha) — Hamiltonian inconsistent, refusing to emit.")

    return {
        "e_rhf": float(e_rhf), "e_casscf": float(e_casscf), "ecore": float(ecore),
        "e_fci_active": float(e_fci),
        "h1e": h1e, "h2e": h2e, "ncas": int(ncas), "nelecas": int(nelecas),
        "n_ao": int(n_ao), "converged": bool(mc.converged),
    }


# ── Jordan-Wigner build (generalised to N orbitals) ──────────────────────────
# Convention identical to generate_expansion_jw.build_jw_terms, with n_orb=ncas.
# Spin-orbital mapping: index 2*p = orbital p (alpha), 2*p+1 = orbital p (beta).

def build_jw_terms(h1e, h2e, ncas, exact=True):
    from openfermion import InteractionOperator, jordan_wigner
    from openfermion.linalg import get_sparse_operator
    import scipy.sparse.linalg as spla

    n_orb = ncas
    dim = 2 * n_orb
    one_body = np.zeros((dim, dim))
    two_body = np.zeros((dim, dim, dim, dim))

    # PySCF get_h2eff is CHEMIST-ordered (pq|rs); OpenFermion's two_body_tensor
    # is PHYSICIST-ordered for a†_p a†_q a_r a_s. Correct remap: g = h2e(0,2,3,1),
    # paired with the ab/ba spin pattern below. Verified end-to-end (2e-sector
    # ground state == e_casscf - ecore for both symmetric and asymmetric compounds).
    g = np.transpose(np.asarray(h2e), (0, 2, 3, 1))

    for p in range(n_orb):
        for q in range(n_orb):
            one_body[2*p,   2*q]   = h1e[p, q]      # alpha
            one_body[2*p+1, 2*q+1] = h1e[p, q]      # beta

    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    val = 0.5 * g[p, q, r, s]
                    two_body[2*p,   2*q,   2*r,   2*s]   = val   # aa
                    two_body[2*p+1, 2*q+1, 2*r+1, 2*s+1] = val   # bb
                    two_body[2*p,   2*q+1, 2*r+1, 2*s]   = val   # ab
                    two_body[2*p+1, 2*q,   2*r,   2*s+1] = val   # ba

    ham_op = InteractionOperator(constant=0.0,
                                 one_body_tensor=one_body,
                                 two_body_tensor=two_body)
    jw_op = jordan_wigner(ham_op)

    e_active_exact = None
    if exact:
        # Sparse Lanczos is fine up to ~16 qubits; guard larger sizes.
        if dim <= 16:
            sparse = get_sparse_operator(jw_op)
            e_active_exact = float(spla.eigsh(sparse, k=1, which="SA")[0][0])
        else:
            print(f"  NOTE: skipping exact diagonalisation at {dim} qubits "
                  f"(too large for sparse Lanczos here); rely on VQE/CASSCF.",
                  file=sys.stderr)

    terms = []
    for term, coeff in jw_op.terms.items():
        if abs(coeff) < 1e-12:
            continue
        pauli = "I" if len(term) == 0 else " ".join(
            f"{op}{qubit}" for qubit, op in sorted(term))
        terms.append({"coeff": float(coeff.real), "pauli": pauli})

    return terms, e_active_exact, jw_op


# ── Optional statevector VQE on GPU ──────────────────────────────────────────

def run_vqe(terms, n_qubits, nelec, steps=80):
    """Statevector VQE. Prefers lightning.gpu, falls back to default.qubit.
    Returns dict(energy_ha, variance, elapsed_s, device, gate_count, depth, steps).
    """
    import pennylane as qml
    from pennylane import numpy as pnp

    # Build the PennyLane Hamiltonian from JW terms.
    coeffs, ops = [], []
    pauli_map = {"X": qml.PauliX, "Y": qml.PauliY, "Z": qml.PauliZ}
    for t in terms:
        coeffs.append(t["coeff"])
        if t["pauli"] == "I":
            ops.append(qml.Identity(0))
        else:
            factors = []
            for tok in t["pauli"].split():
                p, wire = tok[0], int(tok[1:])
                factors.append(pauli_map[p](wire))
            op = factors[0]
            for f in factors[1:]:
                op = op @ f
            ops.append(op)
    H = qml.Hamiltonian(coeffs, ops)

    for dev_name in ("lightning.gpu", "lightning.qubit", "default.qubit"):
        try:
            dev = qml.device(dev_name, wires=n_qubits)
            break
        except Exception:
            continue
    print(f"  VQE device: {dev.name}")

    singles, doubles = qml.qchem.excitations(nelec, n_qubits)
    hf = qml.qchem.hf_state(nelec, n_qubits)

    @qml.qnode(dev, diff_method="adjoint")
    def circuit(params):
        qml.AllSinglesDoubles(params, range(n_qubits), hf, singles, doubles)
        return qml.expval(H)

    n_params = len(singles) + len(doubles)
    params = pnp.zeros(n_params, requires_grad=True)
    opt = qml.AdamOptimizer(stepsize=0.1)

    # Adaptive loop: `steps` is a MAX. When the energy stops improving we shrink
    # the stepsize (settles into the minimum → closes the last mHa to chemical
    # accuracy); when the stepsize is tiny AND progress has stalled we stop early
    # (no wasted iterations on a plateau). Fixes both the ~4 mHa plateau and the
    # "400 steps but converged by ~120" waste.
    t0 = time.time()
    energies = []
    best = float("inf")
    best_params = params
    stall = 0
    for i in range(steps):
        params, e = opt.step_and_cost(circuit, params)
        e = float(e)
        energies.append(e)
        if i % 20 == 0:
            print(f"    step {i:3d}  E = {e:.6f} Ha  (stepsize {opt.stepsize:.3g})")
        if e < best - 1e-8:
            best = e
            best_params = params.copy()      # keep the BEST, not the last
            stall = 0
        else:
            stall += 1
            if stall >= 25:                 # progress stalled → refine
                opt.stepsize *= 0.5
                stall = 0
                if opt.stepsize < 5e-4:      # already at the minimum → stop
                    print(f"    converged at step {i} (best E = {best:.8f} Ha)")
                    break
    elapsed = time.time() - t0
    params = best_params                     # report the best variational estimate
    tail = energies[-20:] if len(energies) >= 20 else energies
    variance = float(np.var(tail))

    # Circuit-size metadata for P1 provenance. qml.specs() returns different
    # shapes across PennyLane versions (dict vs Specs/CircuitSpecs object) — keep
    # it best-effort so a metadata quirk never aborts a completed VQE run.
    gate_count, depth = 0, 0
    try:
        specs = qml.specs(circuit)(params)
        res = specs["resources"] if isinstance(specs, dict) else getattr(specs, "resources", None)
        if res is not None:
            gate_count = int(getattr(res, "num_gates", 0))
            depth = int(getattr(res, "depth", 0))
    except Exception:
        pass

    return {
        "energy_ha": float(best), "variance": variance,
        "elapsed_s": round(elapsed, 3), "device": str(dev.name),
        "gate_count": gate_count, "depth": depth, "steps": len(energies),
        "convergence": energies,
    }


# ── P1-P9 provenance (seal matches SOLANGE backend build_p8_seal) ────────────

def build_p8_seal(record: dict) -> str:
    seal_payload = json.dumps(
        {k: v for k, v in record.items()
         if k.startswith(("p1_", "p2_", "p3_", "p4_",
                          "p5_", "p6_", "p7_", "p9_"))},
        sort_keys=True, default=str)
    return hashlib.sha256(seal_payload.encode()).hexdigest()


def build_provenance(args, cas, jw_terms, e_active_exact, vqe, gpu_name, vram_mb):
    now = datetime.now(timezone.utc).isoformat()
    n_qubits = 2 * cas["ncas"]
    energy = vqe["energy_ha"] if vqe else (e_active_exact if e_active_exact is not None
                                           else cas["e_casscf"])
    method = "statevector VQE (AllSinglesDoubles UCCSD)" if vqe else \
             "exact active-space diagonalisation"
    backend = (f"{gpu_name} · {vram_mb} MiB · "
               f"{vqe['device'] if vqe else 'CPU sparse'}") if gpu_name else "CPU"

    record = {
        # identity / non-sealed context
        "mutation_id":       args.key,
        "side":              args.side,
        "provenance_source": "HPC/Laguna",
        "created_at":        now,

        # P1 — circuit
        "p1_qubit_count":  n_qubits,
        "p1_gate_count":   vqe["gate_count"] if vqe else 0,
        "p1_depth":        vqe["depth"] if vqe else 0,
        "p1_ansatz":       ("AllSinglesDoubles UCCSD (HPC statevector)"
                            if vqe else "none — exact diagonalisation"),
        "p1_circuit_hash": hashlib.sha256(json.dumps(
            {"key": args.key, "side": args.side, "ncas": cas["ncas"],
             "nelecas": cas["nelecas"], "basis": args.basis},
            sort_keys=True).encode()).hexdigest(),

        # P2 — compilation
        "p2_compiler":         "PySCF + openfermion",
        "p2_encoding":         f"Jordan-Wigner (CASSCF({cas['nelecas']},{cas['ncas']}) "
                               f"-> {len(jw_terms)} Pauli terms)",
        "p2_basis_set":        f"{args.basis} (CASSCF({cas['nelecas']},{cas['ncas']}))",
        "p2_active_electrons": cas["nelecas"],
        "p2_active_orbitals":  cas["ncas"],
        "p2_model_compound":   args.compound,

        # P3 — backend (transparent about where it ran)
        "p3_backend":           backend,
        "p3_calibration_epoch": now,
        "p3_simulator":         True,

        # P4 — noise (exact statevector -> none)
        "p4_gate_error_rate":    0.0,
        "p4_readout_error_rate": 0.0,
        "p4_note":               "Exact statevector on HPC GPU — noiseless; "
                                 "Phase 3B records real IBM Heron r3 calibration.",

        # P5 — execution
        "p5_shots":            None,
        "p5_raw_energy":       energy,
        "p5_energy_variance":  vqe["variance"] if vqe else 0.0,
        "p5_opt_steps":        vqe["steps"] if vqe else 0,
        "p5_elapsed_s":        vqe["elapsed_s"] if vqe else None,
        "p5_ecore_ha":         cas["ecore"],
        "p5_casscf_ref_ha":    cas["e_casscf"],

        # P6 — mitigation
        "p6_method": "none — exact statevector",
        "p6_note":   "Phase 3B: ZNE + Pauli Twirling on IBM Heron r3",

        # P7 — result
        "p7_energy_ha":  energy,
        "p7_method":     method,
        "p7_confidence": 0.95,

        # P9 — decoder (not applicable to statevector)
        "p9_applicable": False,
        "p9_note":       "P9 applies when the QEC decoder is active (Phase 3B hardware).",
    }
    record["p8_hash"]      = build_p8_seal(record)
    record["p8_algorithm"] = "SHA-256 over P1-P7,P9 (sort_keys, default=str)"
    record["p8_sealed_at"] = now
    return record


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="SOLANGE HPC runner (Laguna).")
    ap.add_argument("--compound", required=True, help="model compound (key in GEOM)")
    ap.add_argument("--basis", default="6-31g", help="sto-3g | 6-31g | cc-pvdz | ...")
    ap.add_argument("--ncas", type=int, default=None, help="active orbitals (auto if omitted)")
    ap.add_argument("--nelecas", type=int, default=None, help="active electrons (defaults to 2*ncas//... )")
    ap.add_argument("--key", required=True, help="jw_hamiltonians key, e.g. ARID2_LOF")
    ap.add_argument("--side", default="native", choices=["native", "mutant"])
    ap.add_argument("--residue", default="", help="residue_note for the JW entry")
    ap.add_argument("--vqe", action="store_true", help="also run statevector VQE (GPU)")
    ap.add_argument("--vqe-steps", type=int, default=80)
    ap.add_argument("--out", default="./out", help="output directory")
    ap.add_argument("--verbose", type=int, default=0)
    args = ap.parse_args()

    Path(args.out).mkdir(parents=True, exist_ok=True)

    gpu_name, vram_mb, max_qubits = detect_gpu()
    print("=" * 68)
    print(f"SOLANGE HPC · key={args.key}/{args.side} · compound={args.compound} "
          f"· basis={args.basis}")
    if gpu_name:
        print(f"GPU: {gpu_name} · {vram_mb} MiB · statevector ceiling ≈ {max_qubits} qubits "
              f"(≈ CAS({max_qubits}) orbitals)")
    else:
        print("GPU: none detected — CPU path (exact diag limited to ~16 qubits)")

    # Choose active space: honour --ncas, else auto from GPU ceiling (capped at 8 for safety).
    if args.ncas is None:
        ceiling_orbs = (max_qubits // 2) if max_qubits else 8
        args.ncas = max(2, min(ceiling_orbs, 8))
        print(f"--ncas not given → auto-selected ncas={args.ncas} "
              f"(override with --ncas for larger runs)")
    if args.nelecas is None:
        args.nelecas = args.ncas          # half-filled default; override for real chemistry
        print(f"--nelecas not given → default nelecas={args.nelecas}")

    n_qubits = 2 * args.ncas
    if gpu_name and max_qubits and n_qubits > max_qubits:
        print(f"ERROR: requested {n_qubits} qubits exceeds GPU ceiling {max_qubits}. "
              f"Reduce --ncas to ≤ {max_qubits // 2}.", file=sys.stderr)
        sys.exit(2)

    print(f"Active space: CAS({args.nelecas},{args.ncas}) → {n_qubits} qubits")
    print("-" * 68)

    cas = run_casscf(args.compound, args.basis, args.ncas, args.nelecas, args.verbose)
    print(f"RHF     E = {cas['e_rhf']:.8f} Ha  ({cas['n_ao']} AOs)")
    print(f"CASSCF  E = {cas['e_casscf']:.8f} Ha  (ecore {cas['ecore']:.8f})")

    terms, e_active_exact, _ = build_jw_terms(cas["h1e"], cas["h2e"], args.ncas,
                                              exact=(n_qubits <= 16))
    print(f"JW: {len(terms)} Pauli terms on {n_qubits} qubits"
          + (f" · exact active E = {e_active_exact:.8f} Ha" if e_active_exact is not None else ""))

    vqe = None
    if args.vqe:
        print("Running statevector VQE ...")
        vqe = run_vqe(terms, n_qubits, args.nelecas, steps=args.vqe_steps)
        print(f"VQE     E = {vqe['energy_ha']:.8f} Ha  ({vqe['elapsed_s']}s, {vqe['device']})")
        # exact reference: openfermion diag (<=16q) else the pyscf FCI from the
        # consistency gate (available at every size) — so ΔE prints beyond 16 qubits.
        exact_ref = e_active_exact if e_active_exact is not None else cas.get("e_fci_active")
        if exact_ref is not None:
            print(f"exact active E = {exact_ref:.8f} Ha  (reference)")
            print(f"ΔE(VQE−exact) = {(vqe['energy_ha']-exact_ref)*1000:.3f} mHa")

    # e_active_rhf: active-space HF reference = total RHF − frozen-core energy.
    e_active_rhf = cas["e_rhf"] - cas["ecore"]

    # ── Artifact 1: jw_hamiltonians.json entry (schema-compatible) ──
    jw_entry = {
        "compound":       args.compound,
        "residue_note":   args.residue,
        "ecore":          cas["ecore"],
        "e_casscf":       cas["e_casscf"],
        "e_active_exact": cas.get("e_fci_active",
                                  e_active_exact if e_active_exact is not None else cas["e_casscf"]),
        "e_active_rhf":   float(e_active_rhf),
        "n_paulis":       len(terms),
        "terms":          terms,
        # provenance stamps (extra keys; SOLANGE reads the ones above)
        "basis":          args.basis,
        "ncas":           args.ncas,
        "nelecas":        args.nelecas,
        "source":         "HPC/Laguna",
    }
    jw_path = Path(args.out) / f"jw_{args.key}_{args.side}.json"
    jw_path.write_text(json.dumps({args.key: {args.side: jw_entry}}, indent=2))
    print(f"WROTE {jw_path}")

    # ── Artifact 2: signed P1-P9 provenance ──
    prov = build_provenance(args, cas, terms, e_active_exact, vqe, gpu_name, vram_mb)
    prov_path = Path(args.out) / f"provenance_{args.key}_{args.side}.json"
    prov_path.write_text(json.dumps(prov, indent=2, default=str))
    print(f"WROTE {prov_path}")
    print(f"P8 seal: {prov['p8_hash'][:16]}…  "
          f"(SOLANGE re-verifies this later by recomputing the hash — that check needs no GPU)")
    print("=" * 68)
    print("NOTE: these two files are on Laguna ONLY — they are NOT in SOLANGE yet and")
    print("nothing has been pushed. To get them into SOLANGE, either:")
    print(f"  • POST them to the /api/simulate/hpc endpoint (dynamic — coming), or")
    print(f"  • send them to sync manually:")
    print(f"      {jw_path}")
    print(f"      {prov_path}")


if __name__ == "__main__":
    main()
