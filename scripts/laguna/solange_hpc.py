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

_SEAL_EXCLUDE = {"p3_calibration_epoch"}   # not round-trip-safe in the DB (timestamptz)


def build_p8_payload(record: dict) -> str:
    """The exact canonical JSON string the P8 seal hashes over (P1–P7 + P9)."""
    return json.dumps(
        {k: v for k, v in record.items()
         if k.startswith(("p1_", "p2_", "p3_", "p4_",
                          "p5_", "p6_", "p7_", "p9_"))
         and k not in _SEAL_EXCLUDE},
        sort_keys=True, default=str)


def build_p8_seal(record: dict) -> str:
    return hashlib.sha256(build_p8_payload(record).encode()).hexdigest()


def build_provenance(args, cas, jw_terms, e_active_exact, vqe, gpu_name, vram_mb,
                     elapsed_total=None):
    now = datetime.now(timezone.utc).isoformat()
    n_qubits = 2 * cas["ncas"]
    # p7 must be a TOTAL energy (ecore + active) to match e_casscf (total) and the
    # live 4q runs' convention — VQE/exact solve the ACTIVE Hamiltonian (ecore=0),
    # so add ecore back. Storing the active energy made ΔE nonsense in the panel.
    if vqe:
        energy = cas["ecore"] + vqe["energy_ha"]
    elif e_active_exact is not None:
        energy = cas["ecore"] + e_active_exact
    else:
        energy = cas["e_casscf"]
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
        "p2_compiler_version": None,
        "p2_encoding":         f"Jordan-Wigner (CASSCF({cas['nelecas']},{cas['ncas']}) "
                               f"-> {len(jw_terms)} Pauli terms)",
        "p2_basis_set":        f"{args.basis} (CASSCF({cas['nelecas']},{cas['ncas']}))",
        "p2_active_electrons": cas["nelecas"],
        "p2_active_orbitals":  cas["ncas"],
        "p2_model_compound":   args.compound,
        # The full source Hamiltonian — every Pauli term + coefficient. P2's own
        # definition ("compilation lineage") already names this; it used to be
        # written only to the local jw_*.json file and never sent to SOLANGE.
        # Named with the p2_ prefix so LEON's existing P8 seal covers it too.
        "p2_jw_terms": jw_terms,

        # P3 — backend (transparent about where it ran)
        "p3_backend":           backend,
        "p3_backend_version":   None,
        "p3_calibration_epoch": now,
        "p3_simulator":         True,

        # P4 — noise (exact statevector -> none)
        "p4_gate_error_rate":    0.0,
        "p4_readout_error_rate": 0.0,
        "p4_t1_us":              None,
        "p4_t2_us":              None,
        "p4_note":               "Exact statevector on HPC GPU — noiseless; "
                                 "Phase 3B records real IBM Heron r3 calibration.",

        # P5 — execution
        "p5_shots":            None,
        "p5_raw_energy":       energy,
        "p5_energy_variance":  vqe["variance"] if vqe else 0.0,
        "p5_opt_steps":        vqe["steps"] if vqe else 0,
        "p5_elapsed_s":        elapsed_total,   # TOTAL compute time = the classical-wall cost
        "p5_ecore_ha":         cas["ecore"],
        "p5_active_energy_ha": None,
        "p5_casscf_ref_ha":    cas["e_casscf"],

        # P6 — mitigation
        "p6_method": "none — exact statevector",
        "p6_note":   "Phase 3B: ZNE + Pauli Twirling on IBM Heron r3",

        # P7 — result
        "p7_energy_ha":  energy,
        "p7_ci_lower":   None,
        "p7_ci_upper":   None,
        "p7_method":     method,
        "p7_confidence": 0.95,

        # P9 — decoder (not applicable to statevector)
        "p9_applicable": False,
        "p9_note":       "P9 applies when the QEC decoder is active (Phase 3B hardware).",
    }
    record["p8_seal_payload"] = build_p8_payload(record)   # stored verbatim for robust re-verify
    record["p8_hash"]      = hashlib.sha256(record["p8_seal_payload"].encode()).hexdigest()
    record["p8_algorithm"] = "SHA-256 over P1-P7,P9 (sort_keys, default=str)"
    record["p8_sealed_at"] = now
    return record


# ── Main ─────────────────────────────────────────────────────────────────────

# ── Pull agent: run INSIDE the Laguna session; pull jobs from SOLANGE's queue ──
# This is the outbound half of the loop. Cluster security (Duo 2FA, no inbound)
# blocks SOLANGE from pushing, so the agent reaches OUT to the queue, runs the job
# here (with full in-session privileges), and --submits the verified result back.
_SUPA_URL  = "https://lzzuxtnubznrkxwxjaab.supabase.co"
_SUPA_ANON = ("eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imx6"
              "enV4dG51Ynpucmt4d3hqYWFiIiwicm9sZSI6ImFub24iLCJpYXQiOjE3Nzg5MTQ0MzgsImV4cCI6"
              "MjA5NDQ5MDQzOH0.fKApXc3ZPHXh4O008A5oFE5vbTNqJ168AI9NzIl4vHA")


def _supabase_login(email, password):
    import urllib.request
    req = urllib.request.Request(
        _SUPA_URL + "/auth/v1/token?grant_type=password",
        data=json.dumps({"email": email, "password": password}).encode(),
        method="POST",
        headers={"apikey": _SUPA_ANON, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())["access_token"]


def _resolve_compound(key, side):
    """Map a mutation key (+side) to its model compound (a GEOM key)."""
    try:
        from generate_expansion_jw import EXPANSION_MODELS
    except Exception:
        EXPANSION_MODELS = {}
    if key in EXPANSION_MODELS:
        entry = EXPANSION_MODELS[key]
        if side in entry:
            return entry[side]["compound"]
        # The key IS a known target — this is a DIFFERENT failure than "unknown
        # gene": the requested side was deliberately left unassigned (e.g. a
        # generic LOF target's mutant compound is scientifically ambiguous), not
        # forgotten. Say so, instead of the same generic error as a typo'd key.
        raise ValueError(
            f"{key} is a known expansion target, but its '{side}' side has no "
            f"model compound defined (only {sorted(set(entry) & {'native', 'mutant'})} "
            f"available) — a generic LOF target often models only the native/"
            f"wild-type residue; a specific point mutation (like TP53_C275F) is "
            f"needed for a scientifically grounded mutant compound.")
    core = {  # core report genes not in EXPANSION_MODELS
        "TP53_C275F": {"native": "methanethiol", "mutant": "toluene"},
        "KEAP1_LOF":  {"native": "formamide",    "mutant": "methanethiol"},
        "STK11_LKB1": {"native": "acetic_acid",  "mutant": "acetamide"},
    }
    if key in core:
        return core[key][side]
    raise ValueError(f"cannot resolve model compound for {key}/{side}")


def _db_status_from(stdout):
    """Pull the db=<status> token the --submit subprocess prints, for diagnostics."""
    for tok in stdout.split():
        if tok.startswith("db="):
            return tok[3:]
    return "?"


def _run_id_from(stdout):
    """Pull the run_id=<uuid> token the --submit subprocess prints, so the agent can
    report it back to /hpc/dispatch/{id}/status — this is what links a dispatch
    queue entry to its resulting simulation_runs row (hpc_dispatch.run_id)."""
    for tok in stdout.split():
        if tok.startswith("run_id="):
            return tok[len("run_id="):]
    return None


def _p8_prefix_from(stdout):
    """Pull the P8 hash prefix the --submit subprocess prints (the '  P8 = xxxx…'
    line), so it can be echoed in the DONE line — lets you cross-reference this
    exact terminal run against the SOLANGE HPC Ladder's P8 column by eye."""
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith("P8 = "):
            return line[len("P8 = "):].split("…")[0].strip()
    return None


def _ts():
    """Wall-clock timestamp for agent log lines — lets you read elapsed time
    directly off a running JupyterLab/terminal session (or a tail -f'd log file)
    without cross-referencing OS timestamps."""
    return datetime.now().strftime("%H:%M:%S")


def _post_status(api, hdr, did, status, note=None, run_id=None):
    import urllib.request
    body = {"status": status}
    if note:   body["note"] = str(note)[:400]
    if run_id: body["run_id"] = run_id
    req = urllib.request.Request(
        api.rstrip("/") + f"/api/simulate/hpc/dispatch/{did}/status",
        data=json.dumps(body).encode(), method="POST",
        headers={**hdr, "Content-Type": "application/json"})
    try:
        urllib.request.urlopen(req, timeout=30).read()
    except Exception as e:
        print(f"[{_ts()}] [agent] status post failed: {e}", file=sys.stderr)


def _post_heartbeat(api, hdr):
    """Tell SOLANGE this agent is alive (best-effort; never blocks the loop)."""
    import urllib.request
    try:
        req = urllib.request.Request(
            api.rstrip("/") + "/api/simulate/hpc/agent/heartbeat",
            data=json.dumps({"agent": "laguna"}).encode(), method="POST",
            headers={**hdr, "Content-Type": "application/json"})
        urllib.request.urlopen(req, timeout=15).read()
    except Exception:
        pass


def run_agent(api, poll_s, token, out_dir):
    import urllib.request
    hdr = {"Authorization": "Bearer " + token}
    print(f"[{_ts()}] [agent] up · polling {api} every {poll_s}s · Ctrl+C to stop")
    while True:
        _post_heartbeat(api, hdr)   # liveness ping so SOLANGE shows the agent online
        try:
            req = urllib.request.Request(
                api.rstrip("/") + "/api/simulate/hpc/dispatch/next", headers=hdr)
            with urllib.request.urlopen(req, timeout=30) as r:
                job = json.loads(r.read().decode()).get("job")
        except Exception as e:
            print(f"[{_ts()}] [agent] poll error: {e}", file=sys.stderr)
            time.sleep(poll_s); continue
        if not job:
            time.sleep(poll_s); continue
        did = job["id"]
        t_job_start = time.time()
        print(f"[{_ts()}] [agent] job {did[:8]} · {job['key']}/{job.get('side','native')} "
              f"CAS({job['nelecas']},{job['ncas']}) vqe={job.get('run_vqe')}")
        try:
            compound = job.get("compound") or _resolve_compound(job["key"], job.get("side", "native"))
            cmd = [sys.executable, str(_HERE),
                   "--compound", compound, "--basis", job.get("basis", "6-31g"),
                   "--ncas", str(job["ncas"]), "--nelecas", str(job["nelecas"]),
                   "--key", job["key"], "--side", job.get("side", "native"),
                   "--residue", job.get("residue") or "agent",
                   "--out", out_dir, "--submit", api]
            if job.get("run_vqe"):
                cmd += ["--vqe", "--vqe-steps", "200"]
            res = subprocess.run(cmd, capture_output=True, text=True)
            # Don't mark DONE unless the backend actually STORED the run — otherwise a
            # run can pass verification yet never appear in SOLANGE (e.g. the
            # p8_seal_payload migration not run). "stored" or "stored_no_payload" both count.
            passed = "status=PASSED" in res.stdout
            stored = ("db=stored" in res.stdout or "db=stored_no_payload" in res.stdout)
            ok = (res.returncode == 0 and passed and stored)
            if passed and not stored:
                note = "verified but NOT stored (db=%s)" % _db_status_from(res.stdout)
            else:
                note = "ok" if ok else (res.stderr[-300:] or res.stdout[-300:])
            run_id = _run_id_from(res.stdout) if ok else None
            p8 = _p8_prefix_from(res.stdout) if ok else None
            _post_status(api, hdr, did, "done" if ok else "failed", note=note, run_id=run_id)
            job_elapsed = time.time() - t_job_start
            # Cross-reference this run against the SOLANGE HPC Ladder: the P8 prefix
            # printed here is the SAME hash shown (truncated) in the ladder's P8
            # column — match them by eye, no need to dig through raw logs.
            print(f"[{_ts()}] [agent] job {did[:8]} → {'DONE' if ok else 'FAILED'}"
                  f"  [{job_elapsed/60:.1f}m]"
                  + (f"  P8={p8}" if ok and p8 else "")
                  + ("" if ok else f"  ({note[:80]})"))
        except KeyboardInterrupt:
            raise
        except Exception as e:
            _post_status(api, hdr, did, "failed", note=str(e))
            print(f"[{_ts()}] [agent] job {did[:8]} error: {e}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="SOLANGE HPC runner (Laguna).")
    ap.add_argument("--compound", default=None, help="model compound (key in GEOM)")
    ap.add_argument("--basis", default="6-31g", help="sto-3g | 6-31g | cc-pvdz | ...")
    ap.add_argument("--ncas", type=int, default=None, help="active orbitals (auto if omitted)")
    ap.add_argument("--nelecas", type=int, default=None, help="active electrons (defaults to 2*ncas//... )")
    ap.add_argument("--key", default=None, help="jw_hamiltonians key, e.g. ARID2_LOF")
    ap.add_argument("--side", default="native", choices=["native", "mutant"])
    ap.add_argument("--residue", default="", help="residue_note for the JW entry")
    ap.add_argument("--vqe", action="store_true", help="also run statevector VQE (GPU)")
    ap.add_argument("--vqe-steps", type=int, default=80)
    ap.add_argument("--out", default="./out", help="output directory")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    default=None, metavar="URL",
                    help="POST result+provenance to the SOLANGE backend for dynamic ingest "
                         "(re-verified there). Optional URL; defaults to the production API.")
    ap.add_argument("--verbose", type=int, default=0)
    # ── Pull-agent mode: run in the Laguna session, execute queued SOLANGE jobs ──
    ap.add_argument("--agent", action="store_true",
                    help="run as the cluster pull-agent: poll SOLANGE's dispatch queue "
                         "and execute queued jobs here, submitting results back")
    ap.add_argument("--api", default="https://qcaihpc-simulation-api.onrender.com",
                    help="SOLANGE backend base URL (agent mode)")
    ap.add_argument("--poll", type=int, default=15, help="agent poll interval (seconds)")
    ap.add_argument("--token", default=None, help="Bearer JWT (agent mode); else log in")
    ap.add_argument("--email", default="guest@solange.bio", help="login email (agent mode)")
    ap.add_argument("--password", default="Solange2026", help="login password (agent mode)")
    args = ap.parse_args()

    Path(args.out).mkdir(parents=True, exist_ok=True)

    if args.agent:
        token = args.token
        if not token:
            print(f"[{_ts()}] [agent] logging in as {args.email} …")
            token = _supabase_login(args.email, args.password)
        run_agent(args.api, args.poll, token, args.out)
        return

    if not args.compound or not args.key:
        ap.error("--compound and --key are required for a direct run (or use --agent)")

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

    def _fmt(s):
        return f"{s:.1f}s" if s < 60 else f"{int(s//60)}m{s%60:04.1f}s"

    # Built-in real-time stage timing — no need for `time` or --verbose. Each stage
    # prints how long it took the moment it finishes; totals are archived (the
    # per-size runtime IS the classical-wall evidence for the dissertation).
    t_run = time.time()
    t0 = time.time()
    cas = run_casscf(args.compound, args.basis, args.ncas, args.nelecas, args.verbose)
    casscf_s = time.time() - t0
    print(f"RHF     E = {cas['e_rhf']:.8f} Ha  ({cas['n_ao']} AOs)")
    print(f"CASSCF  E = {cas['e_casscf']:.8f} Ha  (ecore {cas['ecore']:.8f})  "
          f"[done in {_fmt(casscf_s)}]")

    # The JW Pauli Hamiltonian is only needed to RUN the VQE. For an exact-reference
    # run (no --vqe) it is wasted work — the reference comes straight from PySCF
    # (cas['e_fci_active']) — and building ~O(ncas^4) terms dominates the runtime at
    # large sizes (e.g. ~30k terms / ~2 min at CAS(12,12)). So build it only when needed.
    jw_s = None
    if args.vqe:
        t0 = time.time()
        terms, e_active_exact, _ = build_jw_terms(cas["h1e"], cas["h2e"], args.ncas,
                                                  exact=(n_qubits <= 16))
        jw_s = time.time() - t0
        print(f"JW: {len(terms)} Pauli terms on {n_qubits} qubits"
              + (f" · exact active E = {e_active_exact:.8f} Ha" if e_active_exact is not None else "")
              + f"  [built in {_fmt(jw_s)}]")
    else:
        terms, e_active_exact = [], None
        print(f"Exact-reference run (no VQE) — skipping JW Hamiltonian build "
              f"(reference from PySCF FCI = {cas.get('e_fci_active'):.8f} Ha)")

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
    # Unique, timestamped stem so every run is archived (dissertation reproducibility)
    # instead of overwriting. e.g. ARID2_LOF_native_cas8-8_20260707T1830.
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    stem  = f"{args.key}_{args.side}_cas{args.nelecas}-{args.ncas}_{stamp}"
    jw_path = Path(args.out) / f"jw_{stem}.json"
    jw_path.write_text(json.dumps({args.key: {args.side: jw_entry}}, indent=2))
    print(f"WROTE {jw_path}")

    # ── Artifact 2: signed P1-P9 provenance ──
    compute_s = round(casscf_s + (jw_s or 0) + (vqe["elapsed_s"] if vqe else 0), 2)
    prov = build_provenance(args, cas, terms, e_active_exact, vqe, gpu_name, vram_mb,
                            elapsed_total=compute_s)
    prov_path = Path(args.out) / f"provenance_{stem}.json"
    prov_path.write_text(json.dumps(prov, indent=2, default=str))
    print(f"WROTE {prov_path}")
    print(f"P8 seal: {prov['p8_hash'][:16]}…  "
          f"(SOLANGE re-verifies this by recomputing the hash — that check needs no GPU)")

    # ── Master run log (append-only): one line per run, the local archive index ──
    exact_ref = e_active_exact if e_active_exact is not None else cas.get("e_fci_active")
    total_s = time.time() - t_run
    log_row = {
        "timestamp": stamp, "key": args.key, "side": args.side,
        "compound": args.compound, "basis": args.basis,
        "ncas": args.ncas, "nelecas": args.nelecas, "qubits": n_qubits,
        "e_casscf": cas["e_casscf"], "e_exact_active": exact_ref,
        "e_vqe_ha": (vqe["energy_ha"] if vqe else None),
        "delta_mha": ((vqe["energy_ha"] - exact_ref) * 1000
                      if (vqe and exact_ref is not None) else None),
        "device": (vqe["device"] if vqe else None),
        # stage timings — the per-size runtime is the classical-wall evidence
        "casscf_s": round(casscf_s, 2),
        "jw_build_s": (round(jw_s, 2) if jw_s is not None else None),
        "vqe_s": (vqe["elapsed_s"] if vqe else None),
        "total_s": round(total_s, 2),
        "n_paulis": len(terms), "p8_hash": prov["p8_hash"],
        "jw_file": jw_path.name, "provenance_file": prov_path.name,
    }
    log_path = Path(args.out) / "runs_log.jsonl"
    with open(log_path, "a") as lf:
        lf.write(json.dumps(log_row, default=str) + "\n")
    print(f"LOGGED → {log_path}  (append-only archive of every run)")
    print(f"TOTAL run time: {_fmt(total_s)}  "
          f"(CASSCF {_fmt(casscf_s)}"
          + (f" · JW {_fmt(jw_s)}" if jw_s is not None else "")
          + (f" · VQE {vqe['elapsed_s']}s" if vqe else "") + ")")
    print("=" * 68)

    if args.submit:
        try:
            import urllib.request
            url = args.submit.rstrip("/") + "/api/simulate/hpc/submit"
            body = json.dumps({"jw": jw_entry, "provenance": prov}, default=str).encode()
            req = urllib.request.Request(
                url, data=body, method="POST",
                headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                resp = json.loads(r.read().decode())
            print(f"SUBMITTED → {url}")
            print(f"  status={resp.get('status')}  seal_ok={resp.get('seal_ok')}  "
                  f"consistency_ok={resp.get('consistency_ok')}  db={resp.get('db_status')}  "
                  f"run_id={resp.get('run_id')}")
            # Printed explicitly (not just embedded above) so run_agent() can extract it
            # and so you can cross-reference this exact terminal run against the SOLANGE
            # HPC Ladder's P8 column by eye — same 8-char prefix, same hash.
            print(f"  P8 = {str(resp.get('recomputed_p8'))[:16]}… "
                  f"({'matched → verified' if resp.get('seal_ok') else 'MISMATCH'})")
            print("  → now visible in SOLANGE · phase=3A-HPC · provenance_source=HPC/Laguna")
        except Exception as e:
            print(f"SUBMIT FAILED ({e}). Files are saved locally — retry --submit or sync manually:")
            print(f"  {jw_path}\n  {prov_path}")
    else:
        print("NOTE: files are on Laguna only — NOT in SOLANGE. Add --submit to POST them")
        print("      (re-verified there), or sync manually:")
        print(f"  {jw_path}\n  {prov_path}")


if __name__ == "__main__":
    main()
