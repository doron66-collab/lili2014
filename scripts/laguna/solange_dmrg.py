#!/usr/bin/env python3
"""
solange_dmrg.py — classical DMRG classifier for SOLANGE target A/B/C tiers.

Determines, from FIRST PRINCIPLES (not size), whether a target's active space is
classically solvable to chemical accuracy — the defensible basis for A/B/C:

  C  classical exact reaches chemical accuracy   (≤ ~18 active e⁻; CCSD(T)/FCI)
  B  DMRG (classical, approximate) reaches it     (converges at practical bond dim)
  A  even DMRG cannot                             (strong correlation / high
                                                    entanglement; quantum-necessary)

DMRG is a CLASSICAL algorithm (block2). It runs on CPU (largemem — big bond
dimensions need RAM) and can be GPU-accelerated; either way it is classical, NOT
quantum. The point of THIS script is exactly to find where classical (incl. DMRG)
stops reaching chemical accuracy — because that, not the 18e exact wall, is where
quantum is truly necessary.

Method: run DMRG at increasing bond dimension M. If the energy stops changing by
more than chemical accuracy (1 kcal/mol ≈ 1.6 mHa) at a practical M, DMRG has
delivered → the target is classically tractable (B). If it is still changing, or
the entanglement entropy is high (bond dim would blow up), DMRG has NOT delivered
→ quantum-necessary (A). The max bipartite entanglement entropy S_max is reported
as the physical reason.

USAGE (validate on a small case first):
  python solange_dmrg.py --compound acetamide --basis 6-31g --ncas 8 --nelecas 8 \
      --key ARID2_LOF --out ./out

HONEST SCOPE: this classifies whatever active space you give it. Classifying a
target's FULL functional site (e.g. ARID2 56-qubit site) rigorously requires
building that active space from the PDB structure first — a separate pipeline.
On the tiny model compounds it demonstrates the diagnostic, not a full-site verdict.
"""

import argparse
import hashlib
import json
import sys
import time
import uuid
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))          # for solange_hpc
sys.path.insert(0, str(_HERE.parents[2]))      # repo root, for generate_expansion_jw

CHEM_ACC_MHA = 1.6          # 1 kcal/mol
PRACTICAL_M  = 2000         # bond dimension beyond which DMRG is deemed impractical here
EXACT_WALL_E = 18           # active electrons up to which exact classical (FCI) is fine
S_HARD       = 1.5          # max bipartite entanglement entropy above which the DMRG
                            # bond dimension (~e^S) blows up at the full site -> DMRG-hard.
                            # Calibration: N2 equilibrium S_max=0.42 (easy) vs N2 stretched
                            # S_max=2.86 (strongly correlated). 1.5 sits between them.


def run_dmrg(h1e, h2e, ecore, ncas, nelecas, bond_dims, scratch="./tmp_dmrg",
             n_threads=4, max_minutes=None):
    """Run DMRG at increasing bond dimensions. Returns per-M energies + S_max.

    HPC-ticket-aware: prints live per-M timing (so `tail -f` shows real progress,
    letting you judge whether to Ctrl+C before a fixed-walltime allocation ends),
    and honors max_minutes — once the cumulative wall-clock budget is exceeded, it
    stops requesting new (larger) bond dimensions and returns whatever it has
    rather than being killed mid-sweep with nothing recorded. Reusing the same
    `scratch` directory across separate job submissions lets block2 resume the MPS
    from where a prior run left off instead of restarting from bond_dims[0].
    """
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=n_threads)
    drv.initialize_system(n_sites=ncas, n_elec=nelecas, spin=0)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=h2e, ecore=ecore, iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=min(bond_dims[0], 250), nroots=1)
    energies = []
    t_start = time.time()
    for M in bond_dims:
        elapsed_before = time.time() - t_start
        if max_minutes is not None and elapsed_before > max_minutes * 60 and energies:
            print(f"  [time budget] {elapsed_before/60:.1f}m elapsed > --max-minutes "
                  f"{max_minutes} — stopping before M={M}; returning {len(energies)} "
                  f"completed bond dim(s). Re-run with the same --scratch to resume.",
                  file=sys.stderr)
            break
        t0 = time.time()
        e = drv.dmrg(mpo, ket, n_sweeps=10, bond_dims=[M],
                     noises=[1e-5, 1e-6, 0], thrds=[1e-9] * 3, iprint=0)
        dt = time.time() - t0
        energies.append((M, float(e)))
        print(f"  DMRG M={M:5d}  E={float(e):.8f} Ha  [{dt:.1f}s this M, "
              f"{(time.time()-t_start)/60:.1f}m total]", file=sys.stderr)
    try:
        s_max = float(np.max(drv.get_bipartite_entanglement(ket)))
    except Exception:
        s_max = None
    return energies, s_max


def classify(active_electrons, energies, s_max):
    """Map DMRG behaviour to an A/B/C class with an explicit rationale.

    Two signals:
      • ΔE across the last two bond dims — direct evidence DMRG has (not) converged
        at a practical M. Only bites at large active spaces; small ones are exact.
      • S_max (max bipartite entanglement) — predicts the bond dimension the FULL
        site needs (M ~ e^S). This is the leading indicator, since a small model
        space is always DMRG-exact yet still reveals the correlation strength.
    """
    dE_final = abs(energies[-1][1] - energies[-2][1]) * 1000 if len(energies) >= 2 else None
    converged = (dE_final is not None) and (dE_final < CHEM_ACC_MHA) \
                and (energies[-1][0] <= PRACTICAL_M)
    strong = (s_max is not None) and (s_max > S_HARD)
    smx = "n/a" if s_max is None else f"{s_max:.2f}"

    if active_electrons <= EXACT_WALL_E:
        return "C", (f"{active_electrons}e ≤ {EXACT_WALL_E}e exact-classical wall — "
                     f"CCSD(T)/FCI reaches chemical accuracy; classical sufficient. "
                     f"(S_max={smx})")
    if converged and not strong:
        return "B", (f"DMRG reaches chemical accuracy at practical M={energies[-1][0]} "
                     f"(ΔE={dE_final:.2f} mHa) and entanglement is low (S_max={smx} < "
                     f"{S_HARD}) — classical (DMRG) delivers; quantum-advantaged, not necessary.")
    reasons = []
    if strong:
        reasons.append(f"S_max={smx} > {S_HARD} (strong correlation; DMRG bond dim ~e^S "
                       f"blows up at the full {active_electrons}e site)")
    if not converged:
        reasons.append(f"DMRG not at chemical accuracy by practical M={energies[-1][0]} "
                       f"(ΔE={'n/a' if dE_final is None else round(dE_final,2)} mHa)")
    return "A", "quantum-necessary — " + "; ".join(reasons) + "."


def integrals_from_geometry(xyz_path, basis, avas_aos):
    """Chemist-in-the-loop entry: given a QM-cluster geometry (xyz) and the target
    atomic orbitals, AVAS selects the active space automatically. Returns a dict
    shaped like run_casscf's output. The CLUSTER itself (which residues/atoms/metal,
    H-capping) is the chemist's input — that step is NOT auto-generated here."""
    from pyscf import gto, scf, ao2mo, fci
    from pyscf.mcscf import avas
    geom = Path(xyz_path).read_text()
    # accept a raw xyz (skip the 2 header lines if present)
    lines = [l for l in geom.splitlines() if l.strip()]
    if lines and lines[0].strip().isdigit():
        lines = lines[2:]
    mol = gto.M(atom="\n".join(lines), basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    ncas, nelec, mo = avas.avas(mf, [s.strip() for s in avas_aos.split(",")])
    from pyscf import mcscf
    mc = mcscf.CASSCF(mf, ncas, nelec)
    mc.fix_spin_(ss=0)
    e_casscf = mc.kernel(mo)[0]
    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), ncas)
    na = nelec // 2
    e_fci = fci.direct_spin1.FCI().kernel(h1e, h2e, ncas, (na, nelec - na), ecore=0.0)[0]
    return {"e_casscf": float(e_casscf), "ecore": float(ecore), "e_fci_active": float(e_fci),
            "h1e": h1e, "h2e": h2e, "ncas": int(ncas), "nelecas": int(nelec)}


# ── LEON seal (self-contained — mirrors backend/routes/leon.py's generic seal
# bit-for-bit, deliberately duplicated rather than imported: this script must run
# standalone on the Laguna cluster with no dependency on the SOLANGE backend
# package). LEON re-verifies this at ingestion; a mismatch is REJECTED, not stored.
def _seal_payload(record, exclude):
    return json.dumps({k: v for k, v in record.items() if k not in exclude},
                      sort_keys=True, default=str)


def _seal_hash(record, exclude):
    return hashlib.sha256(_seal_payload(record, exclude).encode()).hexdigest()


def _submit_dmrg(api, out):
    """POST the sealed DMRG classification to SOLANGE. Best-effort: prints the
    outcome but never raises — a submit failure must not discard the local JSON
    already written to --out."""
    import urllib.request
    try:
        url = api.rstrip("/") + "/api/simulate/hpc/dmrg/submit"
        body = json.dumps(out, default=str).encode()
        req = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            resp = json.loads(r.read().decode())
        print(f"SUBMITTED → {url}")
        print(f"  status={resp.get('status')}  seal_ok={resp.get('seal_ok')}  "
              f"db={resp.get('db_status')}  run_id={resp.get('run_id')}")
    except Exception as e:
        print(f"  SUBMIT FAILED (result is still safe locally in --out): {e}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="SOLANGE DMRG A/B/C classifier (Laguna largemem).")
    ap.add_argument("--compound", help="model compound (key in GEOM) — demo mode")
    ap.add_argument("--geometry", help="path to a QM-cluster .xyz (real functional site)")
    ap.add_argument("--avas", help="comma-separated AVAS target AOs, e.g. 'Zn 3d,S 3p'")
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ncas", type=int)
    ap.add_argument("--nelecas", type=int)
    ap.add_argument("--key", required=True, help="target key, e.g. ARID2_LOF")
    ap.add_argument("--bond-dims", default="250,500,1000,2000",
                    help="comma-separated increasing bond dimensions")
    ap.add_argument("--out", default="./out")
    ap.add_argument("--verbose", type=int, default=0)
    ap.add_argument("--threads", type=int, default=4,
                    help="CPU threads for block2 (was hard-coded to 4; raise this to use "
                         "more of a largemem node's cores and finish faster)")
    ap.add_argument("--scratch", default="./tmp_dmrg",
                    help="block2 scratch dir. Reuse the SAME path across separate HPC-ticket "
                         "submissions to resume an interrupted MPS instead of restarting from "
                         "bond_dims[0] each time.")
    ap.add_argument("--max-minutes", type=float, default=None,
                    help="stop requesting larger bond dimensions once this wall-clock budget "
                         "is exceeded, and return/save whatever completed so far — for "
                         "fixed-walltime HPC allocations (e.g. a 2-hour ticket) where getting "
                         "killed mid-sweep would lose everything.")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the sealed classification to SOLANGE so it's LEON-notarized and "
                         "stored immediately — safe even if this Laguna session later becomes "
                         "unreachable. Bare flag uses the default SOLANGE API URL; pass a value "
                         "to override.")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    bond_dims = [int(x) for x in args.bond_dims.split(",")]

    print("=" * 68)
    if args.geometry:
        if not args.avas:
            ap.error("--geometry requires --avas (target AOs for active-space selection)")
        print(f"SOLANGE DMRG classifier · {args.key} · geometry={args.geometry} "
              f"· AVAS[{args.avas}]/{args.basis}")
        cas = integrals_from_geometry(args.geometry, args.basis, args.avas)
        args.ncas, args.nelecas = cas["ncas"], cas["nelecas"]
        print(f"AVAS selected active space: CAS({args.nelecas},{args.ncas})")
    else:
        if not (args.compound and args.ncas and args.nelecas):
            ap.error("provide either --geometry+--avas, or --compound+--ncas+--nelecas")
        from solange_hpc import run_casscf
        print(f"SOLANGE DMRG classifier · {args.key} · {args.compound}/{args.basis} "
              f"· CAS({args.nelecas},{args.ncas})")
        cas = run_casscf(args.compound, args.basis, args.ncas, args.nelecas, args.verbose)
    print(f"CASSCF E = {cas['e_casscf']:.8f} Ha")

    t0 = time.time()
    energies, s_max = run_dmrg(cas["h1e"], cas["h2e"], cas["ecore"],
                               args.ncas, args.nelecas, bond_dims,
                               scratch=args.scratch, n_threads=args.threads,
                               max_minutes=args.max_minutes)
    # (per-M timing is already printed live inside run_dmrg, as each M finishes —
    # so a `tail -f` on a background run shows real progress, not a single dump at exit.)
    print(f"max bipartite entanglement S_max = {s_max}")
    time_budget_hit = len(energies) < len(bond_dims)
    if time_budget_hit:
        print(f"NOTE: stopped early at {len(energies)}/{len(bond_dims)} bond dimensions "
              f"(--max-minutes {args.max_minutes}). The class below is PROVISIONAL — "
              f"re-run with --scratch {args.scratch} to resume toward the full bond-dim list.")

    cls, rationale = classify(args.nelecas, energies, s_max)
    print("-" * 68)
    print(f"CLASS {cls}{' (PROVISIONAL — time budget hit)' if time_budget_hit else ''}")
    print(f"  {rationale}")
    print(f"elapsed {round(time.time()-t0,1)}s")

    out = {
        "id": str(uuid.uuid4()),
        "key": args.key, "compound": args.compound, "basis": args.basis,
        "ncas": args.ncas, "nelecas": args.nelecas,
        "e_casscf": cas["e_casscf"],
        "dmrg_energies": energies, "s_max": s_max,
        "bqp_class": cls, "class_rationale": rationale,
        "time_budget_hit": time_budget_hit, "bond_dims_requested": bond_dims,
        "method": "DMRG (block2, classical) convergence + entanglement diagnostic",
        "provenance_source": "HPC/Laguna (DMRG classifier)",
    }
    # Seal at source (LEON re-verifies at ingestion — a mismatch is rejected, not
    # trusted). dmrg_seal_payload is stored verbatim so re-verification later is
    # exact-string, not float-reconstruction (the same robustness fix the P8 seal
    # got after floats/timestamps proved to reformat across a DB round-trip).
    out["dmrg_seal_payload"] = _seal_payload(out, exclude={"dmrg_hash", "dmrg_seal_payload"})
    out["dmrg_hash"] = hashlib.sha256(out["dmrg_seal_payload"].encode()).hexdigest()

    p = Path(args.out) / f"dmrg_class_{args.key}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"WROTE {p}")
    if args.submit:
        _submit_dmrg(args.submit, out)
    print("=" * 68)
    print("NOTE: DMRG is classical (CPU/largemem or GPU-accelerated), NOT quantum.")
    print("Full-site classification needs the target's active space built from PDB first.")


if __name__ == "__main__":
    main()
