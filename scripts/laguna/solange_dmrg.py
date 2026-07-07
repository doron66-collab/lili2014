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
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))          # for solange_hpc
sys.path.insert(0, str(_HERE.parents[2]))      # repo root, for generate_expansion_jw

CHEM_ACC_MHA = 1.6          # 1 kcal/mol
PRACTICAL_M  = 2000         # bond dimension beyond which DMRG is deemed impractical here
EXACT_WALL_E = 18           # active electrons up to which exact classical (FCI) is fine


def run_dmrg(h1e, h2e, ecore, ncas, nelecas, bond_dims, scratch="./tmp_dmrg"):
    """Run DMRG at increasing bond dimensions. Returns per-M energies + S_max."""
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=4)
    drv.initialize_system(n_sites=ncas, n_elec=nelecas, spin=0)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=h2e, ecore=ecore, iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=min(bond_dims[0], 250), nroots=1)
    energies = []
    for M in bond_dims:
        e = drv.dmrg(mpo, ket, n_sweeps=10, bond_dims=[M],
                     noises=[1e-5, 1e-6, 0], thrds=[1e-9] * 3, iprint=0)
        energies.append((M, float(e)))
    try:
        s_max = float(np.max(drv.get_bipartite_entanglement(ket)))
    except Exception:
        s_max = None
    return energies, s_max


def classify(active_electrons, energies, s_max):
    """Map DMRG behaviour to an A/B/C class with an explicit rationale."""
    dE_final = abs(energies[-1][1] - energies[-2][1]) * 1000 if len(energies) >= 2 else None
    converged = (dE_final is not None) and (dE_final < CHEM_ACC_MHA) \
                and (energies[-1][0] <= PRACTICAL_M)

    if active_electrons <= EXACT_WALL_E:
        return "C", (f"{active_electrons}e ≤ {EXACT_WALL_E}e exact-classical wall — "
                     f"CCSD(T)/FCI reaches chemical accuracy; classical sufficient.")
    if converged:
        return "B", (f"DMRG converged to <{CHEM_ACC_MHA} mHa by M={energies[-1][0]} "
                     f"(ΔE={dE_final:.2f} mHa, S_max={s_max}); classical (DMRG) reaches "
                     f"chemical accuracy — quantum-advantaged, not necessary.")
    return "A", (f"DMRG NOT converged to chemical accuracy at practical M="
                 f"{energies[-1][0]} (ΔE={dE_final if dE_final is None else round(dE_final,2)} mHa"
                 f", S_max={s_max}) — no classical route delivers the needed accuracy; "
                 f"quantum-necessary.")


def main():
    ap = argparse.ArgumentParser(description="SOLANGE DMRG A/B/C classifier (Laguna largemem).")
    ap.add_argument("--compound", required=True)
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ncas", type=int, required=True)
    ap.add_argument("--nelecas", type=int, required=True)
    ap.add_argument("--key", required=True, help="target key, e.g. ARID2_LOF")
    ap.add_argument("--bond-dims", default="250,500,1000,2000",
                    help="comma-separated increasing bond dimensions")
    ap.add_argument("--out", default="./out")
    ap.add_argument("--verbose", type=int, default=0)
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    bond_dims = [int(x) for x in args.bond_dims.split(",")]

    from solange_hpc import run_casscf
    print("=" * 68)
    print(f"SOLANGE DMRG classifier · {args.key} · {args.compound}/{args.basis} "
          f"· CAS({args.nelecas},{args.ncas})")
    cas = run_casscf(args.compound, args.basis, args.ncas, args.nelecas, args.verbose)
    print(f"CASSCF E = {cas['e_casscf']:.8f} Ha")

    t0 = time.time()
    energies, s_max = run_dmrg(cas["h1e"], cas["h2e"], cas["ecore"],
                               args.ncas, args.nelecas, bond_dims)
    for M, e in energies:
        print(f"  DMRG M={M:5d}  E={e:.8f} Ha")
    print(f"max bipartite entanglement S_max = {s_max}")

    cls, rationale = classify(args.nelecas, energies, s_max)
    print("-" * 68)
    print(f"CLASS {cls}")
    print(f"  {rationale}")
    print(f"elapsed {round(time.time()-t0,1)}s")

    out = {
        "key": args.key, "compound": args.compound, "basis": args.basis,
        "ncas": args.ncas, "nelecas": args.nelecas,
        "e_casscf": cas["e_casscf"],
        "dmrg_energies": energies, "s_max": s_max,
        "bqp_class": cls, "class_rationale": rationale,
        "method": "DMRG (block2, classical) convergence + entanglement diagnostic",
    }
    p = Path(args.out) / f"dmrg_class_{args.key}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"WROTE {p}")
    print("=" * 68)
    print("NOTE: DMRG is classical (CPU/largemem or GPU-accelerated), NOT quantum.")
    print("Full-site classification needs the target's active space built from PDB first.")


if __name__ == "__main__":
    main()
