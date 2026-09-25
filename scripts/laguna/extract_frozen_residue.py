#!/usr/bin/env python3
"""
extract_frozen_residue.py — extract a SINGLE residue, at its crystallographic
(frozen) geometry, capped exactly the way build_qm_cluster.py caps a cluster
boundary — as a standalone .xyz "model compound" for isodesmic-correction
work.

WHY THIS EXISTS
Found 2026-09-23 on TP53 Y220C: the isodesmic correction compares a QM
cluster (frozen at its crystallographic geometry) against small reference
compounds (p-cresol, methanethiol, toluene, ...) computed at their OWN
relaxed, energy-minimized geometries. That mismatch is the leading
explanation for why the corrected result was directionally right but ~13x
too large in magnitude. This script tests that hypothesis directly: extract
just the mutation-site residue itself, capped the same way the cluster's own
boundary residues are, and compute ITS energy at the FROZEN geometry — for
comparison against the existing relaxed reference-compound energy for the
same residue type.

This deliberately reuses build_qm_cluster.py's own selection and capping
functions rather than reimplementing them, so the capping convention is
identical to what the full cluster already uses (same CAP_BOND_LENGTH, same
direction-along-original-bond placement) — not a second, possibly
inconsistent scheme.

USAGE:
  python extract_frozen_residue.py --pdb 2OCJ_native_H.pdb --chain A --resi 220 \\
      --out y220_native_frozen_residue.xyz
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_qm_cluster import (
    parse_structure, cap_dangling_bonds, write_xyz, suggest_avas,
)


def select_single_residue(atoms, chain, resi):
    """Select ONLY the named residue's own atoms — no radius, no neighbors.
    Returns the same {(chain, resseq, icode): [atom, ...]} shape
    select_residues() does, so cap_dangling_bonds()/write_xyz() work
    unmodified. Both backbone neighbors (resi-1, resi+1) are guaranteed
    absent from this selection by construction, so cap_dangling_bonds()
    will cap both the N-side and C-side backbone bonds — turning this into
    a standalone, doubly end-capped fragment, same as any cluster boundary
    residue already is."""
    hits = [a for a in atoms if a["chain"] == chain and a["resseq"] == resi]
    if not hits:
        raise SystemExit(f"no atoms found at chain={chain} resi={resi}")
    icode = hits[0]["icode"]
    return {(chain, resi, icode): hits}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="local structure file, .pdb or .cif")
    ap.add_argument("--chain", required=True)
    ap.add_argument("--resi", type=int, required=True)
    ap.add_argument("--out", required=True, help="output .xyz path")
    args = ap.parse_args()

    atoms = parse_structure(args.pdb)
    selected = select_single_residue(atoms, args.chain, args.resi)
    resname = selected[list(selected.keys())[0]][0]["resname"]
    caps = cap_dangling_bonds(selected)
    n_heavy = sum(1 for ats in selected.values() for a in ats)
    comment = (f"Frozen-geometry residue: {args.pdb} chain={args.chain} resi={args.resi} "
               f"({resname}) heavy_atoms={n_heavy} caps={len(caps)} — "
               f"crystallographic geometry, NOT relaxed")
    rows = write_xyz(selected, caps, args.out, include_hetero=True, comment=comment)

    print(f"extracted residue {resname} {args.chain}{args.resi}: "
          f"{n_heavy} heavy atoms, {len(caps)} capping H (both backbone ends, "
          f"since no chain-sequence neighbor was selected)")
    print(f"wrote {args.out}")
    print(f"suggested --avas: \"{suggest_avas(rows)}\"")

    _Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}
    n_elec = sum(_Z.get((r[0] or "").strip().upper(), 0) for r in rows)
    print(f"neutral-atom electron count: {n_elec} ({'even' if n_elec % 2 == 0 else 'ODD'})")
    print("NOTE: this is a bare N-H/C-H capped fragment (this project's own capping "
          "convention), not the standard ACE/NME-capped amino-acid model used in most "
          "QM/MM literature -- consistent with the cluster it's extracted from, which "
          "is what matters for this comparison, but worth stating rather than assuming.")


if __name__ == "__main__":
    main()
