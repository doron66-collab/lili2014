#!/usr/bin/env python3
"""
active_space_population.py -- per-residue Lowdin population of a selected active
space, and an explicit residue-275 orbital count.

Claude Science's second sign-off condition on CAS(36,35)/(36,34) (2026-10-01):
AVAS/MP2 pick orbitals by a site PROJECTOR SCORE and a correlation DEVIATION,
neither of which says, in human terms, "does this orbital actually sit on the
mutated residue." The site_score check (already in avas_mp2_select.py's
output JSON, reproducing plain AVAS weights) rules out gross drift away from
the AO criterion as a whole, but not whether a SPECIFIC residue's orbitals
(e.g. PHE275's ring pi/pi* system) ended up inside the chosen cut.

This needs the atom-index -> residue mapping build_qm_cluster.py's --atom-map-out
now writes (the bare .xyz the rest of the pipeline uses has no residue labels
at all -- that information existed at cluster-build time and was being thrown
away). If a cluster's .xyz predates --atom-map-out, this script cannot run
until the SAME build_qm_cluster.py command (same --target/--chain/--resi/
--radius) is re-run with --atom-map-out added; atom order is deterministic
from the same PDB + selection, so this is a cheap re-derivation, not a QM run.

Usage:
    python active_space_population.py --json cluster_avas_mp2.json \
        --mo cluster_avas_mp2_mo.npy --atom-map cluster_atommap.json \
        --residue 275 --chain A
"""
import argparse
import json
import sys

import numpy


def lowdin_population_per_atom(mf, mo_columns):
    """Lowdin population of each AO on each atom, summed over the given mo
    columns (equal-weighted -- this is a composition check, not a density).
    Returns an (natm,) array.
    """
    import scipy.linalg
    mol = mf.mol
    s = mol.intor_symmetric('int1e_ovlp')
    # Need S^{+1/2} here (psi = chi C = chi' S^{1/2} C, so S^{1/2} C are the MO
    # coefficients in the Lowdin-orthogonalized AO basis chi' = chi S^{-1/2}).
    # NOT lo.orth.lowdin(s), which returns S^{-1/2} (the AO-orthogonalization
    # transform itself, a different matrix) -- using it here silently breaks
    # normalization: caught by the self-check below before this was trusted.
    e, v = scipy.linalg.eigh(s)
    s_half = (v * numpy.sqrt(e)).dot(v.T)
    c_orth = s_half.dot(mo_columns)          # (nao, n_selected)
    col_norms = numpy.sum(c_orth ** 2, axis=0)
    if numpy.max(numpy.abs(col_norms - 1.0)) > 1e-6:
        raise RuntimeError(f"[pop] Lowdin population normalization failed: column norms "
                           f"{col_norms} should all be 1.0 (input orbitals must be "
                           f"orthonormal under S). Refusing to report a population analysis "
                           f"that doesn't sum correctly.")
    weight = numpy.sum(c_orth ** 2, axis=1)  # (nao,) summed over selected orbitals
    ao_atom = [lbl[0] for lbl in mol.ao_labels(fmt=None)]  # ao_labels(fmt=None)[i] = (atom_idx, elem, nlm, pol)
    natm = mol.natm
    per_atom = numpy.zeros(natm)
    for w, iatm in zip(weight, ao_atom):
        per_atom[int(iatm)] += w
    return per_atom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, help="avas_mp2_select.py output JSON")
    ap.add_argument("--mo", required=True, help="matching _mo.npy")
    ap.add_argument("--atom-map", required=True,
                    help="sidecar JSON from build_qm_cluster.py --atom-map-out, for the SAME cluster")
    ap.add_argument("--residue", type=int, required=True, help="residue number to report on, e.g. 275")
    ap.add_argument("--chain", default=None, help="chain ID, if the atom map has more than one")
    a = ap.parse_args()

    meta = json.load(open(a.json))
    amap = json.load(open(a.atom_map))
    atom_residues = amap["atom_residues"]

    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
    from run_gate2_avas import build_mf
    mf = build_mf(meta["xyz"], meta["charge"], meta["spin"], basis=meta["basis"],
                  chkfile=meta.get("chkfile"))

    if mf.mol.natm != len(atom_residues):
        sys.exit(f"[pop] REFUSING: {a.json}'s xyz has {mf.mol.natm} atoms but {a.atom_map} "
                 f"records {len(atom_residues)} -- these do not describe the same cluster. "
                 f"The atom map must come from the SAME build_qm_cluster.py run that produced "
                 f"this xyz (same --target/--chain/--resi/--radius), not a different one.")

    lo_col = meta["active_start_col"]
    hi_col = lo_col + meta["n_occ"] + meta["n_vir"]
    mo = numpy.load(a.mo)
    active_mo = mo[:, lo_col:hi_col]

    per_atom = lowdin_population_per_atom(mf, active_mo)
    per_residue = {}
    for iatm, (pop, res) in enumerate(zip(per_atom, atom_residues)):
        key = (res["chain"], res["resseq"], res["resname"])
        per_residue[key] = per_residue.get(key, 0.0) + pop

    print(f"[pop] Lowdin population of the {meta['n_occ']+meta['n_vir']}-orbital active space, "
          f"by residue (total should sum close to {meta['n_occ']+meta['n_vir']}):")
    total = 0.0
    for key, pop in sorted(per_residue.items(), key=lambda kv: -kv[1]):
        chain, resseq, resname = key
        total += pop
        if pop > 0.01:
            print(f"    {resname}{resseq} (chain {chain}): {pop:.4f}")
    print(f"[pop] sum over all residues: {total:.4f} (vs {meta['n_occ']+meta['n_vir']} orbitals selected)")

    target_pop = sum(pop for (chain, resseq, resname), pop in per_residue.items()
                      if resseq == a.residue and (a.chain is None or chain == a.chain))
    print(f"\n[pop] residue {a.residue}" + (f" chain {a.chain}" if a.chain else "") +
          f": {target_pop:.4f} population-equivalent orbitals inside the "
          f"{meta['n_occ']+meta['n_vir']}-orbital active space")
    if target_pop < 0.5:
        print(f"[pop] WARNING: well under 1 orbital's worth of population on residue {a.residue} -- "
              f"if this is the mutated residue, the active space may not represent the mutation "
              f"at all. This is a flag to inspect, not an automatic refusal: the active space is "
              f"necessarily delocalized across neighbors too (AVAS_threshold={meta.get('site_threshold')}).")


if __name__ == "__main__":
    main()
