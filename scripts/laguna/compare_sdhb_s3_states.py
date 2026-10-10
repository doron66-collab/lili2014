"""Pre-registered test: do the restored N-H...S hydrogen bonds pull spin off the sulfides?

Claude Science's prediction (2026-10-10), fixed BEFORE job 2350825 finished.
Going from the 43-atom model (job 2350823) to the backbone-extended model
(job 2350825), both with ROKS/BP86 + ddCOSMO(eps=4):

    |delta spin(S_mu3)|  <  |delta spin| of every mu2 sulfide
    and   mean spin(Fe) rises above 3.457

The single mu3 sulfide sits buried in the cuboid and gains NO new N-H donor,
so it is an internal control inside the same calculation. The mu2 sulfides
gain 1-2 donors each. This is an ORDER prediction, not a value prediction:
it needs no calibration and it can fail. If the mu3 sulfide moves as much as
the others, the hydrogen-bond diagnosis is wrong and the cause is the
functional or the solvent model.

Diagnostic only -- it does not gate anything (Science: the extended model is
the one to run either way).
"""
import argparse

import numpy as np
from pyscf import lib, scf

from diag_sdhb_s3_spinpop import spin_populations, sulfur_roles


def load(chk):
    mol = lib.chkfile.load_mol(chk)
    mol.verbose = 0
    d = scf.chkfile.load(chk, "scf")
    return mol, spin_populations(mol, d["mo_coeff"], d["mo_occ"])


def cluster_atoms(mol):
    """Fe atoms and Fe-bound S atoms: index -> (symbol, role, xyz in Angstrom)."""
    xyz = mol.atom_coords(unit="Angstrom")
    roles = sulfur_roles(mol)
    out = {}
    for ia in range(mol.natm):
        if mol.atom_symbol(ia) == "Fe":
            out[ia] = ("Fe", "Fe", xyz[ia])
        elif ia in roles:
            out[ia] = ("S", roles[ia], xyz[ia])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--small", default="sdhb_s3_scf_bp86_eps4.chk")
    ap.add_argument("--extended", default="sdhb_s3_scf_bb_mix_bp86_eps4.chk")
    a = ap.parse_args()
    mol_a, pop_a = load(a.small)
    mol_b, pop_b = load(a.extended)
    ca, cb = cluster_atoms(mol_a), cluster_atoms(mol_b)

    rows = []
    for ia, (sym, role, xa) in ca.items():
        match = [ib for ib, (s2, _, xb) in cb.items() if s2 == sym and np.linalg.norm(xa - xb) < 0.01]
        if len(match) != 1:
            raise SystemExit(f"could not match {sym} atom {ia} between the two models")
        ib = match[0]
        rows.append((role, ia, ib, pop_a[ia], pop_b[ib], pop_b[ib] - pop_a[ia]))

    print(f"{'role':9s} {'idx small':>9s} {'idx ext':>8s} {'small':>8s} {'extended':>9s} {'delta':>8s}")
    for role, ia, ib, sa, sb, d in sorted(rows, key=lambda r: (r[0], r[1])):
        print(f"{role:9s} {ia:9d} {ib:8d} {sa:+8.3f} {sb:+9.3f} {d:+8.3f}")

    mu3 = [abs(r[5]) for r in rows if r[0] == "mu3"]
    mu2 = [abs(r[5]) for r in rows if r[0] == "mu2"]
    fe_b = [r[4] for r in rows if r[0] == "Fe"]
    order_ok = bool(mu3) and bool(mu2) and max(mu3) < min(mu2)
    fe_ok = np.mean(fe_b) > 3.457
    print(f"\n[prediction] |delta mu3| = {max(mu3) if mu3 else float('nan'):.3f} < min |delta mu2| = "
          f"{min(mu2) if mu2 else float('nan'):.3f}: {'HELD' if order_ok else 'FAILED'}")
    print(f"[prediction] mean Fe spin (extended) = {np.mean(fe_b):.3f} > 3.457: {'HELD' if fe_ok else 'FAILED'}")
    # The pre-registered prediction is a conjunction of the two parts; report it
    # as such, and say nothing about the cause. An earlier version of this line
    # attributed a FAILED verdict to "the functional or the solvent model" -- an
    # interpretation the prediction itself does not license when one part holds
    # and the other fails (job 2350835: order held, Fe part failed).
    print("[prediction] overall (both parts required): " + ("HELD" if order_ok and fe_ok else "FAILED")
          + f"  [order part: {'held' if order_ok else 'failed'}; Fe part: {'held' if fe_ok else 'failed'}]")
    off_note = ("spin removed from the sulfides did not return to the irons" if not fe_ok else
                "spin removed from the sulfides returned to the irons")
    print(f"[observation] {off_note}; compare the off-cluster sums printed by "
          f"diag_sdhb_s3_spinpop.py for both chkfiles to see where it went.")

if __name__ == "__main__":
    main()
