"""Map the F3S303 atom names (FE1/FE3/FE4, S1-S4) of 8GS8 to atom indices and spins.

Science's dose-response check (third reply, 2026-10-10): the shortest new
N-H...S donor is S4 -> 2.94 A (ASN248), S2 -> 3.14 A (MET247), S1 -> 3.69 A
(HIS244); S3 is the buried mu3 sulfide. If the spin lost by each mu2 sulfide
is ordered S4 > S2 > S1, the loss is monotonic in donor distance -- a second,
independent confirmation of the hydrogen-bond diagnosis. Atoms are matched
by coordinates against the deposited mmCIF, never assumed from list order.
"""
import argparse
import os
import sys

import numpy as np
from pyscf import lib, scf

sys.path.insert(0, os.getcwd())
from diag_sdhb_s3_spinpop import spin_populations  # noqa: E402
from run_gate2_sdhb_s3 import read_f3s_from_cif, PDB_ID, CHAIN  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chkfiles", nargs="+")
    a = ap.parse_args()
    f3s = read_f3s_from_cif(PDB_ID, CHAIN, "303")
    for chk in a.chkfiles:
        mol = lib.chkfile.load_mol(chk)
        mol.verbose = 0
        d = scf.chkfile.load(chk, "scf")
        pop = spin_populations(mol, d["mo_coeff"], d["mo_occ"])
        xyz = mol.atom_coords(unit="Angstrom")
        print(f"=== {chk}")
        for atm in f3s:
            hits = [i for i in range(mol.natm) if np.linalg.norm(xyz[i] - np.array(atm["xyz"])) < 0.01]
            if len(hits) != 1:
                print(f"  {atm['name']:4s} no unique coordinate match ({len(hits)})")
                continue
            print(f"  {atm['name']:4s} atom {hits[0]:4d}  spin {pop[hits[0]]:+.4f}")


if __name__ == "__main__":
    main()
