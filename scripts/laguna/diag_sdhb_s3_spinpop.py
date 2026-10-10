"""Which electronic configuration is the SDHB S3 high-spin SCF in?

The S=15/2 ROHF reference has several self-consistent solutions (different
Fe 3d / S 3p occupation patterns). Job 2350822's second-order SCF left the
earlier -7308.790 Ha solution and descended past -7308.96 Ha. <S^2> cannot
tell these apart -- every ROHF determinant with 15 alpha-only open shells
gives exactly 63.75 -- so this reads where the unpaired electrons actually
sit: meta-Lowdin spin populations per atom.

Expected for the oxidized [3Fe-4S]+ cluster (three high-spin Fe(III), d5):
about +4 to +4.5 on each Fe (covalency moves some spin onto the sulfurs),
summing to 15 over the whole cluster. Spin concentrated on S atoms with an
Fe well below that would mean a ligand-hole / charge-transfer configuration,
a different state from the one the active space was designed for.

Reads only the chkfile; safe to run while a job is still writing it.
"""
import argparse
from collections import defaultdict

import numpy as np
from pyscf import lib, scf
from pyscf.lo import orth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chkfile", default="sdhb_s3_scf.chk")
    a = ap.parse_args()

    mol = lib.chkfile.load_mol(a.chkfile)
    mol.verbose = 0
    data = scf.chkfile.load(a.chkfile, "scf")
    C, occ = data["mo_coeff"], data["mo_occ"]
    Ca, Cb = C[:, occ > 0], C[:, occ > 1]
    spin_dm = Ca @ Ca.T - Cb @ Cb.T

    # Meta-Lowdin orthogonalized AOs: populations that don't depend on the
    # basis-set overlap partitioning the way plain Mulliken does.
    S = mol.intor_symmetric("int1e_ovlp")
    U = orth.orth_ao(mol, "meta_lowdin", s=S)
    Uinv = np.linalg.solve(U, np.eye(U.shape[0]))
    pop_ao = np.diag(Uinv @ spin_dm @ Uinv.T)

    per_atom = np.zeros(mol.natm)
    for i, (ia, *_rest) in enumerate(mol.ao_labels(fmt=False)):
        per_atom[ia] += pop_ao[i]

    print(f"[scf] saved E={float(data['e_tot']):.8f}  total spin population="
          f"{per_atom.sum():+.4f} (must be +15)")
    by_elem = defaultdict(float)
    for ia in range(mol.natm):
        sym = mol.atom_symbol(ia)
        by_elem[sym] += per_atom[ia]
        if sym in ("Fe", "S") and abs(per_atom[ia]) > 0.05:
            print(f"[spin] atom {ia:3d} {sym:2s}  {per_atom[ia]:+.4f}")
    print("[spin] by element: " + "  ".join(f"{k}={v:+.3f}" for k, v in sorted(by_elem.items())))


if __name__ == "__main__":
    main()
