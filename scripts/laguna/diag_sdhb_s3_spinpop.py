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


def spin_populations(mol, C, occ):
    """Meta-Lowdin spin population per atom for an ROHF determinant."""
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
    return per_atom


def state_check(mol, per_atom, fe_window=(3.2, 4.4), fe_spread_max=0.5,
                off_cluster_abs_max=0.3, bridging_s_window=(0.2, 0.8)):
    """Is this the intended three-high-spin-Fe(III) [3Fe-4S]+ state? (ok, reasons)

    Thresholds from Claude Science's 2026-10-10 reply (replacing our first
    guess of "Fe >= 3.5, off-cluster <= 1.0"):
      * each Fe inside 3.2-4.4 -- covalency moves 15-30% of the formal five
        unpaired electrons onto the sulfides, so a clean high-spin Fe(III) in
        an Fe-S cluster sits near 3.5-4.2, not 5; a lower bound alone would
        also pass an unphysical ionic state;
      * spread between the three Fe <= 0.5 -- the three irons of [3Fe-4S]+
        are formally equivalent, so this is the sharpest single test
        (the ligand-radical state found in job 2350822 had a spread of 2.89).
        NOT valid for S2 [4Fe-4S]2+, which is mixed-valence by nature;
      * all Fe the same sign (ferromagnetic high-spin reference);
      * sum of |spin| on atoms that are neither Fe nor S <= 0.3;
      * each bridging sulfide (S within 2.6 A of two or more Fe) inside
        0.2-0.8 -- spin there is normal covalency, not an error.
    Found 2026-10-10: job 2350822's ROHF minimum had one Fe at +0.93 and about
    4 unpaired electrons on backbone N/C/H -- not the state the active space
    is built for.
    """
    reasons = []
    xyz = mol.atom_coords(unit="Angstrom")
    fe = [ia for ia in range(mol.natm) if mol.atom_symbol(ia) == "Fe"]
    for ia in fe:
        if not (fe_window[0] <= per_atom[ia] <= fe_window[1]):
            reasons.append(f"Fe atom {ia} spin {per_atom[ia]:+.3f} outside {fe_window}")
    if fe:
        spread = max(per_atom[fe]) - min(per_atom[fe])
        if spread > fe_spread_max:
            reasons.append(f"Fe spin spread {spread:.3f} > {fe_spread_max}")
        if len({np.sign(per_atom[ia]) for ia in fe}) > 1:
            reasons.append("Fe spins differ in sign")
    off = sum(abs(per_atom[ia]) for ia in range(mol.natm) if mol.atom_symbol(ia) not in ("Fe", "S"))
    if off > off_cluster_abs_max:
        reasons.append(f"sum |spin| on non-Fe/S atoms {off:.3f} > {off_cluster_abs_max}")
    for ia in range(mol.natm):
        if mol.atom_symbol(ia) != "S":
            continue
        n_fe = sum(1 for jf in fe if np.linalg.norm(xyz[ia] - xyz[jf]) < 2.6)
        if n_fe >= 2 and not (bridging_s_window[0] <= per_atom[ia] <= bridging_s_window[1]):
            reasons.append(f"bridging S atom {ia} spin {per_atom[ia]:+.3f} outside {bridging_s_window}")
    return (not reasons), reasons


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chkfile", default="sdhb_s3_scf.chk")
    a = ap.parse_args()

    mol = lib.chkfile.load_mol(a.chkfile)
    mol.verbose = 0
    data = scf.chkfile.load(a.chkfile, "scf")
    per_atom = spin_populations(mol, data["mo_coeff"], data["mo_occ"])

    print(f"[scf] saved E={float(data['e_tot']):.8f}  total spin population="
          f"{per_atom.sum():+.4f} (must be +15)")
    by_elem = defaultdict(float)
    for ia in range(mol.natm):
        sym = mol.atom_symbol(ia)
        by_elem[sym] += per_atom[ia]
        if sym in ("Fe", "S") and abs(per_atom[ia]) > 0.05:
            print(f"[spin] atom {ia:3d} {sym:2s}  {per_atom[ia]:+.4f}")
    print("[spin] by element: " + "  ".join(f"{k}={v:+.3f}" for k, v in sorted(by_elem.items())))
    ok, reasons = state_check(mol, per_atom)
    print("[state] intended 3x high-spin Fe(III) state: " + ("YES" if ok else "NO -- " + "; ".join(reasons)))


if __name__ == "__main__":
    main()
