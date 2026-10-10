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


def sulfur_roles(mol, cutoff=2.6):
    """{atom index: 'mu3' | 'mu2' | 'thiolate'} for every S bonded to an Fe."""
    xyz = mol.atom_coords(unit="Angstrom")
    fe = [ia for ia in range(mol.natm) if mol.atom_symbol(ia) == "Fe"]
    roles = {}
    for ia in range(mol.natm):
        if mol.atom_symbol(ia) != "S":
            continue
        n = sum(1 for jf in fe if np.linalg.norm(xyz[ia] - xyz[jf]) < cutoff)
        if n:
            roles[ia] = {1: "thiolate", 2: "mu2"}.get(n, "mu3")
    return roles


def _off_cluster(mol):
    return [ia for ia in range(mol.natm) if mol.atom_symbol(ia) not in ("Fe", "S")]


def state_check(mol, per_atom, fe_spread_max=0.5, offcluster_atom_to_fe_max=0.10,
                ligand_to_fe_max=0.35, fe_fraction_window=(0.60, 0.85), two_s=15):
    """Is this the intended three-high-spin-Fe(III) [3Fe-4S]+ state? (ok, reasons)

    Gate as revised by Claude Science, 2026-10-10 (second reply). Their first
    version had an absolute 0.2-0.8 bound on each bridging sulfide; Science
    retracted it as wrong in FORM, not only in value. It was absolute (a bound
    on the partition scheme, not on the physics), it applied one bound to mu2
    and mu3 sulfides that receive spin from two vs three irons, and it was
    blind to the functional. It also double-counted: spin that leaves the
    irons must land on the ligands. Job 2350823 (BP86 + ddCOSMO) showed the
    failure mode: three equivalent Fe (+3.470/+3.454/+3.448, spread 0.022),
    refused on two sulfides at +0.821/+0.845.

    Kept:      Fe spread <= 0.5; all Fe same sign.
    Replaced (Science, third reply, 2026-10-10): off-cluster sum|s| <= 0.3 was
               EXTENSIVE in atom count -- the same defect in form as the
               withdrawn per-sulfur bound. The extended model (123 atoms)
               was refused at 0.451 while its per-atom off-cluster spin had
               HALVED (0.00767 -> 0.00399). Now: max(single non-Fe/non-S
               atom) / min(Fe) <= 0.10 -- intensive; asks only whether any
               backbone atom is a radical centre. The sum and the per-atom
               mean are reported, not gated.
    Replaced:  max(ligand spin) / min(Fe spin) <= 0.35 -- scale-free: 0.245
               for the BP86 state, 1.118 for the ROHF ligand-radical state
               (a ligand outranked an iron there);
               sum spin(Fe) / 2S inside 0.60-0.85 -- 0.691 for the BP86 state.
    Reported, not gated: every Fe-bound sulfur by role (mu3 / mu2 / thiolate).
    Asserted:  total spin population == 2S.
    The thresholds are valid for S3 [3Fe-4S]+ (and S1); NOT for S2
    [4Fe-4S]2+, which is mixed-valence by nature.
    """
    reasons = []
    total = float(np.sum(per_atom))
    if abs(total - two_s) > 1e-3:
        reasons.append(f"total spin population {total:.4f} != 2S = {two_s}")
    fe = [ia for ia in range(mol.natm) if mol.atom_symbol(ia) == "Fe"]
    if not fe:
        return False, ["no Fe atoms"]
    fe_spins = per_atom[fe]
    spread = float(fe_spins.max() - fe_spins.min())
    if spread > fe_spread_max:
        reasons.append(f"Fe spin spread {spread:.3f} > {fe_spread_max}")
    if len({np.sign(x) for x in fe_spins}) > 1:
        reasons.append("Fe spins differ in sign")
    off_idx = _off_cluster(mol)
    off_max = max((abs(per_atom[ia]) for ia in off_idx), default=0.0)
    if off_max / float(fe_spins.min()) > offcluster_atom_to_fe_max:
        reasons.append(f"max(single non-Fe/S atom)/min(Fe) {off_max / float(fe_spins.min()):.3f} "
                       f"> {offcluster_atom_to_fe_max}")
    lig_max = max((abs(per_atom[ia]) for ia in range(mol.natm) if mol.atom_symbol(ia) != "Fe"), default=0.0)
    ratio = lig_max / float(fe_spins.min())
    if ratio > ligand_to_fe_max:
        reasons.append(f"max(ligand spin)/min(Fe spin) {ratio:.3f} > {ligand_to_fe_max}")
    frac = float(fe_spins.sum()) / two_s
    if not (fe_fraction_window[0] <= frac <= fe_fraction_window[1]):
        reasons.append(f"sum(Fe spin)/2S {frac:.3f} outside {fe_fraction_window}")
    return (not reasons), reasons


def state_report(mol, per_atom, two_s=15):
    """Lines to print: every Fe, every Fe-bound S by role, the scale-free numbers."""
    fe = [ia for ia in range(mol.natm) if mol.atom_symbol(ia) == "Fe"]
    lines = [f"Fe atom {ia} spin {per_atom[ia]:+.3f}" for ia in fe]
    # Reported as sulfide / thiolate: the mu2/mu3 split was withdrawn in the
    # pre-registration (§6) -- the measured mu3 value sits inside the mu2 range.
    # The bridging count is kept in brackets because compare_sdhb_s3_states.py's
    # pre-registered prediction is stated in terms of the mu3 sulfide.
    for ia, role in sorted(sulfur_roles(mol).items(), key=lambda kv: (kv[1] == "thiolate", kv[0])):
        kind = "thiolate" if role == "thiolate" else f"sulfide [{role}]"
        lines.append(f"S atom {ia} ({kind}) spin {per_atom[ia]:+.3f}")
    off_idx = _off_cluster(mol)
    off = sum(abs(per_atom[ia]) for ia in off_idx)
    lig_max = max((abs(per_atom[ia]) for ia in range(mol.natm) if mol.atom_symbol(ia) != "Fe"), default=0.0)
    if off_idx:
        imax = max(off_idx, key=lambda ia: abs(per_atom[ia]))
        lines += [f"max single non-Fe/S atom: {mol.atom_symbol(imax)} atom {imax} {per_atom[imax]:+.4f} "
                  f"(/min Fe = {abs(per_atom[imax]) / min(per_atom[fe]):.4f}, gate <= 0.10)",
                  f"sum |spin| on non-Fe/S atoms {off:.3f} over {len(off_idx)} atoms "
                  f"(mean {off / len(off_idx):.5f}; reported, not gated)"]
    lines += [
              f"max(ligand)/min(Fe) {lig_max / min(per_atom[fe]):.3f}",
              f"sum(Fe)/2S {sum(per_atom[fe]) / two_s:.3f}",
              f"total spin population {float(np.sum(per_atom)):.4f} (2S = {two_s})"]
    return lines


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
    for line in state_report(mol, per_atom):
        print("[state] " + line)
    ok, reasons = state_check(mol, per_atom)
    print("[state] intended 3x high-spin Fe(III) state: " + ("YES" if ok else "NO -- " + "; ".join(reasons)))


if __name__ == "__main__":
    main()
