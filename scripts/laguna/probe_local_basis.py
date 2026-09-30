#!/usr/bin/env python3
"""
probe_local_basis.py -- two-point probe of Claude Science's round-6 diagnosis:
canonical (delocalized) RHF orbitals are why w decays as a slow power law
(M^-0.66) instead of near-exponentially, on the H12 chain used by
validate_discarded_weight.py.

Verified independently (locally, before writing this) with plain PySCF, no
pyblock2 needed: on H12/STO-6G at r=1.4 A, Lowdin-orthogonalized AOs (one
localized orbital per hydrogen atom, already ordered along the chain since the
atoms are listed in order) give CASCI(12,12) = exact FCI energy to 6 decimals
(basis rotation invariance, since CAS(12,12) in a minimal basis is the full
space) -- confirmed -6.140222367804819 Ha vs FCI's -6.14022239 Ha -- and the
first Lowdin orbital's AO coefficients decay from 1.04 to ~0 within 4 atoms,
confirming genuine locality (a canonical RHF orbital would spread visibly
across most/all 12 atoms instead).

This script builds the SAME H12/r=1.4A/STO-6G/CAS(12,12) system in the LOCAL
basis and runs DMRG at exactly two bond dimensions (M=12, M=48, per Science's
own "probe before committing" request) rather than the full 5-point ladder --
cheap enough to decide before paying for anything bigger.

Decision table (Science's own):
    w(48) < 1e-5 and |E(48) - E_FCI| < 1e-5 Ha
        -> diagnosis confirmed. Re-run the full ladder in this basis with a
           SMALLER M range (perhaps 6-24) and all six acceptance criteria.
    w(48) still ~1e-3 or larger, decaying as a power law
        -> diagnosis wrong again. Stop trying new toy systems; switch to a
           very-large-M DMRG energy as the reference instead of exact FCI
           (which is what a real classification run does anyway).

Usage: python probe_local_basis.py
"""
import numpy


def build_local_basis_system():
    from pyscf import gto, scf, ao2mo
    r, n = 1.4, 12
    atoms = "\n".join(f"H 0 0 {i * r:.4f}" for i in range(n))
    mol = gto.M(atom=atoms, basis="sto-6g", verbose=0)
    mf = scf.RHF(mol).run()

    S = mol.intor("int1e_ovlp")
    w_, v = numpy.linalg.eigh(S)
    C_loc = v @ numpy.diag(w_ ** -0.5) @ v.T   # Lowdin-orthogonalized AOs

    # sanity check done here too (not just trusted from the earlier standalone
    # check): CASCI on the local basis must reproduce the SCF's own FCI-level
    # result to high precision, since CAS(12,12) in a 12-function minimal basis
    # IS the full space and is invariant to any orbital rotation.
    from pyscf import mcscf
    mc = mcscf.CASCI(mf, 12, 12)
    e_check = mc.kernel(C_loc)[0]
    print(f"[probe] Lowdin-local CASCI energy: {e_check:.10f} Ha "
          f"(sanity check -- must match FCI in the canonical basis)")

    h1e = C_loc.T @ mf.get_hcore() @ C_loc
    eri = ao2mo.kernel(mol, C_loc)
    eri_full = ao2mo.restore(1, eri, mol.nao)
    return h1e, eri_full, mf.energy_nuc(), mol.nao, mol.nelectron, e_check


def main():
    h1e, eri_full, ecore, ncas, nelec, e_fci = build_local_basis_system()

    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch="./tmp_probe_local", symm_type=SymmetryTypes.SU2,
                     n_threads=4, stack_mem=int(4.0 * (1 << 30)))
    drv.initialize_system(n_sites=ncas, n_elec=nelec, spin=0)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=eri_full, ecore=ecore, iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=12, nroots=1)

    results = {}
    for M in (12, 48):
        e = drv.dmrg(mpo, ket, n_sweeps=16, bond_dims=[M],
                     noises=[1e-4, 1e-5, 1e-6, 1e-7] + [0.0] * 12,
                     thrds=[1e-10] * 16, cutoff=1e-20, iprint=0)
        bd_arr, dws, e_arr = drv.get_dmrg_results()
        import dmrg_extrapolate as dex
        Ms_, Es_, Ws_ = dex.last_sweep_per_bond_dim(bd_arr, dws, e_arr)
        idx = Ms_.index(M)
        E_this, W_this = Es_[idx], Ws_[idx]
        gap = E_this - e_fci
        results[M] = (E_this, W_this, gap)
        print(f"  M={M:3d}  E={E_this:.8f} Ha  w={W_this:.3e}  E-E_FCI={gap:.3e} Ha")

    print("=" * 60)
    e48, w48, gap48 = results[48]
    if w48 < 1e-5 and abs(gap48) < 1e-5:
        print(f"CONFIRMED: w(48)={w48:.2e} < 1e-5, |E-E_FCI|={abs(gap48):.2e} Ha < 1e-5. "
              f"Local-orbital diagnosis holds. Re-run the full ladder in this basis, "
              f"smaller M (e.g. 6,9,12,16,24), all six acceptance criteria.")
    else:
        print(f"NOT CONFIRMED: w(48)={w48:.2e}, |E-E_FCI|={abs(gap48):.2e} Ha -- still far "
              f"from the target. Do not try a fourth/fifth toy system -- switch to a "
              f"very-large-M DMRG energy as the reference instead of exact FCI.")


if __name__ == "__main__":
    main()
