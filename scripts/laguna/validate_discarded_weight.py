#!/usr/bin/env python3
"""
validate_discarded_weight.py -- confirm the discarded-weight extrapolation
against a KNOWN exact answer, using pyblock2's DOCUMENTED accessor.

History (2026-09-30, same day, three rounds):
  Round 1: probed drv._dmrg for anything matching disc/weight/error, found two
           candidates (discarded_weights, sweep_discarded_weights) with no
           documentation to say which was real.
  Round 2: tested both (plus max-over-array and second-to-last-sweep variants)
           against H8/CAS(8,8) vs exact FCI. ALL FOUR FAILED -- R^2 as low as
           0.0004, fitted slopes with the wrong sign or absurd magnitude.
  Round 3 (this version): Claude Science read pyblock2/driver/core.py directly
           and found the documented accessor, DMRGDriver.get_dmrg_results():
           returns (bond_dims, dws, energies), all per-SWEEP, with dws
           "the maximal discarded weight (sum of discarded eigenvalues) for
           each sweep". sweep_discarded_weights is untouched by the Python
           driver -- a C++-internal field the wrapper never populates, and
           round 2 should never have read it. discarded_weights (round 2's
           "discarded_weights_last_bond") WAS the right field and the right
           element all along -- round 2 failed for a DIFFERENT reason: H8 at
           CAS(8,8) has no asymptotic linear-in-w window at all (M=4/8 are
           qualitatively wrong, M=32/64 are already numerically exact -- a fit
           across that gap cannot be linear by construction).

This version fixes the TEST SYSTEM, not the field: a linear H14 chain (still
r=1.8 A, STO-6G, CAS(14,14)) has a much larger exact bond dimension, so
M=16..256 gives a real asymptotic window while FCI (dim = C(14,7)^2 ~ 1.18e7)
is still a cheap, exact PySCF reference.

Acceptance criteria (five now -- window stability added):
  1. w > 0 and monotonically DEcreasing as M increases
  2. E(M) - E_FCI > 0 and monotonically decreasing as M increases
  3. linear fit of E against w: positive slope, R^2 > 0.99
  4. |E0_extrapolated - E_FCI| < 0.1 mHa
  5. window_stability(): the extrapolated intercept must not move by more than
     chemical accuracy when the fit window is shortened from the large-M end --
     otherwise the ladder isn't actually in the linear regime the fit assumes,
     which is exactly how H8 gave a deceptively plausible R^2 in round 2's
     "sweep_discarded_weights_prelast" case (0.928, on values at machine
     epsilon -- see dmrg_extrapolate.extrapolate's own machine-epsilon guard,
     added because of that specific false positive).

Usage: python validate_discarded_weight.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BOND_DIMS = [16, 32, 64, 128, 256]
N_SITES = 14
BOND_LENGTH = 1.8
BASIS = "sto-6g"


def build_h14_fci():
    from pyscf import gto, scf, fci, ao2mo
    atoms = "\n".join(f"H 0 0 {i * BOND_LENGTH:.4f}" for i in range(N_SITES))
    mol = gto.M(atom=atoms, basis=BASIS, verbose=0)
    mf = scf.RHF(mol).run()
    h1e = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.kernel(mol, mf.mo_coeff)
    eri_full = ao2mo.restore(1, eri, mol.nao)
    na = mol.nelectron // 2
    print(f"[h14] running exact FCI at CAS({mol.nelectron},{mol.nao}) -- "
          f"dim=C({mol.nao},{na})^2, this is the expensive-but-exact step, "
          f"give it a few minutes...", flush=True)
    e_fci, _ = fci.direct_spin1.FCI().kernel(h1e, eri, mol.nao, (na, na), ecore=mf.energy_nuc())
    print(f"[h14] E_scf={mf.e_tot:.8f}  E_fci={e_fci:.8f}  CAS({mol.nelectron},{mol.nao})",
          flush=True)
    return h1e, eri_full, mf.energy_nuc(), mol.nao, mol.nelectron, e_fci


def run_ladder(h1e, eri_full, ecore, ncas, nelec):
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch="./tmp_validate_dw14", symm_type=SymmetryTypes.SU2,
                     n_threads=4, stack_mem=int(4.0 * (1 << 30)))
    drv.initialize_system(n_sites=ncas, n_elec=nelec, spin=0)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=eri_full, ecore=ecore, iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=min(BOND_DIMS[0], 100), nroots=1)

    energies, weights = [], []
    for M in BOND_DIMS:
        # Last sweep(s) at exactly noise=0 (Claude Science, round 3): with
        # noises=None block2 picks its own schedule, which can still carry
        # noise into the final sweep -- inflating both the reported energy
        # (sits above the true variational minimum) and the discarded weight
        # (computed from a noise-perturbed density matrix). Explicit here so
        # the LAST sweep's (E, w) pair is the clean, noise-free one the fit
        # assumes.
        e = drv.dmrg(mpo, ket, n_sweeps=16, bond_dims=[M],
                     noises=[1e-4, 1e-5, 1e-6, 1e-7] + [0.0] * 12,
                     thrds=[1e-10] * 16, cutoff=1e-20, iprint=0)
        bd_arr, dws, e_arr = drv.get_dmrg_results()
        M_this, E_this, W_this = None, None, None
        import dmrg_extrapolate as dex
        Ms_, Es_, Ws_ = dex.last_sweep_per_bond_dim(bd_arr, dws, e_arr)
        # last_sweep_per_bond_dim reduces the WHOLE accumulated history each
        # call passes it; since this loop calls dmrg() once per M with a fresh
        # single-M request, only M's own sweeps are present each time.
        idx = Ms_.index(M)
        E_this, W_this = Es_[idx], Ws_[idx]
        energies.append(E_this)
        weights.append(W_this)
        print(f"  M={M:4d}  E={E_this:.8f} Ha (dmrg() returned {float(e):.8f})  "
              f"w(last sweep)={W_this:.6e}  n_sweeps_actually_run={len(dws)}", flush=True)
    return energies, weights


def main():
    h1e, eri_full, ecore, ncas, nelec, e_fci = build_h14_fci()
    energies, weights = run_ladder(h1e, eri_full, ecore, ncas, nelec)

    import dmrg_extrapolate as dex
    print("=" * 72)
    fails = []
    if not all(weights[i] > weights[i + 1] for i in range(len(weights) - 1)):
        fails.append("w not monotonically decreasing with M")
    resid = [e - e_fci for e in energies]
    if not all(r > 0 for r in resid):
        fails.append("E - E_FCI not positive for all M")
    if not all(resid[i] > resid[i + 1] for i in range(len(resid) - 1)):
        fails.append("E - E_FCI not monotonically decreasing with M")

    r = dex.extrapolate(BOND_DIMS, energies, weights, min_w_decades=0.0)
    print(f"weights={weights}")
    print(f"residuals to FCI (Ha)={resid}")
    if r.get("ok"):
        gap_mha = abs(r["e_extrapolated_Ha"] - e_fci) * 1000.0
        print(f"extrapolated E0={r['e_extrapolated_Ha']:.8f} Ha  gap to FCI={gap_mha:.4f} mHa  "
              f"slope={r['slope_Ha_per_w']:.3e}  R2={r['r2']:.5f}")
        if r["slope_Ha_per_w"] <= 0:
            fails.append(f"slope not positive ({r['slope_Ha_per_w']:.3e})")
        if r["r2"] <= 0.99:
            fails.append(f"R^2 = {r['r2']:.4f}, not > 0.99")
        if gap_mha >= 0.1:
            fails.append(f"|E0 - E_FCI| = {gap_mha:.4f} mHa, not < 0.1 mHa")
        for w in r.get("warnings", []):
            print(f"  extrapolate() warning: {w}")
    else:
        fails.append(f"extrapolate() refused: {r.get('reason')}")

    stab = dex.window_stability(BOND_DIMS, energies, weights)
    if stab.get("ok"):
        print(f"window stability: spread {stab['intercept_spread_mHa']:.4f} mHa -> "
              f"{stab['verdict']}")
        if not stab["stable"]:
            fails.append(f"window NOT stable (spread {stab['intercept_spread_mHa']:.4f} mHa)")
    else:
        fails.append(f"window_stability() refused: {stab.get('reason')}")

    print("=" * 72)
    if fails:
        print("VERDICT: FAILS. Do not trust this ladder/field combination for a real "
              "classification yet:")
        for f in fails:
            print(f"   * {f}")
    else:
        print("VERDICT: PASSES all five acceptance criteria. discarded_weights via "
              "get_dmrg_results(), with the last-sweep reduction, is validated for "
              "real classification runs.")


if __name__ == "__main__":
    main()
