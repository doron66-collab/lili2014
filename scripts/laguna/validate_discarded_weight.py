#!/usr/bin/env python3
"""
validate_discarded_weight.py -- decide which pyblock2 attribute is the real
per-bond-dimension discarded weight, by reproducing a KNOWN answer.

Same discipline as avas_by_count.py's own validation (reproduce a known result
before trusting a new code path on real data). Here the known result is exact
FCI on a small but genuinely-truncating system, so a candidate discarded-weight
series either passes four independent acceptance checks against it or it
doesn't -- no guessing from attribute names.

Test system (Claude Science, 2026-09-30): a linear H8 chain, r=1.8 A (stretched,
strongly correlated -- CCSD(T)-breaking regime), STO-6G, CAS(8,8). FCI dimension
is C(8,4)^2 = 4900 -- exact and instant -- while M=4..8 genuinely truncates the
2^8-dimensional Hilbert space, so this is small enough to check exactly and
still has real truncation error to measure, unlike the earlier 2-orbital probe
(which was already exact at M=4 and could not discriminate anything).

Acceptance criteria, ALL four required for a candidate field to pass:
  1. w > 0 and monotonically DEcreasing as M increases
  2. E(M) - E_FCI > 0 and monotonically decreasing as M increases
  3. linear fit of E against w (dmrg_extrapolate.extrapolate): positive slope,
     R^2 > 0.99
  4. |E0_extrapolated - E_FCI| < 0.1 mHa -- far tighter than the project's
     1.6 mHa chemical-accuracy bar, because this checks the METHOD, not
     chemistry: it should reproduce a number we already know exactly.

Usage: python validate_discarded_weight.py
"""
import sys

CANDIDATE_FIELDS = ["discarded_weights", "sweep_discarded_weights"]
BOND_DIMS = [4, 8, 16, 32, 64]


def build_h8_fci():
    from pyscf import gto, scf, fci, ao2mo
    r = 1.8
    atoms = "\n".join(f"H 0 0 {i * r:.4f}" for i in range(8))
    mol = gto.M(atom=atoms, basis="sto-6g", verbose=0)
    mf = scf.RHF(mol).run()
    h1e = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.kernel(mol, mf.mo_coeff)
    eri_full = ao2mo.restore(1, eri, mol.nao)
    na = mol.nelectron // 2
    e_fci, _ = fci.direct_spin1.FCI().kernel(h1e, eri, mol.nao, (na, na), ecore=mf.energy_nuc())
    print(f"[h8] E_scf={mf.e_tot:.8f}  E_fci={e_fci:.8f}  CAS({mol.nelectron},{mol.nao})  "
          f"FCI dim=4900 (exact reference)")
    return h1e, eri_full, mf.energy_nuc(), mol.nao, mol.nelectron


def run_ladder(h1e, eri_full, ecore, ncas, nelec):
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch="./tmp_validate_dw", symm_type=SymmetryTypes.SU2,
                     n_threads=2, stack_mem=int(2.0 * (1 << 30)))
    drv.initialize_system(n_sites=ncas, n_elec=nelec, spin=0)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=eri_full, ecore=ecore, iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=min(BOND_DIMS[0], 50), nroots=1)

    per_field = {f: [] for f in CANDIDATE_FIELDS}
    energies = []
    for M in BOND_DIMS:
        e = float(drv.dmrg(mpo, ket, n_sweeps=12, bond_dims=[M],
                           noises=[1e-4, 1e-5, 1e-6, 0], thrds=[1e-10] * 4, iprint=0))
        energies.append(e)
        for field in CANDIDATE_FIELDS:
            try:
                val = getattr(drv._dmrg, field)
                per_field[field].append(float(val[-1]))
            except Exception as exc:
                per_field[field].append(None)
        print(f"  M={M:4d}  E={e:.8f} Ha  " +
              "  ".join(f"{f}={per_field[f][-1]}" for f in CANDIDATE_FIELDS))
    return energies, per_field


def check(field, bond_dims, energies, weights, e_fci):
    import dmrg_extrapolate as dex
    fails = []
    if any(w is None for w in weights):
        return None, ["field missing/unreadable"]
    if not all(weights[i] > weights[i + 1] for i in range(len(weights) - 1)):
        fails.append("w not monotonically decreasing with M")
    if any(w <= 0 for w in weights):
        fails.append("w not strictly positive")
    resid = [e - e_fci for e in energies]
    if not all(r > 0 for r in resid):
        fails.append("E - E_FCI not positive for all M")
    if not all(resid[i] > resid[i + 1] for i in range(len(resid) - 1)):
        fails.append("E - E_FCI not monotonically decreasing with M")
    r = dex.extrapolate(bond_dims, energies, weights, min_w_decades=0.0)
    if not r.get("ok"):
        fails.append(f"extrapolate() refused: {r.get('reason')}")
        return r, fails
    if r["slope_Ha_per_w"] <= 0:
        fails.append(f"slope not positive ({r['slope_Ha_per_w']:.3e})")
    if r["r2"] <= 0.99:
        fails.append(f"R^2 = {r['r2']:.4f}, not > 0.99")
    gap_mha = abs(r["e_extrapolated_Ha"] - e_fci) * 1000.0
    if gap_mha >= 0.1:
        fails.append(f"|E0 - E_FCI| = {gap_mha:.4f} mHa, not < 0.1 mHa")
    r["gap_to_fci_mHa"] = gap_mha
    return r, fails


def main():
    h1e, eri_full, ecore, ncas, nelec = build_h8_fci()
    from pyscf import fci as fci_mod, ao2mo
    # recompute e_fci alone for clarity of what's being checked against
    na = nelec // 2
    e_fci = fci_mod.direct_spin1.FCI().kernel(h1e, eri_full, ncas, (na, na), ecore=ecore)[0]

    energies, per_field = run_ladder(h1e, eri_full, ecore, ncas, nelec)

    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
    print("=" * 72)
    passed = []
    for field in CANDIDATE_FIELDS:
        weights = per_field[field]
        r, fails = check(field, BOND_DIMS, energies, weights, e_fci)
        print(f"[{field}]  weights={weights}")
        if r and r.get("ok"):
            print(f"    extrapolated E0={r['e_extrapolated_Ha']:.8f} Ha  "
                  f"gap to FCI={r.get('gap_to_fci_mHa', float('nan')):.4f} mHa  "
                  f"slope={r['slope_Ha_per_w']:.3e}  R2={r['r2']:.5f}")
        if fails:
            print(f"    FAILS: {fails}")
        else:
            print(f"    PASSES all four acceptance criteria")
            passed.append(field)
        print("-" * 72)

    print("=" * 72)
    if len(passed) == 1:
        print(f"VERDICT: use `{passed[0]}` (the last element per bond dimension) "
              f"as the discarded weight in solange_dmrg.py's run_dmrg().")
    elif len(passed) > 1:
        print(f"VERDICT: {passed} both pass — prefer 'sweep_discarded_weights' "
              f"(per-sweep semantics match what the extrapolation needs directly).")
    else:
        print("VERDICT: NEITHER candidate field passed. Do not wire either into "
              "run_dmrg() — fall back to dmrg_extrapolate.parse_block2_log() "
              "against a --verbose/iprint>=1 run instead, and verify that path "
              "the same way (against this same H8/FCI reference) before trusting it.")


if __name__ == "__main__":
    main()
