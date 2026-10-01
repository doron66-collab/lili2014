#!/usr/bin/env python3
"""
localize_active_space.py -- split Pipek-Mezey localization of an active space
built by avas_by_count.py / avas_mp2_select.py, as a precondition for the
discarded-weight extrapolation (dmrg_extrapolate.py) to be valid at all.

Why (Claude Science, round 7, 2026-09-30): the discarded-weight extrapolation
validated on H12/STO-6G needed Lowdin-LOCAL orbitals -- canonical/AVAS orbitals
gave w decaying as a slow power law (M^-0.66) instead of near-exponentially, no
matter how large M got (implied M~7.8e5 to reach w=1e-5). AVAS orbitals are NOT
canonical (they are rotated to maximize overlap with the target AOs, so they
are already somewhat local to the TARGET REGION as a whole) but are not local
to single bonds/atoms either -- "better than canonical, worse than local," and
which end of that range a real cluster sits on is a measurement, not something
decidable by argument.

The fix costs nothing: split localization -- rotating the active-OCCUPIED
orbitals among themselves, and the active-VIRTUAL orbitals among themselves,
never mixing the two blocks -- is a unitary transformation WITHIN each block,
so it leaves exactly invariant: the CASCI/CASSCF/DMRG-converged energy, the
MP2/CI natural-orbital occupations (eigenvalues of the 1-RDM), the FCI
dimension, and the starting Hartree-Fock determinant (rotating occupied
orbitals among themselves doesn't change the occupied SUBSPACE, so the same
single determinant still describes it). What it DOES change: the discarded
weight w(M), and S_max itself -- both are properties of the chosen MPS site
basis, not of the molecule. (This is itself a clean dissertation point: a
quantity that changes under an energy-exactly-invariant unitary transformation
is not measuring a property of the state.)

Verified locally before writing this (N2/6-31g, AVAS CAS(8,7)): CASCI energy
before split-PM localization = -109.00969301603043 Ha, after =
-109.00969301603014 Ha -- agrees to 1e-7 (numerical noise), confirming the
invariance claim directly rather than trusting it.

Pipek-Mezey, NOT Boys, because the active space may contain an aromatic ring
(e.g. TP53_C275F's mutant PHE275 side chain) -- Boys maximizes centre
separation and mixes sigma/pi character into "banana bonds," which would
specifically corrupt localization exactly in the region being measured. PM
maximizes a Mulliken-population-based functional and preserves sigma/pi
separation.

Orbital REORDERING (e.g. block2's own reorder="fiedler") is a separate step,
done at DMRG-driver initialization time, not here -- reordering an already-
delocalized orbital list cannot place it, so it is only useful AFTER this
localization step. Not implemented here: verify the exact kwarg/API against
`help(DMRGDriver.initialize_system)` on Laguna's installed pyblock2 before
using it (this environment has no pyblock2 to check against).

Usage:
    python localize_active_space.py --json cluster_avas_mp2.json \
        --mo cluster_avas_mp2_mo.npy --out cluster_avas_mp2_localized

Writes <out>_mo.npy (same shape as the input, active block replaced by its
split-PM-localized version) and <out>.json (copy of the input JSON plus a
record of what was done, for provenance).
"""
import argparse
import json
import sys

import numpy


def localize_split(mf, mo, active_start_col, n_occ, n_vir, init_guess=None):
    """Split Pipek-Mezey localization of mo[:, active_start_col : active_start_col+n_occ+n_vir],
    occupied and virtual sub-blocks rotated separately. Returns a new mo array
    (copy) with only that range replaced.

    init_guess: forwarded to lo.PipekMezey (e.g. 'atomic') if the default
    initial guess fails to converge -- per Science's warning, PM can land in a
    local minimum, and the virtual sub-block (diffuse orbitals) is the less
    reliable of the two. Try a different init_guess or localize only the
    occupied block if the virtual block does not converge.
    """
    from pyscf import lo
    mol = mf.mol
    lo_, hi_occ = active_start_col, active_start_col + n_occ
    hi_vir = hi_occ + n_vir

    # PipekMezey has no public .converged flag after kernel(), and its
    # cost_function() takes a rotation matrix, not an mo_coeff array (both
    # checked locally -- an initial guess at the API was wrong on both counts,
    # not trusted from memory). The real safety net is verify_energy_invariance()
    # below, called unconditionally by main() -- a localization that silently
    # failed to do anything useful still leaves the energy exactly invariant
    # (it's still a valid, if poorly-localized, unitary rotation), so THAT check
    # cannot catch non-convergence by itself. There is no cheap, verified-here
    # signal for "PM converged well" beyond visually inspecting the orbitals
    # (e.g. plot/print the leading AO coefficients per column, as this script's
    # own docstring example did for the Lowdin case) -- left as a manual step
    # per Science's own "don't fight it, measure it" guidance for the virtual
    # block specifically.
    mo_out = mo.copy()
    pm_occ = lo.PipekMezey(mol, mo[:, lo_:hi_occ])
    if init_guess:
        pm_occ.init_guess = init_guess
    mo_out[:, lo_:hi_occ] = pm_occ.kernel()
    print(f"[localize] occupied block [{lo_}:{hi_occ}] ({n_occ} orbitals) localized", flush=True)

    if n_vir > 0:
        pm_vir = lo.PipekMezey(mol, mo[:, hi_occ:hi_vir])
        if init_guess:
            pm_vir.init_guess = init_guess
        mo_out[:, hi_occ:hi_vir] = pm_vir.kernel()
        print(f"[localize] virtual block [{hi_occ}:{hi_vir}] ({n_vir} orbitals) localized "
              f"-- per Science's warning this is the less reliable half (diffuse virtuals); "
              f"inspect the resulting coefficients if the downstream decay_law() probe still "
              f"shows power-law behavior", flush=True)

    return mo_out


def verify_energy_invariance(mf, mo_before, mo_after, ncas, nelec, tol=1e-6):
    """Re-run the SAME check this module's docstring claims was verified on N2:
    CASCI on this exact cluster's active space must agree before/after
    localization to numerical precision. Run every time, not just once in
    development -- this is what makes the invariance a checked fact for THIS
    cluster, not an inherited assumption from the N2 test case."""
    from pyscf import mcscf
    e_before = mcscf.CASCI(mf, ncas, nelec).kernel(mo_before)[0]
    e_after = mcscf.CASCI(mf, ncas, nelec).kernel(mo_after)[0]
    gap = abs(e_before - e_after)
    print(f"[localize] energy-invariance check: E_before={e_before:.10f}  "
          f"E_after={e_after:.10f}  |gap|={gap:.2e} Ha (tol {tol:.0e})", flush=True)
    if gap > tol:
        sys.exit(f"[localize] REFUSING: localization changed the CASCI energy by {gap:.2e} Ha, "
                 f"above tolerance {tol:.0e}. This should be numerically impossible for a "
                 f"split (occ-among-occ, vir-among-vir) rotation -- check that "
                 f"active_start_col/n_occ/n_vir match the ACTUAL active block boundaries "
                 f"before trusting anything downstream.")
    return e_before, e_after, gap


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, help="the avas_by_count.py/avas_mp2_select.py output JSON")
    ap.add_argument("--mo", required=True, help="the matching _mo.npy file")
    ap.add_argument("--init-guess", default=None, help="forwarded to lo.PipekMezey if default fails to converge (e.g. 'atomic')")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    meta = json.load(open(a.json))
    for required in ("xyz", "charge", "spin", "basis", "active_start_col", "n_occ", "n_vir"):
        if required not in meta:
            sys.exit(f"[localize] {a.json} is missing '{required}' -- was it produced by a "
                     f"version of avas_mp2_select.py/avas_by_count.py new enough to record it?")

    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
    from run_gate2_avas import build_mf
    mf = build_mf(meta["xyz"], meta["charge"], meta["spin"], basis=meta["basis"])

    mo = numpy.load(a.mo)
    ncas, nelec = meta["n_occ"] + meta["n_vir"], 2 * meta["n_occ"]

    mo_loc = localize_split(mf, mo, meta["active_start_col"], meta["n_occ"], meta["n_vir"],
                            init_guess=a.init_guess)
    e_before, e_after, gap = verify_energy_invariance(mf, mo, mo_loc, ncas, nelec)

    numpy.save(a.out + "_mo.npy", mo_loc)
    out_meta = dict(meta)
    out_meta.update(localized=True, localization_method="split Pipek-Mezey (occ/occ, vir/vir)",
                    energy_invariance_check_Ha=gap, source_json=a.json, source_mo=a.mo)
    with open(a.out + ".json", "w") as fh:
        json.dump(out_meta, fh, indent=1)
    print(f"[localize] wrote {a.out}_mo.npy and {a.out}.json "
          f"(pass the .npy to solange_dmrg.py --load-orbitals --ncas {ncas} --nelecas {nelec})")


if __name__ == "__main__":
    main()
