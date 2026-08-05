#!/usr/bin/env python
"""Calibrate the prefactor C in  M_required ≈ C · e^(S_max).

Why this exists
---------------
The platform reasons about classical tractability through two quantities: the
entanglement S_max a state carries, and the bond dimension M a tensor-network
method needs to represent it. The relation between them,

    M  ≈  C · e^S

is what turns a measured entropy into a predicted cost. Everywhere it is used
today, only *ratios* are legitimate, because C has never been measured — the
prefactor cancels in a ratio and does not cancel in an absolute statement. That
is why the affordability ceiling on the tractability plane is drawn with a
derived shape and an uncalibrated height, and why a published bond dimension
cannot be inverted to recover the entropy behind it.

This script measures C.

The design
----------
Every DMRG run already yields both numbers from one calculation: S_max is read
off the converged MPS, and the bond dimension actually required for chemical
accuracy is read off the energy-vs-M ladder. One run is one (S_max, M_required)
pair. To fit a line through such pairs, S has to vary while everything else is
held still — otherwise the orbital-count dependence contaminates the fit.

Stretching N2 is the cleanest available knob for exactly that. Same molecule,
same active space, same basis, same orbital count; only the correlation strength
moves, sweeping S from weakly correlated at equilibrium to strongly correlated
at dissociation. It is also the system the platform's own S_HARD threshold was
calibrated against (S_max = 0.42 at equilibrium, 2.86 stretched), so the fit
lands on ground the rest of the codebase already stands on.

Taking logs of the relation gives a straight line:

    ln M_required  =  ln C  +  S_max

so the fit returns two things, and the first matters more than the second:

  * the SLOPE, which should come out at 1. This tests whether the exponential
    relation holds at all. A slope far from 1 would mean the model is wrong, and
    no value of C would rescue it.
  * the INTERCEPT, which is ln C. This is the calibration.

Scope
-----
C is not a universal constant. It absorbs the shape of the entanglement
spectrum's decaying tail, which differs between chemistries, and the accuracy
demanded (here: chemical accuracy, 1.6 mHa). One molecular family therefore
fixes a value and an order of magnitude, and checks the slope — it does not
license transferring C to an arbitrary target. Fitting a second, chemically
different family is what would make a transferability claim defensible.

Usage
-----
    bash scripts/laguna/with_block2.sh \
        python scripts/laguna/calibrate_bond_dimension.py --out out/calibration.json

Small systems throughout — this does not need a large allocation and is not
affected by the per-user memory ceiling that blocks the full-site runs.
"""
import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))

from solange_dmrg import run_dmrg, CHEM_ACC_MHA          # noqa: E402

# Geometric ladder. Fine enough that ln(M_required) is resolved rather than
# quantised into a few steps, which a coarse 250/500/1000/2000 ladder would do —
# the fit is in log space, so resolution there is what the slope estimate needs.
# The floor starts low because the first calibration run found every geometry
# converged at or below M=28: the requirement was pinned against the bottom of
# the ladder, not measured.
DEFAULT_BOND_DIMS = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536]

# N2 bond lengths in Angstrom. 1.098 is equilibrium; the rest walk out toward
# dissociation, where the sigma bond breaks and static correlation takes over.
DEFAULT_LENGTHS = [1.098, 1.30, 1.60, 2.00, 2.50, 3.00]


_VALIDATED = [False]


def build_integrals(r, basis, ncas, nelecas, verbose=0, dmrg_scf=False,
                    dmrg_scf_maxm=250, dmrg_scf_scratch="./tmp_calib_orb", n_threads=4,
                    max_cycle_macro=None):
    """CASSCF at bond length r, returning the active-space integrals.

    CASSCF rather than CASCI on RHF orbitals: at stretched geometries the RHF
    reference is qualitatively wrong, and orbitals optimised in the presence of
    the active space are what make the S_max at each point comparable to the
    others. Convergence is reported per point rather than assumed.
    """
    from pyscf import gto, scf, mcscf, ao2mo

    mol = gto.M(atom=f"N 0 0 0; N 0 0 {r}", basis=basis, verbose=verbose, spin=0, charge=0)
    mf = scf.RHF(mol).run()
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    if dmrg_scf:
        # Past ~16 active orbitals CASSCF's FCI solver is combinatorially out of
        # reach, and a large active space is exactly what this calibration needs:
        # while the space is small enough for DMRG to be near-exact at trivial
        # bond dimension, M_required is set by the Hilbert space and carries no
        # information about S. Same adapter the classifier rung uses, and the
        # same tolerances — pyscf's defaults sit below DMRG's own RDM noise floor
        # and can never be met.
        from dmrgscf_block2 import Block2FCISolver, validate
        if not _VALIDATED[0]:
            validate(scratch=dmrg_scf_scratch + "_validate", n_threads=n_threads,
                     verbose=True)
            _VALIDATED[0] = True
        mc.fcisolver = Block2FCISolver(maxM=dmrg_scf_maxm, scratch=dmrg_scf_scratch,
                                       n_threads=n_threads)
        mc.internal_rotation = True
        mc.conv_tol = 1e-6
        mc.conv_tol_grad = 1e-3
    else:
        mc.fix_spin_(ss=0)
    # CASSCF calls its solver roughly four times per macro-iteration. With FCI
    # that is cheap and a high cap costs nothing; with DMRG each call is a full
    # sweep set, so the same cap means up to four hundred sweeps per geometry and
    # six geometries behind it. Default accordingly, and report convergence either
    # way rather than inferring it from the run having stopped.
    mc.max_cycle_macro = max_cycle_macro if max_cycle_macro is not None else (30 if dmrg_scf else 100)
    mc.run()

    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), ncas)
    return {
        "h1e": h1e, "h2e": h2e, "ecore": float(ecore),
        "e_casscf": float(mc.e_tot), "casscf_converged": bool(mc.converged),
    }


def required_bond_dim(energies, tol_mha=CHEM_ACC_MHA):
    """Smallest M whose energy is within chemical accuracy of the converged value.

    The reference is the largest M actually run. Two failure modes are reported
    rather than silently folded into the fit:

      at_floor   — even the smallest M on the ladder already qualified, so the
                   true requirement lies below the ladder and this point only
                   bounds it from above.
      at_ceiling — nothing but the reference itself qualified, so the run never
                   converged and the point carries no requirement at all.
    """
    if len(energies) < 2:
        return None, "insufficient"
    e_ref = energies[-1][1]
    for i, (M, e) in enumerate(energies):
        if abs(e - e_ref) * 1000.0 < tol_mha:
            if i == 0:
                return M, "at_floor"
            if i == len(energies) - 1:
                return M, "at_ceiling"
            return M, "ok"
    return None, "at_ceiling"


def fit(points):
    """Least squares of ln(M_required) on S_max. Returns slope, C, R²."""
    S = np.array([p["s_max"] for p in points], dtype=float)
    lnM = np.array([math.log(p["m_required"]) for p in points], dtype=float)
    slope, intercept = np.polyfit(S, lnM, 1)
    pred = slope * S + intercept
    ss_res = float(np.sum((lnM - pred) ** 2))
    ss_tot = float(np.sum((lnM - np.mean(lnM)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), math.exp(float(intercept)), r2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lengths", default=",".join(str(x) for x in DEFAULT_LENGTHS),
                    help="comma-separated N2 bond lengths in Angstrom")
    ap.add_argument("--basis", default="cc-pVDZ")
    ap.add_argument("--ncas", type=int, default=12,
                    help="active orbitals (default 12). NOT the full-valence 8: at "
                         "eight orbitals the middle cut spans only 4^4 states, so "
                         "the bond dimension DMRG needs is set by the size of the "
                         "Hilbert space rather than by the entanglement, and the "
                         "exponential relation this script measures cannot appear. "
                         "The first run confirmed exactly that — every geometry "
                         "landed on the same M. Raise further if M_required still "
                         "does not vary across the series.")
    ap.add_argument("--nelecas", type=int, default=10,
                    help="active electrons (default 10 — full valence for N2)")
    ap.add_argument("--bond-dims", default=",".join(str(x) for x in DEFAULT_BOND_DIMS))
    ap.add_argument("--dmrg-scf", action="store_true",
                    help="optimise the orbitals with DMRG instead of FCI, so --ncas "
                         "can exceed the ~16-orbital exact-solver wall. Needed for a "
                         "clean fit: at small --ncas the bond dimension DMRG requires "
                         "is capped by the Hilbert space rather than driven by S.")
    ap.add_argument("--dmrg-scf-maxm", type=int, default=250)
    ap.add_argument("--max-cycle-macro", type=int, default=None,
                    help="CASSCF orbital macro-iteration cap (default: 30 under "
                         "--dmrg-scf, 100 otherwise). Each macro-iteration is about "
                         "four solver calls, which is free under FCI and a full sweep "
                         "set each under DMRG.")
    ap.add_argument("--dmrg-scf-scratch", default="./tmp_calib_orb")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--scratch", default="./tmp_calib")
    ap.add_argument("--out", default="./out/calibration.json")
    ap.add_argument("--verbose", type=int, default=0)
    args = ap.parse_args()

    lengths = [float(x) for x in args.lengths.split(",") if x.strip()]
    bond_dims = [int(x) for x in args.bond_dims.split(",") if x.strip()]

    print("=" * 74)
    print("Bond-dimension calibration — N2 stretch series")
    print(f"  active space CAS({args.nelecas},{args.ncas}) / {args.basis}")
    print(f"  bond dims    {bond_dims[0]}..{bond_dims[-1]} ({len(bond_dims)} steps)")
    print(f"  lengths      {lengths}")
    print("=" * 74)

    points, t0 = [], time.time()
    for r in lengths:
        print(f"\n--- N2 at r = {r:.3f} A " + "-" * 40)
        try:
            cas = build_integrals(r, args.basis, args.ncas, args.nelecas, args.verbose,
                                  dmrg_scf=args.dmrg_scf, dmrg_scf_maxm=args.dmrg_scf_maxm,
                                  dmrg_scf_scratch=args.dmrg_scf_scratch,
                                  n_threads=args.threads,
                                  max_cycle_macro=args.max_cycle_macro)
        except Exception as exc:                          # noqa: BLE001
            print(f"  SKIP — CASSCF failed: {exc}")
            continue
        if not cas["casscf_converged"]:
            print("  WARNING — orbital optimisation hit the macro-iteration cap "
                  "without converging. The point is flagged, not dropped: raise "
                  "--max-cycle-macro if several geometries report this, since "
                  "unconverged orbitals make this geometry's S_max not strictly "
                  "comparable to the others.")

        # early_stop must be OFF: the whole measurement is where on the ladder the
        # energy settles, and early stop would truncate the ladder at exactly that
        # point, leaving no reference beyond it to measure the requirement against.
        # One scratch directory per (active space, geometry), so block2 never
        # resumes an MPS built for a different system. Keying on the geometry
        # alone is not enough: block2 is designed to pick up a saved state from a
        # reused scratch path, so a CAS(10,16) run would land on whatever a
        # CAS(10,12) run left behind at the same bond length. It must be created
        # here too — block2 creates only the final component of a path, not the
        # parent chain.
        scratch = os.path.join(args.scratch, f"cas{args.nelecas}_{args.ncas}", f"r{r:.3f}")
        os.makedirs(scratch, exist_ok=True)

        energies, s_max, stop_reason = run_dmrg(
            cas["h1e"], cas["h2e"], cas["ecore"], args.ncas, args.nelecas, bond_dims,
            scratch=scratch, n_threads=args.threads, early_stop=False)

        if s_max is None:
            print("  SKIP — S_max unavailable from the converged MPS")
            continue
        m_req, quality = required_bond_dim(energies)
        print(f"  S_max = {s_max:.4f}   M_required = {m_req}   ({quality})"
              f"   E = {energies[-1][1]:.8f} Ha")
        points.append({
            "r_angstrom": r, "s_max": float(s_max), "m_required": m_req,
            "quality": quality, "e_casscf": cas["e_casscf"],
            "casscf_converged": cas["casscf_converged"],
            "e_dmrg_final": float(energies[-1][1]),
            "energies": [[int(M), float(e)] for M, e in energies],
            "stop_reason": stop_reason,
        })

    usable = [p for p in points if p["quality"] == "ok"]
    print("\n" + "=" * 74)
    print(f"{len(points)} point(s) run, {len(usable)} usable for the fit "
          f"({len(points) - len(usable)} at the ladder's floor or ceiling)")

    if points:
        s_vals = [p["s_max"] for p in points]
        m_vals = [p["m_required"] for p in points if p["m_required"]]
        print(f"  S_max spanned      {min(s_vals):.3f} .. {max(s_vals):.3f}")
        print(f"  M_required spanned {min(m_vals)} .. {max(m_vals)} "
              f"({len(set(m_vals))} distinct value(s))")

    # A fit needs the dependent variable to actually move. If every geometry
    # converged at the same bond dimension, the regression returns slope 0 and an
    # undefined R^2 — which reads like a failed model but is really a failed
    # measurement, and the two must not be confused.
    distinct = len({p["m_required"] for p in usable})

    # A series can clear the "3 distinct values" bar and still be worthless: the
    # requirement rises over the easy geometries, then pins against the largest
    # bond dimension the active space can even use and stops responding to S
    # entirely. Fitting through that plateau averages a real slope together with
    # a flat tail and reports neither. Detect it, and fit the responsive part.
    plateau = []
    if usable:
        m_top = max(p["m_required"] for p in usable)
        plateau = [p for p in usable if p["m_required"] == m_top]
    if usable and distinct >= 3 and len(plateau) >= 3:
        responsive = [p for p in usable if p not in plateau] + plateau[:1]
        print(f"\n  M_required PINNED at {m_top} for {len(plateau)} of {len(usable)} "
              "geometries. Past that point the active space, not the entanglement, "
              "is setting the requirement: DMRG is already exact there and cannot "
              "need more than the Hilbert space at the cut allows. A fit across the "
              "plateau mixes a real slope with a flat tail.")
        if len(responsive) >= 3:
            r_slope, _, r_C, r_r2 = fit(responsive)
            print(f"  Over the {len(responsive)} responsive point(s) only: "
                  f"slope = {r_slope:.3f}, C = {r_C:.2f}, R^2 = {r_r2:.4f}")
        print("  To measure this properly the space must be large enough that "
              "M_required never reaches its own ceiling: re-run with a bigger "
              "--ncas (16, then 20) together with --dmrg-scf, which is what lets "
              "the orbital optimisation past the exact solver's ~16-orbital wall.")
        result_note = "requirement_pinned_active_space_bounds_it"
    elif usable and distinct < 3:
        print("\n  MEASUREMENT FAILED, NOT THE MODEL. Only "
              f"{distinct} distinct M_required value(s) across the series, so there "
              "is nothing for a slope to be fitted to. This is the signature of an "
              "active space too small for the relation to appear: when the middle "
              "cut spans few enough states, DMRG is near-exact at trivial bond "
              "dimension no matter how entangled the state is, and the requirement "
              "tracks the Hilbert space instead of the entropy. Raise --ncas (12 -> "
              "16 -> 20) until M_required varies across the geometries.")
        result_note = "measurement_failed_active_space_too_small"
    else:
        result_note = None

    result = {
        "system": "N2 stretch series",
        "active_space": f"CAS({args.nelecas},{args.ncas})",
        "basis": args.basis,
        "bond_dims": bond_dims,
        "tol_mha": CHEM_ACC_MHA,
        "points": points,
        "wall_seconds": round(time.time() - t0, 1),
    }

    if len(usable) >= 3 and distinct >= 3 and len(plateau) < 3:
        slope, intercept, C, r2 = fit(usable)
        result["fit"] = {"slope": slope, "intercept": intercept, "C": C, "r2": r2,
                         "n_points": len(usable)}
        print(f"\n  ln(M_required) = {intercept:.4f} + {slope:.4f} * S_max     (R^2 = {r2:.4f})")
        print(f"  slope = {slope:.3f}   — the model predicts 1; this is the test that it holds")
        print(f"  C     = {C:.2f}      — the calibrated prefactor")
        if abs(slope - 1.0) > 0.35:
            print("\n  WARNING: slope is far from 1. The exponential relation is not "
                  "holding over this range, and C should NOT be used until that is "
                  "understood — a fitted intercept under a wrong model is meaningless.")
    else:
        result["fit"] = None
        print("\n  NOT ENOUGH USABLE POINTS TO FIT (need 3). Widen --bond-dims at the "
              "end that failed: 'at_floor' means the ladder starts too high, "
              "'at_ceiling' means it stops too low.")

    if result_note:
        result["note"] = result_note
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(f"\n  wrote {out}")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
