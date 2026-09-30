#!/usr/bin/env python3
"""
dmrg_extrapolate.py -- discarded-weight extrapolation as the convergence criterion.

Why this replaces an energy difference between two bond dimensions:
a DeltaE between M=250 and M=500 conflates two things - the variational error of
the truncated MPS, and the sweep optimiser not having reached the minimum at that
M. The first shrinks predictably with the discarded weight; the second does not
shrink at all and can be arbitrarily large from a cold random start. That is the
documented R175H artifact. The discarded weight measures only the first, so
extrapolating against it separates them.

The relation used is the standard one: to leading order the truncation error in
the energy is linear in the discarded weight w,

    E(w) ~= E_0 + c * w ,   c > 0

so an ordinary least-squares fit of E against w over the largest-M points, taken
to w -> 0, estimates the untruncated energy. Convergence is then declared when
the best computed energy and the extrapolated one agree inside chemical accuracy.

The fit is only meaningful when w spans enough range; this module refuses to
report a verdict when it does not, rather than returning a confident number from
three points sitting on top of each other.

Nothing here depends on how w is obtained -- pass it in. Two routes, in order of
preference: a numeric attribute on the driver (check help() on your installed
pyblock2), or parse_block2_log() below, which reads the sweep output and always
works.
"""
import re
import sys

CHEM_ACC_HA = 0.0016          # 1 kcal/mol, the bar used throughout this project


def _ols(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    if sxx == 0.0:
        return None, None, None
    slope = sxy / sxx
    icept = my - slope * mx
    ss_tot = sum((yi - my) ** 2 for yi in y)
    ss_res = sum((yi - (icept + slope * xi)) ** 2 for xi, yi in zip(x, y))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return icept, slope, r2


def extrapolate(bond_dims, energies, discarded_weights,
                n_points=None, chem_acc=CHEM_ACC_HA, min_w_decades=1.0):
    """Extrapolate E(w) -> w=0.

    bond_dims, energies, discarded_weights : equal-length sequences, one entry per
        bond dimension, each taken from the LAST sweep at that bond dimension
        (a mid-sweep value is not a converged energy at that M).
    n_points : how many of the largest-M points to fit. Default: all of them if
        4 or fewer, otherwise the largest 4.
    min_w_decades : refuse a verdict unless w spans at least this many decades.
    """
    pts = sorted(zip(bond_dims, energies, discarded_weights), key=lambda t: t[0])
    if len(pts) < 3:
        return {"ok": False, "reason": "need at least 3 bond dimensions, got %d" % len(pts)}
    if any(w is None or w <= 0 for _, _, w in pts):
        return {"ok": False, "reason": "non-positive or missing discarded weight in the series"}

    k = n_points or min(4, len(pts))
    use = pts[-k:]
    M = [p[0] for p in use]
    E = [p[1] for p in use]
    W = [p[2] for p in use]

    import math
    decades = math.log10(max(W) / min(W))
    e0, slope, r2 = _ols(W, E)
    if e0 is None:
        return {"ok": False, "reason": "all discarded weights identical; cannot fit"}

    e_best = E[-1]
    resid = abs(e_best - e0)
    warn = []
    if max(W) < 1e-12:
        # Caught from a real failed validation: a field returning values at machine
        # epsilon can still produce a high R^2 by coincidence, with a fitted slope of
        # order 1e16 Ha per unit weight. A good fit on numerically absent truncation
        # is an artifact, not evidence.
        warn.append("largest discarded weight is %.2e -- truncation is numerically absent, "
                    "so there is nothing to extrapolate against. Either the bond dimensions "
                    "already span the exact MPS, or this field is not the discarded weight "
                    "(a fitted slope far above ~1e2 Ha per unit weight is the tell)" % max(W))
    if decades < min_w_decades:
        warn.append("discarded weight spans only %.2f decades over the fitted points; "
                    "extend the bond-dimension ladder before trusting the intercept" % decades)
    if slope <= 0:
        warn.append("fitted slope is %.3e, not positive -- the energy is not decreasing with "
                    "decreasing truncation, which means the sweeps are not converged at fixed M "
                    "(raise --n-sweeps) rather than that truncation is under control" % slope)
    if r2 < 0.95:
        warn.append("R^2 = %.3f over %d points; the linear-in-w regime has not been reached" % (r2, k))
    if E != sorted(E, reverse=True):
        warn.append("energies are not monotonically decreasing with increasing M -- a variational "
                    "method cannot do that at fixed sweep convergence; suspect an unconverged sweep")

    converged = (resid < chem_acc) and not warn
    return {"ok": True,
            "fitted_points": k, "bond_dims_used": M,
            "w_min": min(W), "w_max": max(W), "w_decades": decades,
            "e_extrapolated_Ha": e0, "slope_Ha_per_w": slope, "r2": r2,
            "e_best_computed_Ha": e_best, "bond_dim_best": M[-1],
            "residual_Ha": resid, "residual_mHa": resid * 1000.0,
            "chem_acc_Ha": chem_acc,
            "converged": converged,
            "verdict": ("CONVERGED: best computed energy is within chemical accuracy of the "
                        "w->0 limit" if converged else
                        "NOT CONVERGED at this ladder" if not warn else
                        "INCONCLUSIVE: see warnings, do not classify on this fit"),
            "warnings": warn}


def last_sweep_per_bond_dim(bond_dims, dws, energies, root=0):
    """Reduce the three PER-SWEEP arrays from DMRGDriver.get_dmrg_results() to one
    (M, E, w) triple per bond dimension, taking the LAST sweep at each M.

    Settled against the block2 driver source (pyblock2/driver/core.py,
    get_dmrg_results): all three arrays are per SWEEP, and dws is documented as
    "the maximal discarded weight (sum of discarded eigenvalues) for each sweep".
    The energies entries may each be a list (one per root), hence `root`.

    Taking the last sweep at each M - rather than the maximum weight over the
    sweeps at that M - is what keeps E and w on the same sweep. The early sweeps
    after a bond-dimension increase carry a larger weight because the MPS has not
    adapted yet; pairing that weight with the converged energy mixes two states.
    """
    def scalar(e):
        try:
            return float(e[root])
        except (TypeError, IndexError):
            return float(e)

    n = min(len(bond_dims), len(dws), len(energies))
    last = {}
    for i in range(n):
        last[int(bond_dims[i])] = (scalar(energies[i]), float(dws[i]))
    out = sorted((m, e, w) for m, (e, w) in last.items())
    return [m for m, _, _ in out], [e for _, e, _ in out], [w for _, _, w in out]


def window_stability(bond_dims, energies, discarded_weights, chem_acc=CHEM_ACC_HA):
    """Is the ladder actually in the asymptotic linear-in-w regime?

    Extrapolation to w -> 0 is a LOCAL statement about the tail, not a global fit.
    Fitting points where the MPS is still qualitatively wrong cannot be linear and
    will return a confident, wrong intercept. The test: refit over progressively
    shorter windows from the large-M end and see whether the intercept moves.
    If dropping the smallest-M point shifts it by more than chemical accuracy, the
    ladder has not reached the regime - extend it upward rather than trusting it.
    """
    rows, n = [], len(bond_dims)
    for k in range(3, n + 1):
        r = extrapolate(bond_dims, energies, discarded_weights, n_points=k,
                        chem_acc=chem_acc)
        if r.get("ok"):
            rows.append({"n_points": k, "bond_dims": r["bond_dims_used"],
                         "e_extrapolated_Ha": r["e_extrapolated_Ha"],
                         "r2": r["r2"], "slope": r["slope_Ha_per_w"]})
    if len(rows) < 2:
        return {"ok": False, "reason": "need at least 4 bond dimensions to test stability"}
    es = [r["e_extrapolated_Ha"] for r in rows]
    spread = max(es) - min(es)
    return {"ok": True, "windows": rows, "intercept_spread_Ha": spread,
            "intercept_spread_mHa": spread * 1000.0,
            "stable": spread < chem_acc,
            "verdict": ("intercept is stable across fit windows; the tail is in the "
                        "linear regime" if spread < chem_acc else
                        "intercept moves by %.2f mHa across fit windows -- the ladder is "
                        "NOT in the asymptotic regime; extend it to larger M rather than "
                        "trusting any single fit" % (spread * 1000.0))}


def parse_block2_log(path, dw_token=r'DW'):
    """Fallback route to the discarded weight: read it out of the sweep output.

    block2 prints a per-sweep summary carrying the bond dimension, the energy and
    the discarded weight. The exact field name has changed between releases, so
    VERIFY the token against your own log before relying on this -- print a few
    sweep lines first. Returns [(bond_dim, energy, dw), ...], one per bond
    dimension, taking the LAST sweep at each.
    """
    num = r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?'
    re_bd = re.compile(r'[Bb]ond\s*dimension\s*=\s*(\d+)')
    re_e = re.compile(r'\bE\s*=\s*(' + num + r')')
    re_dw = re.compile(dw_token + r'\s*=\s*(' + num + r')')
    last = {}
    for ln in open(path, encoding='utf-8', errors='replace'):
        mb, me, mw = re_bd.search(ln), re_e.search(ln), re_dw.search(ln)
        if mb and me and mw:
            last[int(mb.group(1))] = (float(me.group(1)), float(mw.group(1)))
    return [(bd, e, w) for bd, (e, w) in sorted(last.items())]


def _selftest():
    # synthetic series obeying E = E0 + c*w exactly
    E0_true, c_true = -5333.128000, 4.75
    W = [1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5]
    Ms = [250, 500, 1000, 2000]
    E = [E0_true + c_true * w for w in W]
    r = extrapolate(Ms, E, W)
    assert r["ok"], r
    assert abs(r["e_extrapolated_Ha"] - E0_true) < 1e-9, r
    assert abs(r["slope_Ha_per_w"] - c_true) < 1e-6, r
    assert r["r2"] > 0.999999
    print("selftest exact-linear   : E0 recovered to %.2e Ha, slope to %.2e, R2=%.6f"
          % (abs(r["e_extrapolated_Ha"] - E0_true), abs(r["slope_Ha_per_w"] - c_true), r["r2"]))
    print("                          residual %.4f mHa -> %s" % (r["residual_mHa"], r["verdict"]))

    # a cold-start-artifact series: energy rises with M, which is variationally impossible
    Ebad = [-5333.1280, -5333.1275, -5333.1290, -5333.1230]
    rb = extrapolate(Ms, Ebad, W)
    assert not rb["converged"] and rb["warnings"], rb
    print("selftest cold-start fake: caught -> %s" % rb["warnings"][0][:88])

    # too narrow a w range
    rn = extrapolate([1000, 1500, 2000], [-1.0001, -1.00005, -1.00002],
                     [5.0e-5, 4.0e-5, 3.0e-5])
    assert any('decades' in w for w in rn["warnings"]), rn
    print("selftest narrow ladder  : caught -> %s" % [w for w in rn["warnings"] if 'decades' in w][0][:88])

    # per-sweep reduction: 3 sweeps at each of 3 bond dimensions, energies as
    # one-element lists (the multi-root shape the driver documents)
    bd = [16, 16, 16, 32, 32, 32, 64, 64, 64]
    dw = [9e-4, 6e-4, 5e-4, 3e-4, 2e-4, 1.5e-4, 8e-5, 5e-5, 4e-5]
    en = [[-1.10], [-1.12], [-1.1300], [-1.14], [-1.15], [-1.1550], [-1.158], [-1.159], [-1.1595]]
    M, E, W = last_sweep_per_bond_dim(bd, dw, en)
    assert M == [16, 32, 64] and W == [5e-4, 1.5e-4, 4e-5], (M, E, W)
    assert E == [-1.13, -1.155, -1.1595], E
    print("selftest per-sweep glue : reduced 9 sweeps -> M=%s, last-sweep w=%s" % (M, W))

    # window stability: a ladder whose small-M end is outside the linear regime
    Mq = [8, 16, 32, 64, 128, 256]
    Wq = [2.0e-2, 4.0e-3, 8.0e-4, 1.6e-4, 3.2e-5, 6.4e-6]
    # asymptotic part obeys E0 + 3*w; the two smallest M are pulled far off the line
    Eq = [(-1.2000 + 3.0 * w) for w in Wq]
    Eq[0] += 0.35
    Eq[1] += 0.06
    st = window_stability(Mq, Eq, Wq)
    assert st["ok"] and not st["stable"], st
    print("selftest window stability: spread %.2f mHa across windows -> not asymptotic (correct)"
          % st["intercept_spread_mHa"])
    st2 = window_stability(Mq[2:], Eq[2:], Wq[2:])
    assert st2["ok"] and st2["stable"], st2
    print("                          : dropping the two bad points -> stable (%.4f mHa)"
          % st2["intercept_spread_mHa"])


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
        _selftest()
    else:
        print(__doc__)
        print("run with --selftest to verify the fit and the guards")
