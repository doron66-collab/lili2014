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

CHEM_ACC_HA = 0.0016
MAX_W_LINEAR = 1e-3   # ceiling on the discarded weight for the linear-in-w regime          # 1 kcal/mol, the bar used throughout this project

# Below this weight there is nothing left to extrapolate and nothing needs
# extrapolating: c * w is already far under chemical accuracy for any physical c.
# The right output then is the computed energy with c*w as a bound, NOT a failed fit.
W_CONVERGED = 1e-8

# Conservative upper bound on the proportionality constant c in E ~= E_0 + c*w,
# used only to turn a small w into an error bar without relying on the fit.
# Measured on localised H12/STO-6G at r=1.4 A: c ~= 8.45 Ha per unit weight.
# 100 is a deliberate order-of-magnitude margin over that.
C_BOUND_HA_PER_W = 100.0

# Largest bond dimension it is worth planning a production run around, for an
# active space of a few tens of orbitals. Used only to decide whether a fitted
# decay law can ever reach the linear regime.
M_FEASIBLE = 4000


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


def sweep_convergence(bond_dims, energies, root=0, target_mha=0.1, tail=4):
    """Are the sweeps converged at each fixed M? Reads the PER-SWEEP arrays you
    already have, so it costs nothing to run after the fact.

    This exists because a negative fitted slope in extrapolate() has TWO causes
    that need OPPOSITE fixes, and the slope alone cannot tell them apart:

      (a) too few sweeps  -- the energy at the larger M simply has not come down
          yet. The per-sweep energy is still falling at the end of each M. Fix:
          more sweeps, same M, cold start. The drift here predicts how many.

      (b) different states -- the sweeps ARE flat at each M, but the optimiser
          converged to different local minima at different M, so the across-M
          sequence is non-monotonic between converged answers. More sweeps will
          not help at all. Fix: warm-start each M from the previous one.

    Pass the per-sweep bond_dims and energies straight from
    DMRGDriver.get_dmrg_results() (energies entries may be lists over roots).

    Returns per-M drift over the last `tail` sweeps, a geometric estimate of the
    sweeps needed to reach target_mha, and a verdict naming which cause fits.
    """
    import math

    def scalar(e):
        try:
            return float(e[root])
        except (TypeError, IndexError):
            return float(e)

    n = min(len(bond_dims), len(energies))
    by_m = {}
    for i in range(n):
        by_m.setdefault(int(bond_dims[i]), []).append(scalar(energies[i]))

    rows, flat = [], True
    for m in sorted(by_m):
        es = by_m[m]
        t = es[-tail:] if len(es) >= 2 else es
        steps = [abs(t[j + 1] - t[j]) for j in range(len(t) - 1)]
        drift = steps[-1] if steps else None
        ratio = None
        if len(steps) >= 2 and steps[-2] > 0:
            ratio = steps[-1] / steps[-2]
        need = None
        if drift is not None and drift * 1000 > target_mha:
            flat = False
            if ratio is not None and 0 < ratio < 1:
                need = math.ceil(math.log((target_mha / 1000) / drift) / math.log(ratio))
            else:
                need = -1          # not decaying: more sweeps are not converging
        rows.append({"bond_dim": m, "n_sweeps": len(es),
                     "last_energy": es[-1],
                     "last_step_mha": None if drift is None else drift * 1000,
                     "step_ratio": ratio,
                     "converged_at_this_M": drift is not None and drift * 1000 <= target_mha,
                     "extra_sweeps_needed": need})

    e_last = [r["last_energy"] for r in rows]
    monotone = e_last == sorted(e_last, reverse=True)
    worst = max((r["extra_sweeps_needed"] or 0) for r in rows) if rows else 0

    if flat and not monotone:
        verdict = ("DIFFERENT STATES. The sweeps are flat at every M (drift under "
                   "%.2f mHa) yet the energies are not monotone in M, so these are "
                   "converged answers to different local minima. More sweeps will not "
                   "help -- warm-start each bond dimension from the previous one."
                   % target_mha)
    elif not flat and worst < 0:
        verdict = ("NOT CONVERGING. The per-sweep step is not shrinking geometrically "
                   "at one or more bond dimensions, so more sweeps of the same kind "
                   "will not land it. Suspect the noise schedule or a bad initial "
                   "state rather than the sweep count.")
    elif not flat:
        verdict = ("TOO FEW SWEEPS. The per-sweep energy is still falling at the end "
                   "of at least one bond dimension; about %d more sweeps there reaches "
                   "%.2f mHa. Raise --n-sweeps at the SAME bond dimensions and cold "
                   "start again -- do not change M in the same run." % (worst, target_mha))
    else:
        verdict = ("CONVERGED AT EACH M and monotone in M: the ladder is diagnostic, "
                   "so a negative slope from extrapolate() would now be real and not "
                   "a sweep artifact.")

    return {"ok": True, "rows": rows, "all_flat": flat, "monotone_in_M": monotone,
            "target_mha": target_mha, "extra_sweeps_needed": worst, "verdict": verdict}


SMAX_BOND_FRACTION = 0.7   # above this fraction of ln(M), S_max is measuring M


def smax_bond_limited(s_max, bond_dim, frac=SMAX_BOND_FRACTION):
    """Is a reported S_max a property of the molecule, or of the bond dimension?

    An MPS with bond dimension M cannot carry more than ln(M) nats of entropy
    across any cut, so a measured S_max approaching ln(M) is reporting the
    truncation, not the state. Exactly the same disease as a verdict fixed by
    the active-space dimensions -- the answer is set by the configuration before
    the physics gets a say -- with the bond dimension as the parameter instead of
    the orbital count.

    Caught on a real run: S_max = 2.2448 at M = 16, where ln(16) = 2.7726, i.e.
    81% of the ceiling, tagged Class A against an S_HARD of 1.5 while the energy
    was still moving 534 mHa per doubling of M and sat 1665 mHa above the
    reference. Note also that S_max = 2.2448 exceeds ln(8) = 2.0794, so the same
    number could not have arisen at M = 8: the series itself shows the cap.
    """
    import math
    if bond_dim is None or bond_dim < 2 or s_max is None:
        return {"ok": False, "reason": "need S_max and a bond dimension >= 2"}
    ceiling = math.log(bond_dim)
    f = s_max / ceiling
    limited = f >= frac
    return {"ok": True, "s_max": s_max, "bond_dim": bond_dim,
            "ln_M_ceiling": ceiling, "fraction_of_ceiling": f,
            "bond_limited": limited,
            "verdict": ("S_max = %.4f is %.0f%% of the ln(M) = %.4f ceiling at M = %d. "
                        "This is a measurement of the bond dimension, not of the "
                        "molecule -- do not classify on it, and do not report it. "
                        "Raise M until S_max sits well below the ceiling, then read it."
                        % (s_max, 100 * f, ceiling, bond_dim) if limited else
                        "S_max = %.4f is %.0f%% of the ln(M) = %.4f ceiling at M = %d; "
                        "the bond dimension is not the binding constraint."
                        % (s_max, 100 * f, ceiling, bond_dim))}


def decay_law(bond_dims, discarded_weights, w_target=1e-5, m_feasible=M_FEASIBLE):
    """How does the discarded weight fall with bond dimension, and can the ladder
    ever reach the linear regime?

    This is a property of the ORBITAL BASIS AND ORDERING, not of the molecule, and
    it is the cheapest diagnostic available because it needs no reference energy.

    Two forms are fitted,

        log w = a + b*M        (exponential)
        log w = a + b*log M    (power law)

    but the VERDICT is taken from the implied bond dimension, not from which form
    wins. Measured on the same molecule, geometry and active space in two orbital
    bases, H12/STO-6G at r = 1.4 A, CAS(12,12):

        canonical  M = 12..48 : power law, exponent -0.663, R^2 0.941
                                -> M ~ 7.8e5 for w = 1e-5        UNUSABLE
        localised  M = 6..24  : power law, exponent -4.710, R^2 0.992
                                -> M ~ 16 for w = 1e-5           CONVERGED

    So the functional form does not separate the two cases -- both are power laws.
    The exponent does, by a factor of 7.1, and the bond dimension it implies does
    by a factor of 5e4. An earlier version of this function keyed on the form and
    would have been right by accident; it now keys on m_required.

    For reference, on the same pair of runs the quantities that did NOT separate
    them were the R^2 of the E-against-w fit (0.941 failing, 0.992 succeeding --
    both high) and the constancy of c = E_err/w (spread 1.45x in the FAILING run
    against 2.51x in the succeeding one, i.e. the wrong way round). The ones that
    did were the largest w on the ladder, the decades of w spanned, and this.

    Returns a dict; 'basis_suspect' True means fix the orbitals, not the ladder.
    """
    import math
    pts = sorted((int(m), float(w)) for m, w in zip(bond_dims, discarded_weights)
                 if w is not None and w > 0)
    if len(pts) < 3:
        return {"ok": False, "reason": "need at least 3 positive weights, got %d" % len(pts)}
    M = [p[0] for p in pts]
    lw = [math.log(p[1]) for p in pts]

    aE, bE, r2E = _ols(M, lw)
    aP, bP, r2P = _ols([math.log(m) for m in M], lw)
    if bE is None or bP is None:
        return {"ok": False, "reason": "bond dimensions or weights are degenerate"}

    lt = math.log(w_target)
    m_exp = (lt - aE) / bE if bE < 0 else float("inf")
    try:
        m_pow = math.exp((lt - aP) / bP) if bP < 0 else float("inf")
    except OverflowError:
        m_pow = float("inf")

    form = "power" if r2P > r2E else "exponential"
    m_req = m_pow if form == "power" else m_exp
    # Keyed on the IMPLIED BOND DIMENSION, not on the functional form. Settled on
    # the project's own two runs of the same molecule and active space: BOTH are
    # better fit by a power law (R^2 0.992 localised, 0.941 canonical), so the
    # form does not discriminate. What discriminates is the exponent -- -4.71
    # localised against -0.663 canonical, 7.1x steeper -- and therefore the bond
    # dimension each implies: M ~ 16 against M ~ 7.8e5 for w = 1e-5.
    suspect = m_req > m_feasible

    return {"ok": True,
            "bond_dims": M,
            "w": [p[1] for p in pts],
            "w_decades": math.log10(pts[0][1] / pts[-1][1]),
            "form": form,
            "r2_exponential": r2E, "r2_power": r2P, "r2_margin": r2P - r2E,
            "exp_rate_per_M": bE, "power_exponent": bP,
            "w_target": w_target,
            "m_required": m_req,
            "m_feasible": m_feasible,
            "basis_suspect": suspect,
            "note": ("the weight falls as M^%.3f, a power law; reaching w = %.0e would need "
                     "M ~ %.2g against a feasible ceiling of %d. No bond-dimension ladder fixes "
                     "this. The orbitals are delocalised: localise inside the active space "
                     "(split Pipek-Mezey, occupied among occupied and virtual among virtual, "
                     "which leaves the energy in the space exactly invariant) and reorder on "
                     "the exchange matrix, then re-probe."
                     % (bP, w_target, m_req, m_feasible)) if suspect else
                    ("the weight falls as exp(%.4f*M); w = %.0e is reached by M ~ %.0f, which is "
                     "within reach. The basis and ordering are adequate."
                     % (bE, w_target, m_exp)) if form == "exponential" else
                    ("the weight falls as M^%.3f, a power law, but w = %.0e is still reached by "
                     "M ~ %.2g. Usable, though localising the orbitals would shorten the ladder."
                     % (bP, w_target, m_req))}


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

    import math
    w_all = [w for _, _, w in pts]
    if max(w_all) < W_CONVERGED:
        # Not a failure, and not the same thing as the machine-epsilon artifact
        # below. If the largest weight on the whole ladder is already this small,
        # c*w is orders of magnitude under chemical accuracy for any physical c,
        # so the computed energy IS the answer and the honest output is a bound,
        # not an intercept. This is the expected case for a weakly correlated
        # closed-shell cluster: the risk there is too little truncation, not too
        # much. The one thing still worth checking is that the field really is a
        # discarded weight, and the tell for that is a non-physical slope.
        _, s_chk, _ = _ols(w_all, [e for _, e, _ in pts])
        bound = C_BOUND_HA_PER_W * max(w_all)
        w_ch = (math.log10(max(w_all) / min(w_all))
                if min(w_all) > 0 and max(w_all) > min(w_all) else 0.0)
        suspect_field = s_chk is not None and abs(s_chk) > 1e4
        return {"ok": True, "already_converged": True, "converged": True,
                "fitted_points": 0, "bond_dims_used": [p[0] for p in pts],
                "w_min": min(w_all), "w_max": max(w_all), "w_decades": w_ch,
                "e_best_computed_Ha": pts[-1][1], "bond_dim_best": pts[-1][0],
                "e_extrapolated_Ha": None, "slope_Ha_per_w": s_chk, "r2": None,
                "truncation_bound_Ha": bound, "truncation_bound_mHa": bound * 1000.0,
                "chem_acc_Ha": chem_acc,
                "verdict": ("CONVERGED WITHOUT EXTRAPOLATION: largest discarded weight on the "
                            "ladder is %.2e, so the truncation error is bounded by %.2e Ha "
                            "(c <= %.0f Ha per unit weight), %.0fx under chemical accuracy. "
                            "Report the computed energy at M = %d with this bound; an "
                            "extrapolation here would be fitting numerical noise."
                            % (max(w_all), bound, C_BOUND_HA_PER_W,
                               chem_acc / bound if bound > 0 else float('inf'), pts[-1][0])),
                "warnings": ([] if not suspect_field else
                             ["fitted slope is %.3e Ha per unit weight, far above the physical "
                              "range of order 1e0-1e2; the quantity passed in is probably not a "
                              "discarded weight" % s_chk])}

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
    if max(W) > MAX_W_LINEAR:
        # Caught from a real failed validation on a stretched H12 chain, where the
        # largest weight was 1.8e-2. The relation E ~= E0 + c*w is a LEADING-ORDER
        # statement; it needs w small. Published extrapolations work at w ~ 1e-5 to
        # 1e-7. At the percent level the higher-order terms dominate and the fitted
        # intercept is meaningless however good the R^2 looks.
        warn.append("largest discarded weight is %.2e, above the %.0e ceiling for the "
                    "linear regime -- the ladder has not reached the regime the method "
                    "assumes. Raise the bond dimensions, or move to a system where this "
                    "weight is reachable, rather than fitting here" % (max(W), MAX_W_LINEAR))
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

    # Basis-and-ordering guard. Fitted on the WHOLE ladder, not the fitted window,
    # because the decay law is what tells you whether a longer ladder would help
    # at all. A power-law decay means it would not.
    decay = decay_law([p[0] for p in pts], [p[2] for p in pts])
    if decay.get("ok") and decay["basis_suspect"]:
        warn.append("the discarded weight decays as a POWER LAW in the bond dimension "
                    "(R^2 %.3f against %.3f for exponential); %s"
                    % (decay["r2_power"], decay["r2_exponential"], decay["note"]))

    converged = (resid < chem_acc) and not warn
    return {"ok": True, "already_converged": False, "decay": decay,
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
