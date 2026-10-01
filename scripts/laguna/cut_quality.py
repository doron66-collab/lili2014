#!/usr/bin/env python3
"""
cut_quality.py -- is an active-space boundary a measurement or a parameter?

The problem this answers. Selecting an active space means cutting a sorted list
of orbital diagnostics (AVAS singular values, MP2 natural-orbital occupations)
at some rank. If the list has a real gap at that rank, the cut is a property of
the molecule. If the list is locally flat there, the cut is a property of
whatever set the rank -- a threshold, a target size, a fraction cap -- and any
verdict downstream inherits that arbitrariness.

The two cases are told apart by one dimensionless number: the gap at the cut
divided by the typical adjacent spacing around it.

    q(k) = [x_k - x_{k+1}] / median( adjacent spacings in a window around k )

q >> 1 is a cut the spectrum itself chose. q <~ 1 is a cut inside a continuum.

Measured on the project's own C275F run, virtual block, WT side:
    rank 16 (the natural gap)        gap 0.0168   q = 14.5
    rank 26 (where the run cut)      gap 0.0002   q =  0.17
a factor of 84. The rank-26 cut was not chosen by MP2 and not chosen by the gap
structure: n_vir = floor(0.55 * 48) = 26, i.e. --max-vir-frac set it exactly.

Two things follow, and they are separate:

1. The SIZE must stay common to both sides of a comparison. Different (n_elec,
   n_orb) on the two sides makes the energy difference and the entanglement
   measure incomparable; that is precondition R5 in check_active_space.py and it
   is not negotiable for an instrument whose output is a comparison. So the gap
   structure should CHOOSE the common size, not be overridden by it and not be
   applied per side.

2. When no common gap exists, a better-sounding justification does not fix
   anything. What fixes it is showing the verdict does not depend on the cut:
   re-run at +-2 orbitals per block and require the classification to be stable.
   An arbitrary cut inside a flat continuum is harmless exactly when the
   orbitals either side of it are equivalent -- which is a testable claim, not
   an argument.
"""

MIN_Q = 3.0          # below this, a cut is not distinguishable from the continuum
DEFAULT_WINDOW = 5   # ranks either side of k used to estimate the local spacing


def _median(v):
    s = sorted(v)
    n = len(s)
    if n == 0:
        return None
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def adjacent_gaps(x):
    """Spacings of a list sorted in decreasing order. gaps[i] = x[i] - x[i+1]."""
    return [x[i] - x[i + 1] for i in range(len(x) - 1)]


def cut_quality(x, k, window=DEFAULT_WINDOW):
    """Quality of a cut that keeps the first k entries of the sorted list x.

    k is a COUNT, so the relevant gap is between index k-1 and index k.
    Returns None when there is not enough spectrum around k to judge.
    """
    if not (1 <= k < len(x)):
        return None
    g = adjacent_gaps(x)
    i = k - 1                                  # the gap being used as the cut
    lo, hi = max(0, i - window), min(len(g), i + window + 1)
    neighbours = [g[j] for j in range(lo, hi) if j != i]
    if not neighbours:
        return None
    scale = _median(neighbours)
    if scale is None or scale <= 0:
        return float("inf") if g[i] > 0 else 0.0
    return g[i] / scale


def scan_cuts(x, min_q=MIN_Q, window=DEFAULT_WINDOW, k_min=1, k_max=None):
    """Every rank whose cut quality clears min_q, best first."""
    k_max = k_max or len(x) - 1
    out = []
    for k in range(max(1, k_min), min(k_max, len(x) - 1) + 1):
        q = cut_quality(x, k, window)
        if q is not None and q >= min_q:
            out.append({"k": k, "q": q, "gap": x[k - 1] - x[k],
                        "last_kept": x[k - 1], "first_dropped": x[k]})
    return sorted(out, key=lambda d: -d["q"])


def common_cuts(x_a, x_b, tol=1, min_q=MIN_Q, window=DEFAULT_WINDOW, k_max=None):
    """Ranks usable as ONE shared cut for both sides of the comparison.

    tol only widens the SEARCH: the mutation moves the spectrum, so a gap at
    rank 17 on one side and rank 18 on the other is the same chemical feature
    and both ranks become candidates. But the run uses a single k on both sides
    -- R5 requires identical (nelecas, ncas) -- so every candidate k is then
    scored by q evaluated AT THAT SAME k on BOTH sides.

    This matters, and an earlier version of this function got it wrong. If A has
    its gap at 17 (q = 6.4) and B at 18 (q = 21.4), reporting min(6.4, 21.4) =
    6.4 for the shared cut k = 17 is not the quality of that cut: on side B,
    k = 17 is one rank short of B's gap, so the number that governs is q_B(17),
    which was never computed. A cut is only as principled as its weaker side
    AT THE RANK ACTUALLY USED.

    Returns candidates best-first on q_shared_min, each carrying q at the shared
    k for both sides and, for reference, where each side's own gap sits.
    """
    A = {d["k"]: d for d in scan_cuts(x_a, min_q, window, k_max=k_max)}
    B = {d["k"]: d for d in scan_cuts(x_b, min_q, window, k_max=k_max)}
    cand = set()
    for ka in A:
        for kb in B:
            if abs(ka - kb) <= tol:
                cand.add(ka)
                cand.add(kb)

    out = []
    for k in sorted(cand):
        qa = cut_quality(x_a, k, window)
        qb = cut_quality(x_b, k, window)
        if qa is None or qb is None:
            continue
        out.append({"k": k,
                    "q_a_at_k": qa, "q_b_at_k": qb,
                    "q_shared_min": min(qa, qb),
                    "on_own_gap_a": k in A, "on_own_gap_b": k in B,
                    "q_a_own_gap": A[k]["q"] if k in A else None,
                    "q_b_own_gap": B[k]["q"] if k in B else None,
                    "both_clear_at_k": qa >= min_q and qb >= min_q,
                    # kept for backwards compatibility with earlier callers
                    "q_min": min(qa, qb)})
    return sorted(out, key=lambda d: -d["q_shared_min"])


def select_space(occ_wt, vir_wt, occ_mut, vir_mut,
                 min_vir_frac=0.25, min_det=1e9, max_total=None,
                 min_q=MIN_Q, window=DEFAULT_WINDOW):
    """Choose ONE (n_occ, n_vir) for both sides from the gap structure.

    occ_* are the occupied-block diagnostics sorted so that the most correlated
    comes first (for MP2 natural orbitals that means sorting the hole
    occupations 2 - n in decreasing order, NOT the occupations themselves);
    vir_* are the virtual-block occupations in decreasing order.

    Returns the recommendation plus the quality of every cut it considered, so
    the record can state whether the size was measured or imposed.
    """
    import math
    co = common_cuts(occ_wt, occ_mut, min_q=min_q, window=window)
    cv = common_cuts(vir_wt, vir_mut, min_q=min_q, window=window)

    cands = []
    for do in co:
        for dv in cv:
            n_occ, n_vir = do["k"], dv["k"]
            tot = n_occ + n_vir
            frac = n_vir / tot
            dim = math.comb(tot, n_occ) ** 2
            # every q here is evaluated at the SHARED rank on both sides
            q_worst = min(do["q_shared_min"], dv["q_shared_min"])
            ok = (frac >= min_vir_frac and dim >= min_det
                  and (max_total is None or tot <= max_total)
                  and do["both_clear_at_k"] and dv["both_clear_at_k"])
            cands.append({"n_occ": n_occ, "n_vir": n_vir, "n_elec": 2 * n_occ,
                          "total": tot, "vir_frac": frac, "fci_dim": dim,
                          "q_occ_wt": do["q_a_at_k"], "q_occ_mut": do["q_b_at_k"],
                          "q_vir_wt": dv["q_a_at_k"], "q_vir_mut": dv["q_b_at_k"],
                          "q_occ": do["q_shared_min"], "q_vir": dv["q_shared_min"],
                          "q_worst": q_worst,
                          "both_sides_clear": do["both_clear_at_k"] and dv["both_clear_at_k"],
                          "admissible": ok})
    adm = [c for c in cands if c["admissible"]]
    best = max(adm, key=lambda c: c["q_worst"]) if adm else None
    return {"recommended": best,
            "candidates": sorted(cands, key=lambda c: -c["q_worst"]),
            "occupied_common_gaps": co, "virtual_common_gaps": cv,
            "verdict": ("size chosen by the spectrum: both blocks cut on a gap present "
                        "on both sides (worst q = %.1f)" % best["q_worst"]) if best else
                       ("NO common gap clears q = %.1f in both blocks on both sides. The size "
                        "cannot be read off the spectrum here; fix it on computational grounds, "
                        "declare it imposed rather than measured, and report the sensitivity "
                        "ladder below -- the verdict is only reportable if it is stable."
                        % min_q)}


def sensitivity_plan(n_occ, n_vir, q_occ, q_vir, min_q=MIN_Q, step=2):
    """Which extra runs are required before a verdict at this cut is reportable.

    A cut whose quality clears min_q needs nothing: the spectrum separated the
    kept orbitals from the dropped ones. A cut inside a continuum needs proof
    that the verdict does not move when the cut does -- and that proof is cheap,
    because orbitals in a flat region are nearly equivalent by construction.
    """
    runs, why = [], []
    for name, n, q in (("n_occ", n_occ, q_occ), ("n_vir", n_vir, q_vir)):
        if q is None or q < min_q:
            why.append("%s cut has q = %s, below %.1f -- inside a continuum"
                       % (name, "n/a" if q is None else "%.2f" % q, min_q))
            for d in (-step, +step):
                cfg = {"n_occ": n_occ, "n_vir": n_vir}
                cfg[name] = n + d
                if cfg["n_occ"] > 0 and cfg["n_vir"] > 0:
                    runs.append(cfg)
    return {"required_runs": runs, "reasons": why,
            "rule": ("report the classification only if it is identical across all "
                     "%d runs; if it moves, the verdict is an artifact of the cut and "
                     "must be withheld" % (len(runs) + 1)) if runs else
                    "both cuts sit on real gaps; no sensitivity ladder required"}
