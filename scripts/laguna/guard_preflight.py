#!/usr/bin/env python3
"""
guard_preflight.py -- a missing guard must be a HARD failure, never a skipped line.

Why this exists. A production classification run emitted CLASS A with this line
in its own log:

    [discarded-weight extrapolation] skipped (AttributeError: module
    'dmrg_extrapolate' has no attribute 'extrapolate')

The convergence criterion could not run, the exception was caught, the run
printed "skipped", and the verdict was emitted anyway. The provenance record is
complete and truthful: it says skipped. Nothing in it is false. And the class
that came out of that run has no convergence criterion behind it at all.

That is the governed-instrument failure mode in its purest form -- not an
adversary, not a falsified record, just a soft failure in a layer whose whole
purpose is to refuse. So the rule is structural:

    if the guard layer cannot run, no verdict may be emitted.

Call require_guards() before any classification, and classify_preflight() on the
energy that will ACTUALLY be classified -- which on that run was the ladder, not
the reference solve the embedding gate had been given.
"""
import importlib
import sys

REQUIRED = {
    "dmrg_extrapolate": ["extrapolate", "last_sweep_per_bond_dim", "window_stability",
                         "decay_law", "smax_bond_limited", "sweep_convergence",
                         "weight_explains_energy", "same_state_across_M"],
    "basis_check": ["fingerprint", "compare", "check_embedding"],
    "cut_quality": ["cut_quality", "scan_cuts", "common_cuts", "sensitivity_plan"],
}


def require_guards(required=None, verbose=True):
    """Import every guard module and verify every required symbol is present.

    Raises SystemExit(2) on any gap. Prints each module's resolved __file__ so a
    stale or shadowed copy on sys.path is visible -- that is the likeliest cause
    of a missing attribute, and it is invisible from the exception alone.
    """
    required = required or REQUIRED
    # this interpreter may run with a safe path, in which case sys.path[0] is not
    # the script directory and sibling guard modules are invisible. Put our own
    # directory within reach before deciding anything is missing -- otherwise the
    # preflight reports a gap that is really a path problem, which is a different
    # failure and sends the reader looking in the wrong place.
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    found, missing, broken = {}, [], []
    for mod, syms in required.items():
        try:
            m = importlib.import_module(mod)
        except Exception as e:                      # noqa: BLE001 - we report, not handle
            broken.append("%s: %s: %s" % (mod, type(e).__name__, e))
            continue
        path = getattr(m, "__file__", "<unknown>")
        # Byte size is reported because a file truncated to 0 bytes by a bad sync is
        # a VALID Python module: ast.parse("") raises nothing and import succeeds, so
        # a syntax check passes and the only symptom is an AttributeError at call time.
        try:
            nbytes = os.path.getsize(path)
        except OSError:
            nbytes = -1
        absent = [s for s in syms if not hasattr(m, s)]
        found[mod] = {"path": path, "bytes": nbytes, "n_syms": len(syms) - len(absent),
                      "missing": absent}
        if nbytes == 0:
            missing.append("%s (%s) is ZERO BYTES -- truncated, not merely stale" % (mod, path))
        elif absent:
            missing.append("%s (%s, %d bytes) lacks: %s"
                           % (mod, path, nbytes, ", ".join(absent)))

    if verbose:
        for mod, info in found.items():
            mark = "ok  " if not info["missing"] and info["bytes"] else "GAP "
            print("  [guards] %s %-20s %8d B  %2d/%-2d syms  %s"
                  % (mark, mod, info["bytes"], info["n_syms"],
                     len(REQUIRED[mod]), info["path"]))

    if broken or missing:
        print("\n  [guards] PREFLIGHT FAILED -- no verdict may be emitted:", file=sys.stderr)
        for line in broken + missing:
            print("     * %s" % line, file=sys.stderr)
        print("     A missing guard is not a skippable step. Most often this is a stale "
              "or shadowed copy earlier on sys.path than the one you edited -- compare "
              "the paths printed above against the files you changed.", file=sys.stderr)
        raise SystemExit(2)
    if verbose:
        print("  [guards] preflight passed: %d modules, all required symbols present"
              % len(found))
    return found


def classify_preflight(e_to_classify, e_scf, label="", ncore=None,
                       nelectron=None, nelecas=None, bond_dims=None,
                       energies=None, discarded_weights=None, s_max=None,
                       sweep_bond_dims=None, sweep_energies=None):
    """Run every applicable gate on the quantity that will ACTUALLY be classified.

    The gap this closes is specific and was not hypothetical: on the run above
    the embedding gate was given the reference DMRG-SCF solve, where it passed
    (293 mHa of correlation, variational bound holding), and the bond-dimension
    LADDER then ran separately and was never checked -- while every one of its
    seven points sat 180 to 711 mHa ABOVE the SCF energy, which the same gate
    would have refused outright.

    Pass the ladder you are about to read a verdict from, not the best run you
    have lying around.

    bond_dims/energies here are ONE (M, E) pair per bond dimension -- the LAST
    sweep's energy at each M -- because that is what weight_explains_energy()
    and the per-point embedding check need. sweep_convergence() needs the
    opposite shape (every sweep at every M, M repeated), so it is a DIFFERENT
    gap in the same function if given the collapsed per-M data instead: with
    one point per M there is no step between sweeps to measure, sweep_
    convergence() trivially reports every M "flat", and a real non-convergence
    like the real TP53_C275F M=75 (whose last-sweep step GREW, not shrank) is
    invisible. Pass the real per-sweep arrays via sweep_bond_dims/sweep_energies
    (straight from run_dmrg()'s own sweep_history, not energies/bond_dims) to
    get the check it was built for; omitting them falls back to the collapsed
    (uninformative but harmless) check rather than failing.
    """
    import basis_check
    import dmrg_extrapolate as dx

    report = {"label": label, "blocking": [], "advisory": []}

    emb = basis_check.check_embedding(e_to_classify, e_scf, ncore=ncore,
                                      nelectron=nelectron, nelecas=nelecas)
    report["embedding"] = emb
    if not emb["ok"]:
        report["blocking"] += emb["failures"]

    if bond_dims and energies:
        worst = None
        for m, e in zip(bond_dims, energies):
            ev = basis_check.check_embedding(e, e_scf)
            if not ev["ok"]:
                worst = (m, e) if worst is None or e > worst[1] else worst
        if worst is not None:
            report["blocking"].append(
                "at least one ladder point violates the variational bound; the worst is "
                "M = %d at E = %.8f, which is %.1f mHa above E_SCF = %.8f. A ladder cannot "
                "be read for convergence when its points are not valid energies."
                % (worst[0], worst[1], (worst[1] - e_scf) * 1000, e_scf))

    if bond_dims and energies and discarded_weights:
        wee = dx.weight_explains_energy(bond_dims, energies, discarded_weights)
        report["weight_explains_energy"] = wee
        if wee.get("ok") and not wee["usable_as_error_measure"]:
            report["blocking"].append(wee["verdict"])

    if sweep_bond_dims and sweep_energies:
        sc = dx.sweep_convergence(sweep_bond_dims, sweep_energies)
    elif bond_dims and energies:
        # Fallback: one point per M. sweep_convergence() cannot see a within-M
        # step with only one sample, so this trivially reports every M "flat" --
        # informational only, never a substitute for the real per-sweep check.
        sc = dx.sweep_convergence(bond_dims, energies)
    else:
        sc = None
    if sc is not None:
        report["sweep_note"] = sc.get("verdict")
        report["sweep_convergence"] = sc
        # The TWO most recent bond dimensions are what classify() actually reads
        # its dE from -- a non-decaying M anywhere else in the ladder (like M=75
        # when the interval actually used is (125, 250)) does not compromise the
        # number being classified. Block only when the compromised M is one of
        # the last two requested bond dimensions.
        not_decaying = set(sc.get("not_decaying_at") or [])
        if not_decaying and bond_dims:
            tail = set(sorted(bond_dims)[-2:])
            hit = not_decaying & tail
            if hit:
                report["blocking"].append(
                    "the bond dimension(s) %s this verdict's own dE is computed from did not "
                    "converge within their sweeps (sweep_convergence: %s) -- re-run with more "
                    "sweeps at that M before trusting the dE this classification is based on."
                    % (sorted(hit), sc["verdict"]))

    if s_max is not None and bond_dims:
        sb = dx.smax_bond_limited(s_max, max(bond_dims))
        report["smax"] = sb
        if sb.get("ok") and sb["bond_limited"]:
            report["blocking"].append(sb["verdict"])

    report["may_emit_verdict"] = not report["blocking"]
    report["verdict"] = ("preflight clear: a class may be emitted from this run"
                         if report["may_emit_verdict"] else
                         "NO VERDICT MAY BE EMITTED -- %d blocking finding(s)"
                         % len(report["blocking"]))
    return report


if __name__ == "__main__":
    require_guards()
    r = classify_preflight(
        -5332.94554423, -5333.12509349, label="TP53_C275F_LADDER_32_256 as run",
        ncore=329, nelectron=694, nelecas=36,
        bond_dims=[4, 8, 16, 32, 64, 128, 256],
        energies=[-5332.41382913, -5332.69508733, -5332.71219963, -5332.76615041,
                  -5332.84485829, -5332.93523928, -5332.94554423],
        s_max=1.7220954957471821)
    print("\n  %s" % r["verdict"])
    for b in r["blocking"]:
        print("     * %s" % b[:150])
