#!/usr/bin/env python3
"""
check_active_space.py -- blocking precondition before any DMRG/SHCI classification run.

The general rule it enforces:

    Refuse to run a classification whose verdict is already determined by the
    dimensions of the active space, independently of the chemistry.

An active space is a container. Before spending compute on measuring entanglement
inside it, check that it can hold any. Four ways it cannot:

  R1  no virtual orbitals            -> exactly one determinant, S_max == 0 by construction
  R2  virtual fraction too small     -> static correlation is occupied-virtual mixing;
                                        with a few percent virtual there is nowhere to mix
  R3  entropy ceiling below S_HARD   -> the symmetry-aware upper bound on the bipartition
                                        entropy is lower than the Class A trigger, so
                                        Class A is unreachable whatever the chemistry
  R4  FCI dimension below the floor  -> the space is exactly diagonalisable; run FCI and
                                        get the exact answer instead of a DMRG estimate
  R7  the boundary was imposed, not measured (see cut_quality.py) -- a cut inside a flat
                                        continuum of the selection spectrum needs a
                                        sensitivity ladder before its verdict is trusted

Every threshold is a declared parameter, printed with the verdict.

Usage
  check_active_space.py --nelecas 76 --ncas 38
  check_active_space.py --avas-json wt_avas.json --compare-json mutant_avas.json
  check_active_space.py --avas-json wt.json --compare-json mut.json \
      --spectrum-json wt_spectrum.json --compare-spectrum-json mut_spectrum.json
"""
import argparse, json, math, sys
from math import comb, log


def smax_ceiling(K, n_elec):
    """Rigorous ceiling on the bipartition von Neumann entropy of a closed-shell
    singlet CAS(n_elec, K):  S <= ln( sum_sectors min(dim_left, dim_right) ),
    the sum running over the (n_up, n_dn) sectors reachable at fixed particle
    number and Sz. Returns (best_cut, ceiling_in_nats).

    This is strictly tighter than k*ln4 whenever the space is not half filled:
    a nearly full space has few holes, so its blocks have little room for entropy.
    """
    Nu = n_elec // 2
    best = (None, -1.0)
    for k in range(1, K):
        R = K - k
        lo, hi = max(0, Nu - R), min(k, Nu)
        if lo > hi:
            continue
        tot = 0
        for nu in range(lo, hi + 1):
            cl, cr = comb(k, nu), comb(R, Nu - nu)
            for nd in range(lo, hi + 1):
                dl = cl * comb(k, nd)
                dr = cr * comb(R, Nu - nd)
                tot += dl if dl < dr else dr
        if tot > 0:
            b = log(tot)
            if b > best[1]:
                best = (k, b)
    return best


def load(path):
    j = json.load(open(path))
    for ke, kc in (('nelecas', 'ncas'), ('n_elec_act', 'n_orb_act'),
                   ('active_electrons', 'active_orbitals')):
        if ke in j and kc in j:
            return int(j[ke]), int(j[kc])
    raise SystemExit("could not find an (nelecas, ncas) pair in %s -- keys present: %s"
                     % (path, sorted(j)[:20]))


def report(ne, no, a, label=""):
    occ = ne // 2
    vir = no - occ
    fails, warns = [], []
    print("=" * 74)
    print("ACTIVE SPACE %s CAS(%d, %d)   %d qubits under Jordan-Wigner" % (label, ne, no, 2 * no))
    if vir < 0:
        print("  MALFORMED: %d electrons cannot occupy %d orbitals" % (ne, no))
        return ['malformed']
    dets = comb(no, occ) ** 2
    print("  occupied %d | virtual %d | virtual fraction %.1f%% | filling %.3f"
          % (occ, vir, 100.0 * vir / no, ne / (2.0 * no)))
    print("  FCI dimension (closed-shell singlet) = %.3e determinants" % float(dets))

    if vir == 0:
        print("  [R1] FAIL  zero virtual orbitals -> exactly one determinant.")
        print("             S_max == 0.0 is a mathematical certainty here, not a measurement.")
        fails.append('R1 zero virtual orbitals')
    else:
        print("  [R1] ok    %d virtual orbitals present" % vir)

    vf = 100.0 * vir / no
    if vir > 0 and vf < a.min_vir_frac:
        print("  [R2] FAIL  virtual fraction %.1f%% below the %.1f%% floor -- the space has almost"
              % (vf, a.min_vir_frac))
        print("             no room for occupied-virtual mixing, so a low S_max is structural")
        fails.append('R2 virtual fraction %.1f%%' % vf)
    elif vir > 0:
        print("  [R2] ok    virtual fraction %.1f%% >= %.1f%%" % (vf, a.min_vir_frac))

    if vir == 0:
        ceil = 0.0
        print("  [R3] FAIL  entropy ceiling 0.000 nats (single determinant)")
        fails.append('R3 ceiling 0')
    else:
        k, ceil = smax_ceiling(no, ne)
        naive = (no // 2) * math.log(4)
        print("  [R3] entropy ceiling %.3f nats at the %d|%d cut" % (ceil, k, no - k))
        print("       (k*ln4 would say %.3f -- overstated by %.2fx at this filling)"
              % (naive, naive / ceil if ceil else float('inf')))
        if ceil < a.s_hard:
            print("  [R3] FAIL  ceiling %.3f < S_HARD %.3f -> Class A is unreachable by construction"
                  % (ceil, a.s_hard))
            fails.append('R3 ceiling %.3f < S_HARD' % ceil)
        else:
            print("  [R3] ok    ceiling is %.0fx S_HARD; an S_max of %.2f would be %.1f%% of capacity"
                  % (ceil / a.s_hard, a.s_hard, 100.0 * a.s_hard / ceil))

    if dets < a.fci_floor:
        print("  [R4] FAIL  %.2e determinants is below the %.0e floor -- exactly diagonalisable."
              % (float(dets), a.fci_floor))
        print("             Run FCI and report the exact energy; a DMRG estimate here adds nothing,")
        print("             and by the framework's own rule the space is Class C.")
        fails.append('R4 FCI-tractable (%.2e dets)' % float(dets))
    else:
        print("  [R4] ok    %.2e determinants is past the exact-diagonalisation floor" % float(dets))
        if dets > a.fci_ceiling:
            print("  [R4] note  %.2e is very large; check the bond-dimension schedule is adequate"
                  % float(dets))
            warns.append('very large space')
    return fails


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--nelecas', type=int)
    p.add_argument('--ncas', type=int)
    p.add_argument('--avas-json')
    p.add_argument('--compare-json', help='the other side of a WT/mutant pair')
    p.add_argument('--s-hard', type=float, default=1.5)
    p.add_argument('--min-vir-frac', type=float, default=15.0, help='percent')
    p.add_argument('--fci-floor', type=float, default=1e9)
    p.add_argument('--fci-ceiling', type=float, default=1e40)
    p.add_argument('--spectrum-json', help='{"occupied": [...], "virtual": [...]} -- the '
                   'selection diagnostic for every orbital in the POOL, each list sorted '
                   'most-correlated first. For MP2 natural orbitals the occupied list is '
                   'the HOLE occupations 2-n, not the occupations.')
    p.add_argument('--compare-spectrum-json', help='the same for the other side of the pair')
    p.add_argument('--min-q', type=float, default=3.0,
                   help='cut-quality floor: gap at the cut over the local median spacing')
    p.add_argument('--sensitivity-verified', action='store_true',
                   help='assert that the +-2 orbital ladder has been run and the verdict was '
                        'identical across it. Records an override; do not pass it unprovoked.')
    a = p.parse_args()

    if a.avas_json:
        ne, no = load(a.avas_json)
    elif a.nelecas and a.ncas:
        ne, no = a.nelecas, a.ncas
    else:
        raise SystemExit('give --nelecas/--ncas or --avas-json')

    fails = report(ne, no, a, label="A:" if a.compare_json else "")

    if a.compare_json:
        ne2, no2 = load(a.compare_json)
        fails += ['B: ' + f for f in report(ne2, no2, a, label="B:")]
        print("=" * 74)
        if (ne, no) != (ne2, no2):
            print("  [R5] FAIL  the two sides are not the same space: CAS(%d,%d) vs CAS(%d,%d)."
                  % (ne, no, ne2, no2))
            print("             A WT/mutant comparison requires identical (nelecas, ncas).")
            fails.append('R5 mismatched spaces')
        else:
            print("  [R5] ok    both sides are CAS(%d,%d)" % (ne, no))
        print("  [R6] note  a RAW energy difference between a wild-type and a mutant cluster is")
        print("             dominated by the atomic-composition change of the substitution and is")
        print("             not a stability measurement. Only Leu/Ile cancels exactly.")

    # ---- R7: was the boundary of this space measured, or imposed? ----------------
    if a.spectrum_json:
        import json as _json, os as _os
        # this interpreter may run with a safe path, so sys.path[0] is not the
        # script directory; put the sibling module within reach explicitly
        _here = _os.path.dirname(_os.path.abspath(__file__))
        if _here not in sys.path:
            sys.path.insert(0, _here)
        try:
            import cut_quality as _cq
        except ImportError:
            print("  [R7] skip  cut_quality.py not importable; cannot judge the cut")
            _cq = None
        if _cq is not None:
            n_occ, n_vir = ne // 2, no - ne // 2
            sides = [("A", a.spectrum_json)]
            if a.compare_spectrum_json:
                sides.append(("B", a.compare_spectrum_json))
            qs = []
            for tag, path in sides:
                sp = _json.load(open(path))
                qo = _cq.cut_quality(sp["occupied"], n_occ)
                qv = _cq.cut_quality(sp["virtual"], n_vir)
                qs += [qo, qv]
                for blk, k, q in (("occ", n_occ, qo), ("vir", n_vir, qv)):
                    txt = "n/a (pool too small to judge)" if q is None else "%.2f" % q
                    flag = "ok   " if (q is not None and q >= a.min_q) else "LOW  "
                    print("  [R7] %s %s:%s cut at %d -> q = %s" % (flag, tag, blk, k, txt))
            weak = [q for q in qs if q is None or q < a.min_q]
            if weak:
                plan = _cq.sensitivity_plan(n_occ, n_vir,
                                            min([q for q in qs[0::2] if q is not None] or [None]),
                                            min([q for q in qs[1::2] if q is not None] or [None]),
                                            min_q=a.min_q)
                print("             %d of %d cuts sit inside a continuum: the boundary was set by a"
                      % (len(weak), len(qs)))
                print("             parameter, not by the spectrum. That is permitted ONLY with the")
                print("             sensitivity ladder, because a cut in a flat region is harmless")
                print("             exactly when the orbitals either side of it are equivalent --")
                print("             which is testable, not arguable. Required runs, beside this one:")
                for cfg in plan["required_runs"]:
                    print("                CAS(%d,%d)  n_occ=%d n_vir=%d"
                          % (2 * cfg["n_occ"], cfg["n_occ"] + cfg["n_vir"],
                             cfg["n_occ"], cfg["n_vir"]))
                if a.sensitivity_verified:
                    print("             --sensitivity-verified given: override RECORDED, verdict may")
                    print("             be reported together with the ladder it was stable across.")
                else:
                    fails.append('R7 cut inside a continuum and no sensitivity ladder recorded')

    print("=" * 74)
    if fails:
        print("REFUSED -- this run's verdict is determined before it starts:")
        for f in fails:
            print("   * %s" % f)
        sys.exit(2)
    print("PRECONDITION PASSED -- the space can hold the answer being asked of it.")
    sys.exit(0)


if __name__ == '__main__':
    main()
