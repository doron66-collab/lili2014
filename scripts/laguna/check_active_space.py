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

Every threshold is a declared parameter, printed with the verdict.

Origin: written by Claude Science (2026-09-30) after the TP53_C275F CAS(76,38)
WT/mutant run -- occ=38, vir=0, S_max=0.0 was a mathematical certainty (the AVAS
threshold that was computationally tractable for DMRG-SCF happened to also carry
zero virtual orbitals). The check that would have caught this before the run took
hours of Laguna compute: occ = nelecas/2, vir = ncas - occ, both derivable from the
AVAS JSON that already existed. This script generalizes that check to four rules
and is meant to sit alongside the >=90% site-coverage gate, run before
solange_dmrg.py / solange_shci.py rather than after.

Usage
  check_active_space.py --nelecas 76 --ncas 38
  check_active_space.py --avas-json wt_avas.json --compare-json mutant_avas.json
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
