#!/usr/bin/env python3
"""
avas_mp2_select.py -- two-criterion active-space selection: AVAS (site relevance)
AND MP2 natural-orbital occupation (correlation relevance), replacing the
count-only cut in avas_by_count.py.

Why this exists (2026-09-30, Claude Science's review of the CAS(64,48) proposal):

avas_by_count.py cuts the AVAS-ranked occupied/virtual lists by COUNT instead
of by threshold -- correct fix for AVAS applying one threshold to two
differently-scaled sigma distributions (see that script's docstring). But on
the TP53_C275F cluster, Science pointed out a SECOND problem verifiable
directly from the existing threshold sweep: within the occupied block, sigma
saturates almost immediately (38 orbitals above sigma=0.95, a band only 0.05
wide -- average neighbor gap ~0.0013) and there is essentially no discriminating
signal left in sigma near the top of the ranking. Worse: sigma_occ near 1 means
an orbital is nearly PURE target-atom character (a lone pair or a pure bond) --
exactly the LEAST correlated kind of orbital. Cutting top-N by sigma alone
therefore systematically keeps the least-correlated occupied orbitals and
discards the HOMO-region, mixed-character orbitals where correlation actually
lives -- biasing any DMRG run on that space toward a low, uninformative S_max.

The fix (Science's proposed two-stage method, matching the original AVAS paper's
own point that AVAS measures "site relevance," not "correlation relevance"):

  1. AVAS(threshold=0.1) still defines the SITE-relevant subspace (unchanged --
     this is what keeps the >=90% site-coverage gate satisfied).
  2. Within that subspace, run density-fitted MP2 (cheap: the active pool is
     ~100-200 orbitals, not the full ~1000-function basis) and diagonalize its
     one-particle density matrix SEPARATELY in the occupied and virtual blocks
     -- giving MP2 natural-orbital (NO) occupations, which measure correlation
     directly (an occupation far from 2 or 0 IS static/dynamic correlation).
  3. For each MP2 natural orbital (a rotation of the AVAS active orbitals),
     recompute its AVAS/site-relevance score as the SAME projector expectation
     value AVAS itself diagonalizes (c^T @ Sa @ c, where Sa is the target-AO
     overlap projector) -- not the original AVAS sigma, which belonged to a
     now-superseded basis.
  4. Select orbitals passing BOTH a site-relevance floor and a correlation
     floor (occupation deviates from 2 or 0 by more than a declared cutoff),
     and record BOTH numbers per selected orbital -- so the record carries WHY
     each orbital is there, not just that AVAS or a count put it there.

Verified locally before writing the two-criterion machinery (see commit): on
N2/6-31g, this script's projector-score formula reproduces plain avas.avas's
own occ_weights/vir_weights EXACTLY (to float precision) when no MP2 rotation
is applied -- confirming the score formula itself is correct, independent of
whether it agrees with sigma after rotation.

Usage:
    python avas_mp2_select.py --xyz cluster.xyz --charge -1 --spin 0 \
        --basis 6-31g --ao-set "S 3p, N 2p, O 2p" \
        --site-threshold 0.1 --total-active 48 \
        --occ-dev-cutoff 0.02 --site-score-floor 0.05 \
        --out cluster_avas_mp2.json

Prints the sigma spectrum and the NO-occupation-deviation spectrum side by
side for the active pool BEFORE selecting anything -- this diagnostic alone
(per Science's suggestion) settles whether the two-criterion selection changes
anything for a given cluster: if sigma is flat where NO-occupation is not,
the count-only cut in avas_by_count.py was picking the wrong orbitals.
"""
import argparse
import json
import sys
from functools import reduce

import numpy
import scipy.linalg


def _site_projector(mf, aolabels, minao='minao', ncore=0):
    """The exact same projector avas.AVAS diagonalizes (Sa in the AVAS paper's
    notation), returned as a matrix acting on coordinates in the
    mf.mo_coeff[:, ncore:] basis -- i.e. sa[i,j] with i,j indexing that slice.
    """
    from pyscf import gto
    mol = mf.mol
    ovlp = mol.intor_symmetric('int1e_ovlp')
    pmol = mol.copy()
    pmol.atom = mol._atom
    pmol.unit = 'B'
    pmol.symmetry = False
    pmol.basis = minao
    pmol.build(False, False)
    aos = [s.strip() for s in aolabels.split(",")] if isinstance(aolabels, str) else aolabels
    baslst = pmol.search_ao_label(aos)
    s2 = pmol.intor_symmetric('int1e_ovlp')[baslst][:, baslst]
    s21 = gto.intor_cross('int1e_ovlp', pmol, mol)[baslst]
    s21 = numpy.dot(s21, mf.mo_coeff[:, ncore:])
    return s21.T.dot(scipy.linalg.solve(s2, s21, assume_a='pos'))


def two_criterion_select(mf, aolabels, site_threshold=0.1, occ_dev_cutoff=0.02,
                          site_score_floor=0.05, total_active=None,
                          min_vir_frac=0.20, max_vir_frac=0.55, ncore=0):
    """Returns a dict with the diagnostic spectra and the selected (n_occ, n_vir)
    orbital-index lists (indices into the assembled active-pool coordinate
    space -- caller reassembles mo_coeff from these, see build_mo()).
    """
    from pyscf import mp
    from pyscf.mcscf import avas

    nocc_total = numpy.count_nonzero(mf.mo_occ)
    nmo = mf.mo_coeff.shape[1]

    av = avas.AVAS(mf, aolabels, threshold=site_threshold, ncore=ncore)
    ncas0, nelec0, mo0 = av.kernel()
    ncas0_occ = nelec0 // 2
    ncas0_vir = ncas0 - ncas0_occ
    n_discarded_occ = (nocc_total - ncore) - ncas0_occ

    # --- MP2 within the AVAS(site_threshold) active pool only (everything else frozen)
    mf2 = mf.copy()
    mf2.mo_coeff = mo0
    mo_occ_new = numpy.zeros(nmo)
    mo_occ_new[:nocc_total] = 2.0
    mf2.mo_occ = mo_occ_new
    frozen = list(range(0, ncore + n_discarded_occ)) + list(range(ncore + n_discarded_occ + ncas0, nmo))
    # Printed BEFORE the MP2 call (not after, as a first version of this had it) --
    # if get_nocc()'s own `assert nocc > 0` fires, this is the only place the
    # actual counts get seen at all, instead of a bare AssertionError.
    n_occ_active_unfrozen = sum(1 for i in range(ncore + n_discarded_occ, ncore + n_discarded_occ + ncas0_occ)
                                if i not in frozen)
    print(f"[avas-mp2] diagnostic before MP2: nmo={nmo} nocc_total={nocc_total} ncore={ncore} "
          f"ncas0={ncas0} ncas0_occ={ncas0_occ} ncas0_vir={ncas0_vir} "
          f"n_discarded_occ={n_discarded_occ} len(frozen)={len(frozen)} "
          f"unfrozen_active_occ={n_occ_active_unfrozen}", flush=True)
    if n_occ_active_unfrozen <= 0:
        sys.exit(f"[avas-mp2] REFUSING: the active occupied block has {n_occ_active_unfrozen} "
                 f"unfrozen orbitals after the freeze list -- MP2 cannot run. This means "
                 f"ncas0_occ={ncas0_occ} itself is <= 0 (AVAS at threshold={site_threshold} "
                 f"selected no occupied orbitals for this AO set) or n_discarded_occ is "
                 f"negative/wrong. Check the numbers printed above before re-running.")
    pt = mp.MP2(mf2, frozen=frozen).density_fit().run()
    dm1 = pt.make_rdm1()

    lo, hi_occ = ncore + n_discarded_occ, ncore + n_discarded_occ + ncas0_occ
    hi_vir = hi_occ + ncas0_vir
    occ_block = dm1[lo:hi_occ, lo:hi_occ]
    vir_block = dm1[hi_occ:hi_vir, hi_occ:hi_vir]
    occ_no_occ, occ_no_u = numpy.linalg.eigh(occ_block)     # ascending
    vir_no_occ, vir_no_u = numpy.linalg.eigh(vir_block)     # ascending

    # site-relevance score of each NEW natural orbital = projector expectation
    # value in the SAME coordinate space avas.AVAS itself used (mo_coeff[:,ncore:])
    sa = _site_projector(mf, aolabels, ncore=ncore)
    # av.mo_coeff columns [ncore+n_discarded_occ : ncore+n_discarded_occ+ncas0] are
    # the AVAS active pool expressed in AO basis; recover their coordinates in the
    # original mf.mo_coeff[:,ncore:] basis via the overlap (they are orthonormal
    # combinations of it, so this is just an inner product with the ORIGINAL mo set).
    ovlp = mf.mol.intor_symmetric('int1e_ovlp')
    active_pool_ao = mo0[:, lo:hi_vir]              # AO basis, occ-then-vir active pool
    orig_mo_ncore = mf.mo_coeff[:, ncore:]           # AO basis
    coord = orig_mo_ncore.T.dot(ovlp).dot(active_pool_ao)   # coords of active pool in orig-mo space
    coord_occ, coord_vir = coord[:, :ncas0_occ], coord[:, ncas0_occ:]
    no_coord_occ = coord_occ.dot(occ_no_u)           # NEW natural orbitals, same coord space
    no_coord_vir = coord_vir.dot(vir_no_u)
    site_score_occ = numpy.einsum('ij,jk,ik->i', no_coord_occ.T, sa, no_coord_occ.T)
    site_score_vir = numpy.einsum('ij,jk,ik->i', no_coord_vir.T, sa, no_coord_vir.T)

    occ_dev = 2.0 - occ_no_occ     # 0 = uncorrelated, larger = more correlated
    vir_dev = vir_no_occ           # 0 = uncorrelated, larger = more correlated

    print("[avas-mp2] active pool from AVAS(threshold=%.2f): occ=%d vir=%d"
          % (site_threshold, ncas0_occ, ncas0_vir))
    print("[avas-mp2] occ block -- NO occupation (desc) | deviation from 2 | site score:")
    order_o = numpy.argsort(-occ_dev)
    for i in order_o:
        print("    occ_no=%.4f  dev=%.4f  site_score=%.4f" % (occ_no_occ[i], occ_dev[i], site_score_occ[i]))
    print("[avas-mp2] vir block -- NO occupation (desc) | deviation from 0 | site score:")
    order_v = numpy.argsort(-vir_dev)
    for i in order_v:
        print("    vir_no=%.4f  dev=%.4f  site_score=%.4f" % (vir_no_occ[i], vir_dev[i], site_score_vir[i]))

    passes_occ = (occ_dev > occ_dev_cutoff) & (site_score_occ > site_score_floor)
    passes_vir = (vir_dev > occ_dev_cutoff) & (site_score_vir > site_score_floor)
    cand_occ = numpy.where(passes_occ)[0]
    cand_vir = numpy.where(passes_vir)[0]
    print(f"[avas-mp2] pass both criteria: {len(cand_occ)} occ, {len(cand_vir)} vir "
          f"(occ_dev_cutoff={occ_dev_cutoff}, site_score_floor={site_score_floor})")

    if total_active is not None:
        # rank each pool by occ_dev/vir_dev (correlation strength) among those that
        # passed both criteria, then fill toward total_active within the declared
        # virtual-fraction band rather than a pre-fixed split.
        cand_occ = cand_occ[numpy.argsort(-occ_dev[cand_occ])]
        cand_vir = cand_vir[numpy.argsort(-vir_dev[cand_vir])]
        best = None
        for n_vir in range(min(len(cand_vir), total_active),
                            max(min(len(cand_vir), int(total_active * max_vir_frac)) - 1, -1), -1):
            n_occ = total_active - n_vir
            if n_occ < 0 or n_occ > len(cand_occ):
                continue
            vf = n_vir / total_active
            if min_vir_frac <= vf <= max_vir_frac:
                best = (n_occ, n_vir)
                break
        if best is None:
            sys.exit(f"[avas-mp2] REFUSING: no (n_occ,n_vir) summing to {total_active} orbitals "
                     f"with virtual fraction in [{min_vir_frac:.0%},{max_vir_frac:.0%}] survives "
                     f"both criteria (candidates: {len(cand_occ)} occ, {len(cand_vir)} vir). "
                     f"Widen --occ-dev-cutoff/--site-score-floor or --total-active.")
        n_occ, n_vir = best
        sel_occ, sel_vir = sorted(cand_occ[:n_occ]), sorted(cand_vir[:n_vir])
    else:
        sel_occ, sel_vir = sorted(cand_occ.tolist()), sorted(cand_vir.tolist())

    return dict(ncas0_occ=int(ncas0_occ), ncas0_vir=int(ncas0_vir),
                sel_occ=sel_occ, sel_vir=sel_vir,
                occ_no_u=occ_no_u, vir_no_u=vir_no_u,
                coord_occ=coord_occ, coord_vir=coord_vir,
                occ_dev=occ_dev.tolist(), vir_dev=vir_dev.tolist(),
                site_score_occ=site_score_occ.tolist(), site_score_vir=site_score_vir.tolist(),
                mo0=mo0, lo=lo, hi_occ=hi_occ, hi_vir=hi_vir, nmo=nmo, ncore=ncore,
                n_discarded_occ=n_discarded_occ)


def build_mo(mf, result):
    """Assemble the full mo_coeff with the active block replaced by the
    selected MP2 natural orbitals (occ-then-vir), everything else unchanged."""
    mo0, lo, hi_occ, hi_vir = result['mo0'], result['lo'], result['hi_occ'], result['hi_vir']
    no_occ_ao = mo0[:, lo:hi_occ].dot(result['occ_no_u'])[:, result['sel_occ']]
    no_vir_ao = mo0[:, hi_occ:hi_vir].dot(result['vir_no_u'])[:, result['sel_vir']]
    return numpy.hstack((mo0[:, :lo], no_occ_ao, no_vir_ao, mo0[:, hi_vir:]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--charge", type=int, required=True)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ao-set", required=True)
    ap.add_argument("--site-threshold", type=float, default=0.1)
    ap.add_argument("--occ-dev-cutoff", type=float, default=0.02,
                    help="minimum |occupation - formal| to count as correlated (default 0.02e)")
    ap.add_argument("--site-score-floor", type=float, default=0.05,
                    help="minimum projector score to count as site-relevant")
    ap.add_argument("--total-active", type=int, default=None,
                    help="if given, pick the best (n_occ,n_vir) summing to this within "
                         "--min-vir-frac/--max-vir-frac; otherwise report every orbital "
                         "passing both criteria and let the caller decide")
    ap.add_argument("--min-vir-frac", type=float, default=0.20)
    ap.add_argument("--max-vir-frac", type=float, default=0.55)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
    from run_gate2_avas import build_mf

    mf = build_mf(a.xyz, a.charge, a.spin, basis=a.basis)
    result = two_criterion_select(mf, a.ao_set, site_threshold=a.site_threshold,
                                   occ_dev_cutoff=a.occ_dev_cutoff,
                                   site_score_floor=a.site_score_floor,
                                   total_active=a.total_active,
                                   min_vir_frac=a.min_vir_frac, max_vir_frac=a.max_vir_frac)
    n_occ, n_vir = len(result['sel_occ']), len(result['sel_vir'])
    ncas, nelec = n_occ + n_vir, 2 * n_occ
    print(f"[avas-mp2] FINAL selection: CAS({nelec},{ncas})  n_occ={n_occ} n_vir={n_vir} "
          f"({100*n_vir/ncas:.1f}% virtual)")

    mo = build_mo(mf, result)
    out = dict(xyz=a.xyz, charge=a.charge, spin=a.spin, basis=a.basis, ao_set=a.ao_set,
               site_threshold=a.site_threshold, occ_dev_cutoff=a.occ_dev_cutoff,
               site_score_floor=a.site_score_floor,
               n_occ=n_occ, n_vir=n_vir, ncas=int(ncas), nelecas=int(nelec),
               qubits=int(2 * ncas), e_scf=float(mf.e_tot), scf_converged=bool(mf.converged),
               selected_occ_deviations=[result['occ_dev'][i] for i in result['sel_occ']],
               selected_occ_site_scores=[result['site_score_occ'][i] for i in result['sel_occ']],
               selected_vir_deviations=[result['vir_dev'][i] for i in result['sel_vir']],
               selected_vir_site_scores=[result['site_score_vir'][i] for i in result['sel_vir']])
    numpy.save(a.out.rsplit(".", 1)[0] + "_mo.npy", mo)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[avas-mp2] wrote {a.out} and {a.out.rsplit('.', 1)[0]}_mo.npy "
          f"(pass the .npy to solange_dmrg.py --load-orbitals --ncas {ncas} --nelecas {nelec})")


if __name__ == "__main__":
    main()
