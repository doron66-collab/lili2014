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
                          random_seed=None,
                          site_score_floor=0.05, total_active=None,
                          min_vir_frac=0.20, max_vir_frac=0.55, ncore=0,
                          mp2_max_memory=128000, force_n_occ=None, force_n_vir=None):
    """Returns a dict with the diagnostic spectra and the selected (n_occ, n_vir)
    orbital-index lists (indices into the assembled active-pool coordinate
    space -- caller reassembles mo_coeff from these, see build_mo()).
    """
    from pyscf import mp
    from pyscf.mcscf import avas

    nocc_total = numpy.count_nonzero(mf.mo_occ)
    nmo = mf.mo_coeff.shape[1]

    # avas.AVAS/mol.search_ao_label want a LIST of labels, not a comma-joined
    # string -- verified locally: mol.search_ao_label('C 2p, N 2p') returns []
    # (no match at all) while mol.search_ao_label(['C 2p', 'N 2p']) matches
    # correctly. This is exactly why ncas0 came back 0 on the real cluster
    # (178 atoms, "C 2p, N 2p, O 2p, S 3p" passed as one string): AVAS silently
    # found zero target orbitals and never raised, since an empty AO list is
    # not itself an error to avas.avas. _site_projector already split the
    # string correctly; this call did not, and the two need to search the
    # SAME set of AOs for the score-consistency check above to mean anything.
    aos = [s.strip() for s in aolabels.split(",")] if isinstance(aolabels, str) else aolabels
    av = avas.AVAS(mf, aos, threshold=site_threshold, ncore=ncore)
    ncas0, nelec0, mo0 = av.kernel()
    ncas0_occ = nelec0 // 2
    ncas0_vir = ncas0 - ncas0_occ
    n_discarded_occ = (nocc_total - ncore) - ncas0_occ

    # --- MP2 within the AVAS(site_threshold) active pool only (everything else frozen)
    mf2 = mf.copy()
    # pyscf's MP2 checks t2's size (nocc^2 * nvir^2 doubles) against mp.max_memory
    # BEFORE running, not against the machine's real RAM -- inherited from mf's
    # own max_memory, which build_mf() never raises past pyscf's conservative
    # default (~4000 MB). Hit live on the real C275F cluster: AVAS(0.1) selected
    # occ=241/vir=139, so t2 alone is 241^2*139^2*8 bytes =~ 9.0 GB, well past
    # the default -- "Insufficient memory for holding t2 incore" even though
    # Laguna's largemem nodes have far more RAM than that. with_t2=False is NOT
    # the fix: DFMP2.make_rdm1() asserts self.t2 is not None (checked directly
    # in the installed pyscf source), so the 1-RDM this function needs requires
    # t2 to actually be held. Raise the ceiling instead of avoiding the array.
    mf2.max_memory = mp2_max_memory
    mf2.mol.max_memory = mp2_max_memory
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

    if force_n_occ is not None or force_n_vir is not None:
        # Bypass the occ_dev_cutoff/site_score_floor/total_active machinery
        # entirely: take the top-N orbitals by correlation rank (occ_dev /
        # vir_dev descending) straight from the FULL AVAS pool, exactly the
        # ranks cut_quality.select_space() reasons about when it reads the
        # --spectrum-out file (which is this same occ_dev/vir_dev, sorted
        # descending, over the WHOLE pool -- not the "passes both criteria"
        # subset). Using the site-score/deviation floors here would select a
        # DIFFERENT, smaller candidate pool than the one the gap was measured
        # on, silently reintroducing the same "cut inside a continuum"
        # failure cut_quality.py exists to catch -- one layer higher.
        if force_n_occ is None or force_n_vir is None:
            sys.exit("[avas-mp2] REFUSING: --force-n-occ and --force-n-vir must both be given "
                     "together (a spectrum-measured cut fixes both block sizes at once).")
        if not (0 < force_n_occ <= len(occ_dev)) or not (0 < force_n_vir <= len(vir_dev)):
            sys.exit(f"[avas-mp2] REFUSING: --force-n-occ={force_n_occ} / --force-n-vir={force_n_vir} "
                     f"out of range for pool sizes occ={len(occ_dev)} vir={len(vir_dev)}.")
        if random_seed is not None:
            # Scramble control (Claude Science, 2026-10-03): isolate SELECTION
            # pressure from SIZE/FILLING. Same pool, same force_n_occ/force_n_vir
            # count, but drawn uniformly at random instead of ranked by occ_dev --
            # if a random subset of this size from this same pool reproduces the
            # ranked subset's S_max, site_score was never doing any discriminating
            # work and the result is a property of (molecule, active-space size),
            # not of the site. Sampling the FULL candidate index range (not just
            # the "passes both criteria" subset), matching what the ranked branch
            # draws from (order_o/order_v span the whole pool).
            rng = numpy.random.default_rng(random_seed)
            sel_occ = sorted(rng.choice(len(occ_dev), size=force_n_occ, replace=False).tolist())
            sel_vir = sorted(rng.choice(len(vir_dev), size=force_n_vir, replace=False).tolist())
            print(f"[avas-mp2] SCRAMBLE selection (seed={random_seed}): n_occ={force_n_occ} "
                  f"n_vir={force_n_vir} drawn UNIFORMLY AT RANDOM from the full pool "
                  f"(occ_dev/site_score ignored by construction)")
        else:
            sel_occ = sorted(order_o[:force_n_occ].tolist())
            sel_vir = sorted(order_v[:force_n_vir].tolist())
            print(f"[avas-mp2] FORCED selection (spectrum-measured, bypassing occ_dev_cutoff/"
                  f"site_score_floor/total_active): n_occ={force_n_occ} n_vir={force_n_vir} "
                  f"ranked purely by occ_dev/vir_dev over the full pool")
    elif total_active is not None:
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
    """Assemble the FULL-WIDTH mo_coeff (same column count as mo0/nmo -- nothing
    dropped) with the active block holding only the SELECTED natural orbitals,
    and every AVAS-pool candidate that correlation-ranking did NOT select kept
    as an ordinary core (if occupied) or virtual (if virtual) column instead of
    being discarded.

    Found the hard way (2026-10-02): the previous version kept only the
    selected columns -- numpy.hstack((mo0[:, :lo], no_occ_ao[sel_occ],
    no_vir_ao[sel_vir], mo0[:, hi_vir:])) -- silently shrinking a 977-column
    system to 631. That breaks pyscf's own bookkeeping: CASCI/CASSCF compute
    ncore = (mol.nelectron - nelecas) // 2 from the MOLECULE's real electron
    count (694 here), not from whatever this function decided to keep, so
    get_h1eff/get_h2eff sliced mo_coeff[:, 329:363] against a 631-column array
    whose real active columns sat at 106:140 -- entirely the wrong orbitals,
    caught only via a variational-bound check (E_CASCI came out 623 Ha ABOVE
    E_SCF, impossible for a split-rotated active space that still contains the
    HF determinant). The 223 unselected occupied candidates and 123 unselected
    virtual candidates are REAL molecular orbitals of this system; dropping
    them silently changes the implied electron count, it does not shrink the
    problem.

    Returns (mo, ncore) -- ncore is the ACTUAL number of core columns in the
    returned array (lo + unselected-occ count), which the caller must record
    (as active_start_col) and the caller downstream must pass to pyscf AS
    ncore, not treat as a free offset. The split (occ-among-occ, vir-among-vir)
    rotation keeps the HF determinant exactly invariant regardless of where the
    unselected pool columns are placed within the core/virtual groups, so their
    internal order does not matter -- only the group (core vs. virtual) does.
    """
    mo0, lo, hi_occ, hi_vir = result['mo0'], result['lo'], result['hi_occ'], result['hi_vir']
    occ_no_u, vir_no_u = result['occ_no_u'], result['vir_no_u']
    sel_occ, sel_vir = set(result['sel_occ']), set(result['sel_vir'])

    no_occ_ao_full = mo0[:, lo:hi_occ].dot(occ_no_u)      # whole AVAS-pool occ block, NO basis
    no_vir_ao_full = mo0[:, hi_occ:hi_vir].dot(vir_no_u)  # whole AVAS-pool vir block, NO basis

    unsel_occ = [i for i in range(no_occ_ao_full.shape[1]) if i not in sel_occ]
    unsel_vir = [i for i in range(no_vir_ao_full.shape[1]) if i not in sel_vir]
    sel_occ_sorted = sorted(sel_occ)
    sel_vir_sorted = sorted(sel_vir)

    core_block = numpy.hstack((mo0[:, :lo], no_occ_ao_full[:, unsel_occ]))
    active_block = numpy.hstack((no_occ_ao_full[:, sel_occ_sorted], no_vir_ao_full[:, sel_vir_sorted]))
    virt_block = numpy.hstack((no_vir_ao_full[:, unsel_vir], mo0[:, hi_vir:]))

    mo_full = numpy.hstack((core_block, active_block, virt_block))
    ncore = core_block.shape[1]
    assert mo_full.shape[1] == mo0.shape[1], (mo_full.shape, mo0.shape)
    return mo_full, ncore


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
    ap.add_argument("--force-n-occ", type=int, default=None,
                    help="bypass occ_dev_cutoff/site_score_floor/total_active entirely and take "
                         "exactly this many occupied orbitals, ranked by occ_dev over the FULL "
                         "pool -- use this to realize a size cut_quality.select_space() already "
                         "found and verified against the --spectrum-out file, not to pick a new "
                         "size. Must be given together with --force-n-vir.")
    ap.add_argument("--force-n-vir", type=int, default=None,
                    help="counterpart to --force-n-occ for the virtual block.")
    ap.add_argument("--random-seed", type=int, default=None,
                    help="scramble control: with --force-n-occ/--force-n-vir, draw that many "
                         "orbitals UNIFORMLY AT RANDOM from the full pool instead of ranking by "
                         "occ_dev/vir_dev -- isolates selection pressure (site_score) from pure "
                         "size/filling.")
    ap.add_argument("--min-vir-frac", type=float, default=0.20)
    ap.add_argument("--max-vir-frac", type=float, default=0.55)
    ap.add_argument("--mp2-max-memory", type=int, default=128000,
                    help="MB passed as the MP2 object's max_memory (default 128000 -- pyscf's "
                         "own default ~4000 MB is sized for a laptop; the MP2 t2 amplitude "
                         "tensor for a few hundred active orbitals easily needs several GB, "
                         "and pyscf checks THIS ceiling, not the machine's real RAM, before "
                         "running -- raise it if you still see 'Insufficient memory for "
                         "holding t2 incore'")
    ap.add_argument("--chkfile", default=None,
                    help="pyscf SCF checkpoint path -- if it already exists, seeds (not skips) "
                         "SCF from it so a repeat run on the same cluster converges in ~1-2 "
                         "cycles instead of from scratch; if it doesn't exist yet, this run "
                         "creates it for a later localize_active_space.py/solange_dmrg.py call "
                         "on the SAME geometry/charge/spin/basis to reuse")
    ap.add_argument("--out", required=True)
    ap.add_argument("--spectrum-out", default=None,
                    help="also write {\"occupied\": [...], \"virtual\": [...]} with the FULL, "
                         "full-precision occ_dev/vir_dev spectra for the whole AVAS pool (not "
                         "just the selected subset, and not rounded for terminal printing) -- "
                         "feed this to check_active_space.py's --spectrum-json for cut_quality.py's "
                         "R7 check. Rounding to 4 decimals (what the terminal log shows) changes "
                         "q meaningfully when gaps are as small as 1e-4, so compute cut quality "
                         "from this file, never from printed output.")
    a = ap.parse_args()

    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
    from run_gate2_avas import build_mf

    mf = build_mf(a.xyz, a.charge, a.spin, basis=a.basis, chkfile=a.chkfile)
    result = two_criterion_select(mf, a.ao_set, site_threshold=a.site_threshold,
                                   occ_dev_cutoff=a.occ_dev_cutoff,
                                   site_score_floor=a.site_score_floor,
                                   total_active=a.total_active,
                                   min_vir_frac=a.min_vir_frac, max_vir_frac=a.max_vir_frac,
                                   mp2_max_memory=a.mp2_max_memory,
                                   force_n_occ=a.force_n_occ, force_n_vir=a.force_n_vir,
                                   random_seed=a.random_seed)
    n_occ, n_vir = len(result['sel_occ']), len(result['sel_vir'])
    ncas, nelec = n_occ + n_vir, 2 * n_occ
    print(f"[avas-mp2] FINAL selection: CAS({nelec},{ncas})  n_occ={n_occ} n_vir={n_vir} "
          f"({100*n_vir/ncas:.1f}% virtual)")

    if a.spectrum_out:
        # Full pool, full precision, sorted most-correlated-first -- exactly the
        # shape cut_quality.py/check_active_space.py's R7 expect. occ_dev is
        # ALREADY the hole occupation (2 - n), not the occupation itself (see
        # two_criterion_select's own comment) -- do not re-derive it here.
        spectrum = dict(occupied=sorted(result['occ_dev'], reverse=True),
                        virtual=sorted(result['vir_dev'], reverse=True))
        with open(a.spectrum_out, "w") as fh:
            json.dump(spectrum, fh)
        print(f"[avas-mp2] wrote {a.spectrum_out} (full pool: {len(spectrum['occupied'])} occ, "
              f"{len(spectrum['virtual'])} vir, full precision)")

    mo, ncore = build_mo(mf, result)
    nelectron = mf.mol.nelectron
    expected_ncore = (nelectron - nelec) // 2
    if ncore != expected_ncore:
        sys.exit(f"[avas-mp2] REFUSING: build_mo() returned ncore={ncore} but "
                 f"(nelectron - nelecas)//2 = {expected_ncore} (nelectron={nelectron}, "
                 f"nelecas={nelec}). ncore is arithmetic, not a choice -- this mismatch means "
                 f"mo_full is not full-width, or some AVAS-pool orbital was counted twice or "
                 f"dropped. Fix build_mo() before trusting this file's active_start_col.")
    print(f"[avas-mp2] ncore={ncore} (== (nelectron-nelecas)//2={expected_ncore}, verified) "
          f"active columns [{ncore}:{ncore+ncas}]")
    out = dict(xyz=a.xyz, charge=a.charge, spin=a.spin, basis=a.basis, ao_set=a.ao_set,
               site_threshold=a.site_threshold, occ_dev_cutoff=a.occ_dev_cutoff,
               site_score_floor=a.site_score_floor, chkfile=a.chkfile,
               # column where the active block starts in the saved _mo.npy -- this IS
               # pyscf's ncore (verified == (nelectron-nelecas)//2 above), needed by
               # localize_active_space.py to split-localize exactly the active occ/vir
               # blocks, and by solange_dmrg.py's --active-start-col to tell CASCI/CASSCF
               # where the active columns sit (its default ncore formula is only valid
               # for a full-width, standard-order mo_coeff).
               active_start_col=int(ncore),
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
