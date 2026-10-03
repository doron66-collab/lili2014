#!/usr/bin/env python3
"""
fiedler_mean_range.py -- how much long-range coupling survives the best 1D
orbital ordering DMRG can find, for a given active space.

Why this matters (Claude Science, 2026-10-03): a scramble control (random
orbitals drawn from a real AVAS pool, instead of the site-selected subset)
isolates selection pressure from active-space size -- but a random draw from a
spatially delocalized pool gives orbitals scattered in space, and a scattered
set cannot be laid out on a 1D DMRG chain without long-range couplings
surviving. Long-range coupling inflates the bipartite entanglement entropy at
every cut and slows DMRG convergence for a purely GEOMETRIC reason, unrelated
to the real physical correlation this project's S_max is meant to measure. A
scramble control that "lights up" is ambiguous without this diagnostic: it
could mean "this active-space size alone drives Class A" (the real finding) or
just "scattered orbitals are hard for DMRG" (an artifact of the control's own
construction, not a property of the molecule).

fiedler_order(K): the same spectral reordering DMRG codes use to choose a
1D orbital sequence -- order sites by the Fiedler vector (the eigenvector of
the graph Laplacian built from |K| for the second-smallest eigenvalue) of the
orbital coupling/exchange matrix K.

mean_range(K, perm): the coupling-weighted mean |i-j| distance after applying
permutation perm -- small for a compact, well-orderable set; large for a
scattered one no permutation can fix.
"""
import sys
import numpy


def fiedler_order(K):
    """Order sites by the Fiedler vector of the graph Laplacian built from |K|.
    K: real symmetric coupling/exchange matrix, shape (n, n).
    Returns a permutation (array of indices) ordering the sites along the
    algebraic-connectivity axis -- the same spectral heuristic DMRG codes use
    to choose a 1D orbital sequence from a general coupling matrix.
    """
    W = numpy.abs(numpy.asarray(K, dtype=float))
    numpy.fill_diagonal(W, 0.0)
    D = numpy.diag(W.sum(axis=1))
    L = D - W
    # Symmetric Laplacian is PSD; smallest eigenvalue is ~0 (constant vector),
    # the Fiedler vector is the eigenvector of the SECOND smallest eigenvalue.
    vals, vecs = numpy.linalg.eigh(L)
    fiedler_vec = vecs[:, 1]
    return numpy.argsort(fiedler_vec)


def mean_range(K, perm):
    """Coupling-weighted mean |i-j| distance under permutation perm.
    K: coupling/exchange matrix. perm: a permutation of range(len(K)).
    """
    K = numpy.asarray(K, dtype=float)
    Kp = K[numpy.ix_(perm, perm)]
    W = numpy.abs(Kp)
    numpy.fill_diagonal(W, 0.0)
    n = len(W)
    d = numpy.abs(numpy.subtract.outer(numpy.arange(n), numpy.arange(n)))
    return float((W * d).sum() / W.sum())


def exchange_matrix_from_h2e(h2e, ncas):
    """K[i,j] = (ij|ji) in chemist (pq|rs) notation -- h2e indexed [p,q,r,s]."""
    h2e = numpy.asarray(h2e).reshape(ncas, ncas, ncas, ncas)
    K = numpy.zeros((ncas, ncas))
    for i in range(ncas):
        for j in range(ncas):
            K[i, j] = h2e[i, j, j, i]
    return K


def _selftest():
    rng = numpy.random.default_rng(0)
    n = 34

    # "Compact" toy: exponentially decaying coupling with distance in a
    # PLANTED good order (i.e. a set that IS naturally 1D-orderable) --
    # models the site-selected case, where nearby orbitals (in the real,
    # natural sense) are chemically similar and strongly coupled, distant
    # ones weakly so.
    idx = numpy.arange(n)
    dist = numpy.abs(numpy.subtract.outer(idx, idx))
    K_compact = numpy.exp(-dist / 3.0) + 0.01 * rng.standard_normal((n, n))
    K_compact = (K_compact + K_compact.T) / 2

    # "Scattered" toy: coupling assigned with NO relation to any 1D order --
    # models a random draw from a spatially delocalized pool, where orbital
    # index carries no spatial meaning at all.
    K_scattered = rng.random((n, n))
    K_scattered = (K_scattered + K_scattered.T) / 2

    perm_c = fiedler_order(K_compact)
    perm_s = fiedler_order(K_scattered)
    mr_c = mean_range(K_compact, perm_c)
    mr_s = mean_range(K_scattered, perm_s)
    print("selftest compact (planted 1D structure): mean_range = %.4f" % mr_c)
    print("selftest scattered (no structure)       : mean_range = %.4f" % mr_s)
    assert mr_s > mr_c, (mr_c, mr_s)
    print("selftest OK: scattered > compact by %.2fx" % (mr_s / mr_c))

    # Also verify natural (identity) order is clearly worse than Fiedler order
    # on the compact case -- if the Fiedler order recovers something close to
    # the planted identity order (up to reflection), mean_range at Fiedler
    # order should be close to mean_range at identity order, both much better
    # than a random permutation.
    perm_rand = rng.permutation(n)
    mr_rand = mean_range(K_compact, perm_rand)
    print("compact at a RANDOM permutation          : mean_range = %.4f "
          "(sanity: should be >> Fiedler-ordered %.4f)" % (mr_rand, mr_c))
    assert mr_rand > mr_c


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        _selftest()
    else:
        print(__doc__)
        print("run with --selftest to verify fiedler_order/mean_range")
