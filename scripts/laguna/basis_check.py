#!/usr/bin/env python3
"""
basis_check.py -- did the orbital basis you loaded actually reach the Hamiltonian?

The question this answers came from a real null result: two DMRG runs on the same
CAS(36,34), one in canonical and one in split-Pipek-Mezey orbitals, returned
energies identical to 8 decimals at every bond dimension and S_max identical to
13 significant figures -- while the two orbital .npy files provably differed
(max|difference| = 1.24 over the active block).

Verifying that the INPUT FILES differ is not the same as verifying that the
difference survived into h1e/g2e. Anything between loading the coefficients and
building the integrals can silently undo a rotation confined to the active
space: an orbital-optimisation step (the CASSCF energy has no gradient along a
rotation inside the active space, so the optimiser is free to land anywhere in
that flat direction), an active-space canonicalisation, a natural-orbital
re-sort. This module tests the integrals themselves.

WHICH SCALARS CAN TEST IT -- verified numerically on a random orthogonal
rotation of a symmetrised (pq|rs) tensor, n = 6:

    SENSITIVE, usable                      INVARIANT, useless
    sum_p |h1e[p,p]|                       trace(h1e)
    sum_pq |h1e[p,q]|                      ||h1e||_F
    sum_p (pp|pp)   <- locality measure    ||g2e||_F
                                           sum_pq (pp|qq)
                                           sum_pq (pq|pq)

The trap is the right-hand column. A Frobenius norm is the natural thing to
reach for when fingerprinting a tensor, and ||g2e||_F is EXACTLY invariant under
a unitary rotation of all four indices -- so a norm-based fingerprint returns
the same number for two genuinely different bases and appears to confirm that
nothing changed. Use the left column.

sum_p (pp|pp) is the one with physical meaning: the on-site Coulomb repulsion
summed over orbitals is large when electrons are confined to single centres and
small when they are spread out. It is therefore both a validity check and a
direct measure of how localised the basis actually is.
"""
import json


def fingerprint(h1e, g2e, label=""):
    """Rotation-sensitive scalars of the integrals as handed to the solver.

    h1e : (n,n) array-like.  g2e : (n,n,n,n) array-like in chemist (pq|rs) order.
    Call this on the arrays passed INTO the DMRG driver, not on anything upstream.
    """
    import numpy as np
    h = np.asarray(h1e, dtype=float)
    g = np.asarray(g2e, dtype=float)
    n = h.shape[0]
    assert h.shape == (n, n) and g.shape == (n, n, n, n), (h.shape, g.shape)
    onsite = float(np.einsum("pppp->", g))
    return {
        "label": label, "n_orb": n,
        # --- sensitive to rotation: these are the test ---
        "sum_abs_h_diag": float(np.abs(np.diag(h)).sum()),
        "sum_abs_h_all": float(np.abs(h).sum()),
        "sum_onsite_coulomb": onsite,
        "mean_onsite_coulomb": onsite / n,
        # --- invariant: recorded only so a reader can see they are invariant ---
        "INVARIANT_trace_h": float(np.trace(h)),
        "INVARIANT_norm_h": float(np.linalg.norm(h)),
        "INVARIANT_norm_g": float(np.linalg.norm(g)),
    }


def compare(fp_a, fp_b, rtol=1e-10):
    """Did the two runs use different integrals? Decides on the sensitive scalars."""
    sens = ["sum_abs_h_diag", "sum_abs_h_all", "sum_onsite_coulomb"]
    inv = [k for k in fp_a if k.startswith("INVARIANT_")]
    rows, differ = [], False
    for k in sens + inv:
        a, b = fp_a[k], fp_b[k]
        scale = max(1.0, abs(a))
        same = abs(a - b) <= rtol * scale
        if k in sens and not same:
            differ = True
        rows.append({"quantity": k, "a": a, "b": b,
                     "relative_difference": abs(a - b) / scale, "identical": same})
    loc = None
    if fp_a["sum_onsite_coulomb"] > 0 and fp_b["sum_onsite_coulomb"] > 0:
        loc = fp_b["sum_onsite_coulomb"] / fp_a["sum_onsite_coulomb"]
    return {
        "integrals_differ": differ,
        "onsite_coulomb_ratio_b_over_a": loc,
        "rows": rows,
        "verdict": (
            "the two runs used DIFFERENT integrals -- the rotation reached the "
            "Hamiltonian, so an identical energy ladder is not an input problem "
            "and must be explained physically"
            if differ else
            "the two runs used the SAME integrals to within %.0e. The rotation did "
            "NOT reach the Hamiltonian: whatever was loaded, the solver was handed "
            "one basis twice. An identical energy ladder is then expected and says "
            "nothing about localisation. Look for an orbital-optimisation or "
            "canonicalisation step between loading the coefficients and building "
            "the integrals." % rtol),
    }


def _selftest():
    import numpy as np
    rng = np.random.default_rng(7)
    n = 6
    A = rng.normal(size=(n, n)); h = 0.5 * (A + A.T)
    B = rng.normal(size=(n, n, n, n))
    g = B + B.transpose(1, 0, 2, 3) + B.transpose(0, 1, 3, 2) + B.transpose(1, 0, 3, 2)
    g = g + g.transpose(2, 3, 0, 1)
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    hp = Q.T @ h @ Q
    gp = np.einsum("ap,bq,cr,ds,abcd->pqrs", Q, Q, Q, Q, g, optimize=True)

    a, b = fingerprint(h, g, "unrotated"), fingerprint(hp, gp, "rotated")
    r = compare(a, b)
    assert r["integrals_differ"], "a real rotation must be detected"
    for row in r["rows"]:
        if row["quantity"].startswith("INVARIANT_"):
            assert row["identical"], ("this quantity is meant to be invariant", row)
    print("selftest rotation detected   : integrals_differ =", r["integrals_differ"])
    print("selftest invariants held     : trace/norms unchanged under rotation")
    same = compare(a, fingerprint(h, g, "same again"))
    assert not same["integrals_differ"]
    print("selftest identical input     : integrals_differ =", same["integrals_differ"])
    print("selftest on-site ratio       : %.4f (localisation indicator)"
          % r["onsite_coulomb_ratio_b_over_a"])


if __name__ == "__main__":
    _selftest()
    print(json.dumps({"usage": "fingerprint(h1e, g2e) in each run, then compare(fp_a, fp_b)"}))
