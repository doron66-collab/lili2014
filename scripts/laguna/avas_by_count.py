#!/usr/bin/env python3
"""
avas_by_count.py -- AVAS active-space selection by explicit orbital COUNT,
not by projection-score threshold.

Why this exists (2026-09-30, TP53_C275F cluster76 finding):

AVAS (pyscf.mcscf.avas) projects target atomic orbitals onto the occupied and
virtual MO spaces SEPARATELY and diagonalizes each projector, giving two
independent sigma-eigenvalue spectra (sigma_occ, sigma_vir in [0,1]). In a
closed-shell, saturated organic/protein cluster, the target-AO character lives
overwhelmingly in the occupied space (filled bonds), so sigma_occ clusters near
1 and sigma_vir clusters near 0. A single --threshold applied to BOTH spectra
is therefore a broken instrument: it either takes nearly all of the occupied
block and none of the virtual block, or takes everything. Verified on the
TP53_C275F 12-residue cluster: the occupied selection was IDENTICAL (111
orbitals) across thresholds 0.1/0.2/0.4/0.6 -- a factor-of-six change in
threshold produced zero change on the occupied side -- while the virtual
selection collapsed 70 -> 7 -> 0 -> 0 over that same range. Every threshold in
[0.4, 0.95] produced CAS(*, occ) with ZERO virtual orbitals: exactly one
Slater determinant, S_max = 0.0 guaranteed by construction, not measured.
See check_active_space.py for the general precondition this motivated.

Checked empirically before writing this (Laguna PySCF 2.14, project pin 2.13):
avas.avas()/avas.AVAS has no nocc_act/nvir_act parameter in this version --
there is no one-line API fix. This module is the "version-independent
fallback": run AVAS's own projector, but keep occ_weights/vir_weights (the
per-orbital sigma eigenvalues AVAS already computes and normally discards
after thresholding) and slice the top-N of EACH block explicitly, instead of
applying one threshold to both. Mirrors pyscf.mcscf.avas._kernel's own
algorithm (openshell_option=2 path) exactly, replacing its boolean threshold
masks with explicit count-based masks -- so the selected orbitals and their
canonicalization are identical to what AVAS itself would produce, just cut by
COUNT instead of VALUE.

Usage (library):
    from avas_by_count import avas_by_count
    ncas, nelecas, mo, occ_w, vir_w = avas_by_count(mf, aolabels, n_occ=32, n_vir=16)

Usage (CLI, reuses run_gate2_avas.py's cluster/SCF building):
    python avas_by_count.py --target TP53_C275F --xyz tp53_c275f_cluster.xyz \
        --charge -1 --spin 0 --basis 6-31g --n-occ 32 --n-vir 16 \
        --ao-set "S 3p, N 2p, O 2p" --out tp53_c275f_wt_cas64_48.json
"""
import argparse
import json
import sys

import numpy


def avas_by_count(mf, aolabels, n_occ, n_vir, minao='minao', with_iao=False,
                   openshell_option=2, canonicalize=True, ncore=0):
    """Same projector/eigendecomposition as pyscf.mcscf.avas._kernel, but
    selects the top n_occ occupied and top n_vir virtual orbitals BY COUNT
    (highest sigma first) instead of by a shared threshold. Returns
    (ncas, nelecas, mo_coeff, occ_weights_selected, vir_weights_selected).

    occ_weights_selected / vir_weights_selected are the sigma eigenvalues of
    the orbitals actually kept, sorted ascending as pyscf's eigh returns them
    -- report min/max of each to show how far down the ranked list this cut
    reached, since (unlike a threshold cut) a count cut has no natural score
    floor of its own.
    """
    from functools import reduce
    import scipy.linalg
    from pyscf import gto

    if openshell_option != 2:
        raise NotImplementedError("only openshell_option=2 (closed-shell / "
                                   "project singly-occupied as alpha) is implemented here")

    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    mo_energy = mf.mo_energy
    mol = mf.mol

    nocc = numpy.count_nonzero(mo_occ != 0)
    ovlp = mol.intor_symmetric('int1e_ovlp')

    pmol = mol.copy()
    pmol.atom = mol._atom
    pmol.unit = 'B'
    pmol.symmetry = False
    pmol.basis = minao
    pmol.build(False, False)

    aos = [s.strip() for s in aolabels.split(",")] if isinstance(aolabels, str) else aolabels
    baslst = pmol.search_ao_label(aos)

    if with_iao:
        from pyscf.lo import iao
        c = iao.iao(mol, mo_coeff[:, ncore:nocc], minao)[:, baslst]
        s2 = reduce(numpy.dot, (c.T, ovlp, c))
        s21 = reduce(numpy.dot, (c.T, ovlp, mo_coeff[:, ncore:]))
    else:
        s2 = pmol.intor_symmetric('int1e_ovlp')[baslst][:, baslst]
        s21 = gto.intor_cross('int1e_ovlp', pmol, mol)[baslst]
        s21 = numpy.dot(s21, mo_coeff[:, ncore:])
    sa = s21.T.dot(scipy.linalg.solve(s2, s21, assume_a='pos'))

    n_occ_pool = nocc - ncore
    n_vir_pool = mo_coeff.shape[1] - nocc
    if n_occ > n_occ_pool:
        raise ValueError(f"asked for n_occ={n_occ} but only {n_occ_pool} occupied "
                          f"orbitals exist above core")
    if n_vir > n_vir_pool:
        raise ValueError(f"asked for n_vir={n_vir} but only {n_vir_pool} virtual "
                          f"orbitals exist")

    # eigh returns ascending eigenvalues -- the LAST n_occ / n_vir columns are
    # the highest-sigma (most target-AO-like) orbitals in each block.
    wocc, u = numpy.linalg.eigh(sa[:n_occ_pool, :n_occ_pool])
    occ_sel = numpy.zeros(n_occ_pool, dtype=bool)
    occ_sel[-n_occ:] = True if n_occ > 0 else occ_sel
    nelecas = 2 * n_occ
    mocore = mo_coeff[:, ncore:nocc].dot(u[:, ~occ_sel])
    mocas_occ = mo_coeff[:, ncore:nocc].dot(u[:, occ_sel])

    wvir, uv = numpy.linalg.eigh(sa[n_occ_pool:, n_occ_pool:])
    vir_sel = numpy.zeros(n_vir_pool, dtype=bool)
    vir_sel[-n_vir:] = True if n_vir > 0 else vir_sel
    mocas_vir = mo_coeff[:, nocc:].dot(uv[:, vir_sel])
    movir = mo_coeff[:, nocc:].dot(uv[:, ~vir_sel])

    mocas = numpy.hstack((mocas_occ, mocas_vir))
    ncas = mocas.shape[1]
    occ_weights = wocc[occ_sel]
    vir_weights = wvir[vir_sel]

    mofreeze = mo_coeff[:, :ncore]
    if canonicalize:
        from pyscf.mcscf import dmet_cas

        def trans(c):
            if c.shape[1] == 0:
                return c
            csc = reduce(numpy.dot, (c.T, ovlp, mo_coeff))
            fock = numpy.dot(csc * mo_energy, csc.T)
            e, uu = scipy.linalg.eigh(fock)
            return dmet_cas.symmetrize(mol, e, numpy.dot(c, uu), ovlp, None)
        if ncore > 0:
            mofreeze = trans(mofreeze)
        mocore = trans(mocore)
        mocas = trans(mocas)
        movir = trans(movir)

    mo = numpy.hstack((mofreeze, mocore, mocas, movir))
    return ncas, nelecas, mo, occ_weights, vir_weights


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True, help="cluster geometry, e.g. from run_gate2_avas.py --cluster-only")
    ap.add_argument("--charge", type=int, required=True)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ao-set", required=True, help="AVAS AO criterion string, e.g. 'S 3p, N 2p, O 2p'")
    ap.add_argument("--n-occ", type=int, required=True)
    ap.add_argument("--n-vir", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
    from run_gate2_avas import build_mf  # reuse the project's SCF-building step

    mf = build_mf(a.xyz, a.charge, a.spin, basis=a.basis)
    ncas, nelecas, mo, occ_w, vir_w = avas_by_count(
        mf, a.ao_set, n_occ=a.n_occ, n_vir=a.n_vir)

    print(f"[avas-by-count] CAS({nelecas},{ncas})  n_occ={a.n_occ} n_vir={a.n_vir}  "
          f"-> {2*ncas} qubits under Jordan-Wigner")
    print(f"[avas-by-count] occ sigma range kept: [{occ_w.min():.4f}, {occ_w.max():.4f}]")
    print(f"[avas-by-count] vir sigma range kept: [{vir_w.min():.4f}, {vir_w.max():.4f}]  "
          f"(compare to the discarded pool's max, printed above by pyscf's own log if --verbose)")

    out = dict(xyz=a.xyz, charge=a.charge, spin=a.spin, basis=a.basis, ao_set=a.ao_set,
               n_occ=a.n_occ, n_vir=a.n_vir, ncas=int(ncas), nelecas=int(nelecas),
               qubits=int(2 * ncas), e_scf=float(mf.e_tot), scf_converged=bool(mf.converged),
               occ_sigma_min=float(occ_w.min()), occ_sigma_max=float(occ_w.max()),
               vir_sigma_min=float(vir_w.min()), vir_sigma_max=float(vir_w.max()))
    numpy.save(a.out.rsplit(".", 1)[0] + "_mo.npy", mo)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[avas-by-count] wrote {a.out} and {a.out.rsplit('.', 1)[0]}_mo.npy "
          f"(pass the .npy to solange_dmrg.py --load-orbitals)")


if __name__ == "__main__":
    main()
