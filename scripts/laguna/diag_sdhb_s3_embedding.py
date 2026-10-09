"""Diagnose the SDHB S3 high-spin embedding failure (job 2350799).

The S=15/2 DMRG leg converged smoothly (M=250/500/1000) to an energy
3.6 Ha ABOVE E_ROHF. That cannot happen if the ROHF determinant lies inside
the active space, so one of these is true:

  (a) the AVAS core/active split does not contain the ROHF determinant
      (default openshell_option=2 projects singly- and doubly-occupied
      orbitals together, so "core" orbitals can carry SOMO character), or
  (b) h1e/h2e/ecore are built wrong.

This script separates the two without any DMRG, in minutes:
  1. core purity  -- ROHF alpha/beta occupation of every AVAS core orbital
                     (must all be 1.000 if the determinant is in the space);
  2. active count -- ROHF electrons projected onto the active orbitals
                     (must equal nelecas);
  3. E_det        -- energy of the ROHF determinant evaluated with the SAME
                     h1e/h2e/ecore the DMRG legs used (must equal E_ROHF).
Then repeats 1-3 with AVAS openshell_option=3 (SOMOs kept fully active).

Reads only files the run already wrote: the SCF chkfile and the saved AVAS
orbitals. Nothing is overwritten.
"""
import argparse
import json
import os

import numpy as np
from pyscf import ao2mo, lib, mcscf, scf
from pyscf.mcscf import avas


def check(mf, mo, ncas, nelecas, label):
    mol = mf.mol
    S = mf.get_ovlp()
    ncore = (mol.nelectron - nelecas) // 2
    occ = mf.mo_occ
    Ca = mf.mo_coeff[:, occ > 0]
    Cb = mf.mo_coeff[:, occ > 1]
    Da, Db = Ca @ Ca.T, Cb @ Cb.T

    def proj(C):
        return C.T @ S @ Da @ S @ C, C.T @ S @ Db @ S @ C

    core, act, virt = mo[:, :ncore], mo[:, ncore:ncore + ncas], mo[:, ncore + ncas:]
    ca, cb = proj(core)
    aa, ab = proj(act)
    va, vb = proj(virt)
    print(f"\n=== {label}: CAS({nelecas},{ncas}), ncore={ncore} ===")
    print(f"[core]   alpha occ min={np.diag(ca).min():.6f}  beta occ min={np.diag(cb).min():.6f}"
          f"  (both must be 1.000000)")
    print(f"[active] ROHF electrons in active space: alpha={np.trace(aa):.4f} "
          f"beta={np.trace(ab):.4f} total={np.trace(aa) + np.trace(ab):.4f} (must be {nelecas})")
    print(f"[virt]   ROHF electrons leaking into virtuals: {np.trace(va) + np.trace(vb):.2e}")

    mc = mcscf.CASCI(mf, ncas, nelecas)
    h1e, ecore = mc.get_h1eff(mo_coeff=mo)
    h2e = ao2mo.restore(1, mc.get_h2eff(mo), ncas)
    J = np.einsum("pqrs,rs->pq", h2e, aa + ab)
    Ka = np.einsum("psrq,rs->pq", h2e, aa)
    Kb = np.einsum("psrq,rs->pq", h2e, ab)
    e_det = (ecore + np.einsum("pq,pq", h1e, aa + ab)
             + 0.5 * (np.einsum("pq,pq", J, aa + ab)
                      - np.einsum("pq,pq", Ka, aa) - np.einsum("pq,pq", Kb, ab)))
    print(f"[energy] ecore={ecore:.8f}  E_det(ROHF in CAS integrals)={e_det:.8f}  "
          f"E_ROHF={mf.e_tot:.8f}  diff={(e_det - mf.e_tot) * 1000:+.3f} mHa (must be ~0)")
    return e_det


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chkfile", default="sdhb_s3_scf.chk")
    ap.add_argument("--mo", default="sdhb_s3_mo_coeff_highspin.npy")
    ap.add_argument("--json", default="sdhb_s3_both_spins.json")
    ap.add_argument("--threshold", type=float, default=0.2)
    a = ap.parse_args()

    mol = lib.chkfile.load_mol(a.chkfile)
    mol.verbose = 0
    data = scf.chkfile.load(a.chkfile, "scf")
    mf = scf.ROHF(mol)
    mf.mo_coeff, mf.mo_occ, mf.mo_energy = data["mo_coeff"], data["mo_occ"], data["mo_energy"]
    mf.e_tot = mf.energy_tot(mf.make_rdm1())
    mf.converged = True
    print(f"[scf] nelectron={mol.nelectron} spin={mol.spin} nao={mol.nao} "
          f"E_ROHF(recomputed from chkfile)={mf.e_tot:.8f}")

    if os.path.exists(a.json):
        j = json.load(open(a.json))
        ncas, nelecas = j["ncas"], j["nelecas"]
    else:
        ncas = nelecas = None

    if os.path.exists(a.mo) and ncas is not None:
        check(mf, np.load(a.mo), ncas, nelecas, "orbitals the DMRG legs actually used")
    else:
        print(f"[skip] {a.mo} or {a.json} not found -- checking a fresh AVAS only")

    for opt in (2, 3):
        n, ne, mo = avas.avas(mf, ["Fe 3d", "S 3p"], threshold=a.threshold,
                              openshell_option=opt, verbose=0)
        check(mf, mo, int(n), int(ne), f"fresh AVAS openshell_option={opt}")


if __name__ == "__main__":
    main()
