#!/usr/bin/env python3
"""
check_rotation_no_solver.py -- did split-PM localization reach the Hamiltonian?
No kernel(), no fcisolver, no DMRG -- just the AO->active-MO integral transform
called directly on each orbital array. Seconds, not minutes.

Claude Science (2026-10-02), having read pyscf/mcscf/casci.py directly: the
earlier fingerprint script's maxM=250 solver runs were unnecessary and two of
their four cells meaningless -- get_h1eff(mo_coeff)/get_h2eff(mo_coeff), called
with the orbital array passed EXPLICITLY, build the integrals straight from
that array regardless of canonicalization or any solver, so this answers the
same question in seconds: does the physically-different orbital file actually
produce physically-different integrals.
"""
import sys
import numpy as np
from pyscf import mcscf, ao2mo

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from run_gate2_avas import build_mf
import basis_check as bc

XYZ = "tp53_c275f_wt_cluster.xyz"
CHKFILE = "tp53_c275f_wt.chk"
MO_CANONICAL = "tp53_c275f_wt_cas34_measured_mo.npy"
MO_LOCALIZED = "tp53_c275f_wt_cas34_localized_mo.npy"
NCAS, NELEC = 34, 36

mf = build_mf(XYZ, 0, 0, basis="6-31g", chkfile=CHKFILE)

mo_can = np.load(MO_CANONICAL)
mo_loc = np.load(MO_LOCALIZED)
print("mo diff (full matrix):", np.max(np.abs(mo_can - mo_loc)))

mc = mcscf.CASCI(mf, NCAS, NELEC)

h1_can, ecore_can = mc.get_h1eff(mo_can)
h2_can = ao2mo.restore(1, mc.get_h2eff(mo_can), NCAS)
h1_loc, ecore_loc = mc.get_h1eff(mo_loc)
h2_loc = ao2mo.restore(1, mc.get_h2eff(mo_loc), NCAS)

print("ecore: canonical=%.8f  localized=%.8f  diff=%.2e" % (ecore_can, ecore_loc, abs(ecore_can - ecore_loc)))

fp_can = bc.fingerprint(h1_can, h2_can, "canonical")
fp_loc = bc.fingerprint(h1_loc, h2_loc, "localized")
r = bc.compare(fp_can, fp_loc)
print()
print("sum_abs_h_diag:    canonical=%.6f  localized=%.6f" % (fp_can["sum_abs_h_diag"], fp_loc["sum_abs_h_diag"]))
print("sum_abs_h_all:     canonical=%.6f  localized=%.6f" % (fp_can["sum_abs_h_all"], fp_loc["sum_abs_h_all"]))
print("sum_onsite_coulomb: canonical=%.6f  localized=%.6f  ratio=%.4f"
      % (fp_can["sum_onsite_coulomb"], fp_loc["sum_onsite_coulomb"], r["onsite_coulomb_ratio_b_over_a"]))
print("INVARIANT_norm_g:  canonical=%.6f  localized=%.6f (should match -- sanity check on the fingerprint itself)"
      % (fp_can["INVARIANT_norm_g"], fp_loc["INVARIANT_norm_g"]))
print()
print("VERDICT:", r["verdict"])
