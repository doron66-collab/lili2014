#!/usr/bin/env python3
"""
check_rotation_reached_hamiltonian.py -- did split-PM localization survive into
the h1e/h2e handed to the DMRG solver, on the REAL CAS(36,34) target?

No DMRG here at all -- just SCF (seeded from an existing --chkfile, so ~1-2
cycles not 30-45 min) + one CASCI call per orbital file + basis_check's
rotation-sensitive fingerprint. Minutes, decisive.

Tests canonicalization=False explicitly (pyscf's CASCI default canonicalizes
active orbitals post-solve, which is itself a rotation confined to
occ-among-occ/vir-among-vir -- exactly the flat direction that could erase a
prior localization). Reports fingerprints under BOTH settings so the comparison
is not assumed, only verified.
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


ACTIVE_START_COL = 106  # from tp53_c275f_wt_cas34_measured.json's active_start_col -- see
# check_rotation_no_solver.py's comment: pyscf's default ncore assumes the
# active space sits immediately below the Fermi level, which avas_mp2_select.py's
# output does NOT follow. Without this, get_h1eff/get_h2eff silently read the
# wrong 34 columns (confirmed: default ncore=329 here lands in the untouched
# virtual-beyond-AVAS-pool block, not the real active columns 106:140).


def get_fp(mo, label, canonicalization):
    from dmrgscf_block2 import Block2FCISolver
    mc = mcscf.CASCI(mf, NCAS, NELEC)
    mc.ncore = ACTIVE_START_COL
    mc.canonicalization = canonicalization
    # CAS(36,34) has FCI dim 4.86e18 -- pyscf's DEFAULT fcisolver is exact FCI,
    # which is exactly the impossible allocation that killed two earlier jobs.
    # Must use the same Block2FCISolver --dmrg-scf uses in the real probe, or
    # this call either crashes or never returns. maxM kept small (this is a
    # fingerprint check, not an energy measurement) but still has to be a real
    # DMRG solve, not the exact solver.
    # maxM=50 (tried first) failed Block2FCISolver's own RDM-vs-energy
    # self-check (Delta=1.06e-2 Ha) on this CAS(36,34) -- confirmed via
    # dmrgscf_block2.py --diagnose on an independent small system that
    # _BLOCK2_TO_PYSCF_2PDM_AXES=(0,3,1,2) is still correct, so that was
    # maxM=50 being too undertrained for reliable RDM extraction on a space
    # this large, not an axis bug. maxM=250 already passed this exact check
    # in the real probe run (validate() + the dmrg-scf solve itself both
    # succeeded there) -- reuse that value here.
    mc.fcisolver = Block2FCISolver(maxM=250, scratch=f"./tmp_fp_{label}",
                                   n_threads=16, stack_mem_gb=16)
    e = mc.kernel(mo)[0]
    mo_changed = np.max(np.abs(mc.mo_coeff - mo))
    print(f"  {label} (canonicalization={canonicalization}): E={e:.8f}  "
          f"|mo_coeff - input| = {mo_changed:.6e}")
    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), NCAS)
    return bc.fingerprint(h1e, h2e, label)


for canon in (True, False):
    print(f"\n=== canonicalization={canon} ===")
    fp_can = get_fp(mo_can, "canonical", canon)
    fp_loc = get_fp(mo_loc, "localized", canon)
    r = bc.compare(fp_can, fp_loc)
    print("sum_onsite_coulomb: canonical=%.6f  localized=%.6f  ratio=%.4f"
          % (fp_can["sum_onsite_coulomb"], fp_loc["sum_onsite_coulomb"],
             r["onsite_coulomb_ratio_b_over_a"]))
    print("VERDICT:", r["verdict"])
