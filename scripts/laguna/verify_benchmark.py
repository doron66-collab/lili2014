#!/usr/bin/env python3
"""
verify_benchmark.py — recompute every reference value in benchmarks.json from its
own stated specification, and report whether it reproduces.

WHY THIS EXISTS
A reference value nobody can regenerate is an assertion, not a reference. This
script closes that gap: it reads only the specification recorded in the registry
(geometry, basis, charge, spin, active space, convergence tolerances), recomputes
from scratch, and compares. If a stored value cannot be reproduced from its own
stated conditions, the specification is incomplete and the entry is not usable as
a benchmark — which is itself the finding, and worth failing loudly for.

This is the same posture LEON takes toward records, applied to reference data:
verify, don't trust. The registry does not get believed because we wrote it.

USAGE
    python scripts/laguna/verify_benchmark.py                # verify all entries
    python scripts/laguna/verify_benchmark.py --id N2_CAS66_ccpVDZ_re
    python scripts/laguna/verify_benchmark.py --tol-mha 0.01 # tighter agreement

Exit code is non-zero if any entry fails to reproduce, so this can gate CI.
"""
import argparse
import json
import sys
from pathlib import Path

REGISTRY = Path(__file__).with_name("benchmarks.json")

# Reproduction tolerance. Far tighter than chemical accuracy on purpose: this is
# not asking "is the chemistry good", it is asking "does the same specification
# yield the same number". Anything above numerical noise means the spec is
# under-determined or the stored value is wrong.
DEFAULT_REPRO_TOL_MHA = 0.01


def recompute(entry, verbose=0):
    """Recompute a benchmark's reference energy from its stated specification alone."""
    from pyscf import gto, scf, mcscf

    mol = gto.Mole()
    mol.atom    = entry["geometry"].replace("/", "\n")
    mol.basis   = entry["basis_set"]
    mol.charge  = entry.get("charge", 0)
    mol.spin    = entry.get("spin", 0)
    mol.unit    = entry.get("geometry_units", "angstrom")
    mol.verbose = verbose
    mol.build()

    conv = entry.get("convergence", {})
    mf = scf.RHF(mol)
    mf.conv_tol  = conv.get("rhf_conv_tol", 1e-10)
    mf.max_cycle = 400
    e_rhf = mf.kernel()
    if not mf.converged:
        raise RuntimeError("RHF did not converge from the stated specification")

    mc = mcscf.CASSCF(mf, ncas=entry["ncas"], nelecas=entry["nelecas"])
    mc.conv_tol        = conv.get("casscf_conv_tol", 1e-9)
    mc.conv_tol_grad   = conv.get("casscf_conv_tol_grad", 1e-5)
    mc.max_cycle_macro = 200
    mc.verbose         = verbose
    e_casscf = mc.kernel()[0]
    if not mc.converged:
        raise RuntimeError("CASSCF did not converge from the stated specification")

    return {"e_rhf": e_rhf, "e_casscf": e_casscf,
            "n_basis": mol.nao_nr(), "converged": True}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--id", default=None, help="verify only this benchmark id")
    ap.add_argument("--tol-mha", type=float, default=DEFAULT_REPRO_TOL_MHA,
                    help=f"reproduction tolerance in mHa (default {DEFAULT_REPRO_TOL_MHA})")
    ap.add_argument("--verbose", type=int, default=0)
    args = ap.parse_args()

    reg = json.loads(REGISTRY.read_text())
    entries = reg["benchmarks"]
    if args.id:
        entries = [e for e in entries if e["id"] == args.id]
        if not entries:
            print(f"No benchmark with id '{args.id}'. Available: "
                  f"{', '.join(e['id'] for e in reg['benchmarks'])}", file=sys.stderr)
            return 2

    print(f"Verifying {len(entries)} benchmark(s) against their own stated specifications")
    print(f"Reproduction tolerance: {args.tol_mha} mHa\n")

    failures = 0
    for e in entries:
        print(f"── {e['id']}  ({e['system']}, CAS({e['nelecas']},{e['ncas']})/{e['basis_set']})")
        try:
            got = recompute(e, verbose=args.verbose)
        except Exception as exc:
            print(f"   FAILED to recompute: {exc}\n")
            failures += 1
            continue

        stored = e["reference_energy_ha"]
        delta_mha = (got["e_casscf"] - stored) * 1000.0
        ok = abs(delta_mha) <= args.tol_mha

        print(f"   stored     : {stored:.10f} Ha")
        print(f"   recomputed : {got['e_casscf']:.10f} Ha")
        print(f"   delta      : {delta_mha:+.6f} mHa   -> {'REPRODUCES' if ok else 'MISMATCH'}")

        # Cross-check the supporting values too — they are part of the claim.
        sup = e.get("supporting_values", {})
        if "e_rhf_ha" in sup:
            d_rhf = (got["e_rhf"] - sup["e_rhf_ha"]) * 1000.0
            print(f"   RHF delta  : {d_rhf:+.6f} mHa")
            ok = ok and abs(d_rhf) <= args.tol_mha
        if e.get("n_basis_functions") and got["n_basis"] != e["n_basis_functions"]:
            print(f"   basis size MISMATCH: stated {e['n_basis_functions']}, got {got['n_basis']}")
            ok = False

        if not ok:
            failures += 1
        print()

    if failures:
        print(f"{failures} benchmark(s) did NOT reproduce. The stored value or the "
              f"specification is wrong — do not use these as references until resolved.")
        return 1
    print("All benchmarks reproduce from their stated specifications.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
