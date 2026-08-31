#!/usr/bin/env python3
"""
solange_shci.py — SHCI cross-validation for a SOLANGE DMRG classification.

DMRG's own Class A verdict is DMRG-specific evidence, not a cross-method
exclusion (solange_dmrg.py's own comment above classify()): a single method's
non-convergence proves that method struggles, not that classical computation
as a whole fails. This script closes that gap with a second, independent
solver — SHCI (Selected CI with perturbative correction, via the Dice
implementation), whose failure mechanism is determinant-sparsity growth,
unrelated to DMRG's bond-dimension growth across a one-dimensional chain.

It builds the IDENTICAL active-space Hamiltonian (same h1e/h2e/ecore) a DMRG
classification used, runs SHCI over it, and submits the result to SOLANGE
tagged against that DMRG record's id. The backend — not this script — computes
the agreement verdict, by fetching the referenced DMRG record's own stored
energy at ingestion time (DP1, "verify, don't trust": a submitting script's own
claim of agreement is never taken on faith).

Does NOT decide the active space, the charge, the basis, or the AVAS
criterion — it takes exactly the same cluster and criterion the DMRG run
used, so the only variable is the solver. Requires a Dice build (see
scripts/laguna/RUN_GUIDE.md for the from-source build recipe) and a
--dmrg-classification-id naming the DMRG record this run validates — a
cross-validation with nothing to validate against is not one.

USAGE
  python3 solange_shci.py --geometry cluster.xyz --charge 1 --basis sto-3g \\
      --avas "N 2p, O 2p" --threshold 0.2 --key DNMT3A_R882 \\
      --dmrg-classification-id <uuid-from-dmrg-submit-response> \\
      --dice-scripts ~/lili2014/Dice/scripts \\
      --sweep-eps 1e-3,5e-4,1e-4 \\
      --submit
"""
import argparse
import hashlib
import json
import sys
import time
import uuid
from pathlib import Path


# ── LEON seal (self-contained — mirrors backend/routes/leon.py's generic seal
# bit-for-bit and solange_dmrg.py's own _seal_payload, deliberately duplicated
# rather than imported: this script must run standalone on Laguna with no
# dependency on the SOLANGE backend package). LEON re-verifies this at
# ingestion; a mismatch is REJECTED, not stored.
def _seal_payload(record, exclude):
    return json.dumps({k: v for k, v in record.items() if k not in exclude},
                      sort_keys=True, default=str)


def _submit_shci(api, out):
    """POST the sealed SHCI cross-validation to SOLANGE. Best-effort: prints the
    outcome but never raises — a submit failure must not discard the local JSON
    already written to --out."""
    import urllib.request
    try:
        url = api.rstrip("/") + "/api/simulate/hpc/shci/submit"
        body = json.dumps(out, default=str).encode()
        req = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            resp = json.loads(r.read().decode())
        print(f"SUBMITTED → {url}")
        print(f"  status={resp.get('status')}  seal_ok={resp.get('seal_ok')}  "
              f"db={resp.get('db_status')}  run_id={resp.get('run_id')}")
        if "agreement" in resp:
            print(f"  delta_mha={resp.get('delta_mha')}  agreement={resp.get('agreement')}"
                  f"  (computed server-side against the referenced DMRG record)")
    except Exception as e:
        print(f"  SUBMIT FAILED (result is still safe locally in --out): {e}", file=sys.stderr)
        if hasattr(e, "read"):
            try:
                print(f"  server said: {e.read().decode()}", file=sys.stderr)
            except Exception:
                pass


def read_xyz(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    if lines and lines[0].strip().isdigit():
        lines = lines[2:]
    return "\n".join(lines)


def detect_hardware(n_threads):
    """Best-effort hostname/CPU tag — Dice/SHCI is CPU (+MPI) only, no GPU path."""
    import platform
    return f"CPU · {platform.node()} · {n_threads} threads"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geometry", required=True,
                    help="the SAME .xyz the DMRG classification used")
    ap.add_argument("--avas", required=True,
                    help="the SAME AVAS criterion the DMRG classification used")
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0,
                    help="2S for the active space; must match the DMRG run")
    ap.add_argument("--threshold", type=float, default=0.2)
    ap.add_argument("--max-memory", type=int, default=16000)
    ap.add_argument("--df-auxbasis", default="def2-universal-jkfit")
    ap.add_argument("--dice-scripts", required=True,
                    help="path to Dice/scripts (contains shci.py)")
    ap.add_argument("--sweep-eps", default="1e-3,5e-4,1e-4",
                    help="comma-separated epsilon schedule, largest (loosest) first")
    ap.add_argument("--num-thrds", type=int, default=8)
    ap.add_argument("--key", required=True, help="target key, e.g. DNMT3A_R882")
    ap.add_argument("--dmrg-classification-id", required=True,
                    help="id of the DMRG record (from its /hpc/dmrg/submit response) "
                         "this run cross-validates against")
    ap.add_argument("--out", default=".")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    default=None, metavar="API_URL",
                    help="POST the sealed result to SOLANGE (default URL if given bare)")
    args = ap.parse_args()

    sys.path.insert(0, str(Path(args.dice_scripts).expanduser()))
    try:
        from shci import SHCI
    except ImportError as e:
        sys.exit(f"could not import SHCI from {args.dice_scripts}: {e}\n"
                 f"  check --dice-scripts points at the directory containing shci.py")

    from pyscf import gto, scf, ao2mo, mcscf
    from pyscf.mcscf import avas

    aolabels = [s.strip() for s in args.avas.split(",")]
    mol = gto.M(atom=read_xyz(args.geometry), basis=args.basis, charge=args.charge,
                spin=args.spin, verbose=0, max_memory=args.max_memory)
    print(f"cluster: {mol.natm} atoms, {mol.nao} basis functions, "
          f"{mol.nelectron} electrons, charge={args.charge}, spin={args.spin}")

    mf = scf.RHF(mol).density_fit(auxbasis=args.df_auxbasis)
    mf.max_memory = args.max_memory
    mf.kernel()
    if not mf.converged:
        sys.exit("RHF did not converge.")
    print(f"RHF E = {mf.e_tot:.8f} Ha")

    ncas, nelec, mo = avas.avas(mf, aolabels, threshold=args.threshold)
    occ = nelec // 2
    print(f"AVAS(threshold={args.threshold}) -> CAS({nelec},{ncas})  "
          f"{occ} occupied, {ncas - occ} virtual")

    # Same effective-Hamiltonian construction the DMRG path uses (solange_hpc.py's
    # run_casscf / integrals_from_geometry) -- this is the point: identical
    # h1e/h2e/ecore handed to a different solver.
    mc = mcscf.CASCI(mf, ncas, nelec)
    h1e, ecore = mc.get_h1eff(mo_coeff=mo)
    h2e = ao2mo.restore(1, mc.get_h2eff(mo), ncas)

    eps_schedule = [float(x) for x in args.sweep_eps.split(",")]
    solver = SHCI(mol)
    solver.sweep_iter = list(range(len(eps_schedule)))
    solver.sweep_epsilon = eps_schedule
    solver.num_thrds = args.num_thrds
    print(f"SHCI epsilon schedule: {eps_schedule}")

    # shci.py's kernel expects (n_alpha, n_beta), not a total electron count.
    nelec_ab = ((nelec + args.spin) // 2, (nelec - args.spin) // 2)
    t0 = time.time()
    result = solver.kernel(h1e, h2e, ncas, nelec_ab, ecore=ecore)
    elapsed_s = round(time.time() - t0, 1)
    # pyscf FCI-solver convention returns (e, civec); some external solvers
    # return e alone -- handle both without guessing which silently.
    e_shci = result[0] if isinstance(result, (tuple, list)) else result
    print(f"SHCI E = {e_shci:.8f} Ha")
    print(f"elapsed {elapsed_s}s")

    out = {
        "id": str(uuid.uuid4()),
        "key": args.key,
        "dmrg_classification_id": args.dmrg_classification_id,
        "ncas": ncas, "nelec": nelec,
        "e_shci": e_shci,
        "sweep_eps": ",".join(str(e) for e in eps_schedule),
        "method": "SHCI (Dice, semistochastic heat-bath CI) cross-validation",
        "elapsed_s": elapsed_s,
        "provenance_source": "HPC/Laguna (SHCI cross-validation)",
        "hardware": detect_hardware(args.num_thrds),
    }
    # Seal at source (LEON re-verifies at ingestion — a mismatch is rejected, not
    # trusted). e_dmrg_ref/delta_mha/agreement are NOT computed here on purpose —
    # the backend computes them from the DMRG record's own current stored value
    # (DP1), so this script cannot assert an agreement the server hasn't checked.
    out["shci_seal_payload"] = _seal_payload(out, exclude={"shci_hash", "shci_seal_payload"})
    out["shci_hash"] = hashlib.sha256(out["shci_seal_payload"].encode()).hexdigest()

    Path(args.out).mkdir(parents=True, exist_ok=True)
    p = Path(args.out) / f"shci_crossval_{args.key}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"WROTE {p}")
    if args.submit:
        _submit_shci(args.submit, out)
    print("=" * 68)
    print("NOTE: SHCI is classical (CPU/MPI), NOT quantum — same as DMRG.")
    print("This run only cross-validates DMRG's own energy at the SAME active space;")
    print("it does not re-derive the active space, the charge, or the AVAS criterion.")


if __name__ == "__main__":
    main()
