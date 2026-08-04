#!/usr/bin/env python3
"""
resubmit_hpc.py — restore already-completed classical HPC/DMRG-SCF runs to
SOLANGE WITHOUT re-running them.

Same idea as resubmit_qpu.py (which does this for real quantum-hardware runs),
generalized: any solange_hpc.py run that was NOT started with --submit still
wrote a full, LEON-sealed P1-P9 record to its --out directory. That record is
complete and already sealed — the P8 hash is already computed, the CASSCF/
DMRG-SCF already ran — so getting it onto the platform is a matter of POSTing
the same JSON solange_hpc.py's own --submit would have sent, not recomputing
anything. This is exactly the situation a long DMRG-SCF run (past the FCI
wall, an hour or more) creates if you want to see the result before deciding
whether it is worth spending real time re-running with --submit attached.

  * Nothing is re-executed — this only re-POSTs what is already on disk.
  * Safe to run even if the row already exists: the backend upserts by id.

USAGE:
  python scripts/laguna/resubmit_hpc.py                 # scan ./out, re-submit HPC rows
  python scripts/laguna/resubmit_hpc.py --out ./out --dry-run   # list, don't POST
  python scripts/laguna/resubmit_hpc.py --api https://qcaihpc-simulation-api.onrender.com
"""
import argparse
import json
import sys
import urllib.request
from pathlib import Path

DEFAULT_API = "https://qcaihpc-simulation-api.onrender.com"


def is_hpc_record(rec):
    """Recognize a solange_hpc.py provenance record.

    NOT a phase check: unlike solange_qpu.py (which stamps "phase" into its
    own local record before ever reaching the backend), solange_hpc.py writes
    NO "phase" field locally — the backend stamps phase=3A-HPC itself, only at
    ingestion (backend/routes/simulate.py's submit_hpc_run). Checking for
    phase=="3A-HPC" here would silently match zero files, always — a real bug
    caught by checking this function's own premise before shipping it, not by
    it ever raising an error. The reliable local signal is p8_hash +
    mutation_id (every P1-P9 record has both), combined with the filename glob
    below (provenance_*.json — solange_qpu.py names its files qpu_*.json, a
    different prefix, so this never picks one of those up by accident either)."""
    return (isinstance(rec, dict) and rec.get("p8_hash") and rec.get("mutation_id")
            and not str(rec.get("phase", "")).startswith("3B-QPU"))


def find_matching_jw(out_dir, provenance_path, rec):
    """solange_hpc.py writes jw_<stem>.json and provenance_<stem>.json as a
    pair with the same stem. The backend's /hpc/submit wants both — an empty
    jw (no --vqe run) is fine, but a real one should be sent if it exists."""
    stem = provenance_path.name
    if stem.startswith("provenance_"):
        jw_path = out_dir / ("jw_" + stem[len("provenance_"):])
        if jw_path.exists():
            try:
                return json.loads(jw_path.read_text())
            except Exception:
                pass
    return {}


def resubmit(api, rec, jw):
    body = json.dumps({"provenance": rec, "jw": jw}, default=str).encode()
    req = urllib.request.Request(
        api.rstrip("/") + "/api/simulate/hpc/submit",
        data=body, method="POST", headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def main():
    ap = argparse.ArgumentParser(description="Re-submit completed classical HPC/DMRG-SCF "
                                             "records to SOLANGE (nothing re-executed).")
    ap.add_argument("--out", default="./out", help="directory holding the saved record JSONs")
    ap.add_argument("--api", default=DEFAULT_API, help="SOLANGE backend base URL")
    ap.add_argument("--dry-run", action="store_true", help="list what would be sent, POST nothing")
    ap.add_argument("--filter", default=None,
                    help="only match provenance filenames containing this substring — an --out "
                         "directory accumulates every run ever done there, and re-submitting all "
                         "of them (even idempotently) is rarely what a specific --submit-less run "
                         "calls for. E.g. --filter cas20-20_20260804 to target one run.")
    args = ap.parse_args()

    out = Path(args.out)
    if not out.is_dir():
        print(f"ERROR: --out directory not found: {out.resolve()}", file=sys.stderr)
        sys.exit(1)

    found = []
    for f in sorted(out.glob("provenance_*.json")):
        if args.filter and args.filter not in f.name:
            continue
        try:
            rec = json.loads(f.read_text())
        except Exception:
            continue            # not a record file (or unreadable) — skip quietly
        if is_hpc_record(rec):
            found.append((f, rec))

    if not found:
        scope = f" matching --filter {args.filter!r}" if args.filter else ""
        print(f"No classical HPC/DMRG-SCF records{scope} found in {out.resolve()}.")
        return

    print(f"Found {len(found)} HPC/DMRG-SCF record(s) in {out.resolve()}:")
    for f, rec in found:
        method = rec.get("orbital_optimization_method") or rec.get("p1_ansatz") or "?"
        print(f"  {f.name}  ·  {rec.get('mutation_id')}/{rec.get('side')}  ·  "
              f"E={rec.get('p7_energy_ha')}  ·  {method}  ·  "
              f"seal={str(rec.get('p8_hash') or '')[:12]}…")

    if args.dry_run:
        print("\n--dry-run: nothing submitted.")
        return

    ok, fail = 0, 0
    for f, rec in found:
        jw = find_matching_jw(out, f, rec)
        try:
            resp = resubmit(args.api, rec, jw)
            status = resp.get("status") or resp.get("db_status") or "?"
            print(f"  SUBMITTED {f.name} -> status={status} run_id={resp.get('run_id')}")
            ok += 1
        except Exception as e:
            print(f"  FAILED {f.name}: {e}", file=sys.stderr)
            fail += 1

    print(f"\n{ok} submitted, {fail} failed.")


if __name__ == "__main__":
    main()
