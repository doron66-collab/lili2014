#!/usr/bin/env python3
"""
resubmit_qpu.py — restore already-completed QPU runs to SOLANGE WITHOUT touching
the quantum computer again.

Every real QPU run made by solange_qpu.py wrote a full, LEON-sealed P1-P9 record
to its --out directory (default ./out) as a JSON file. Those files are the ground
truth: the P8 seal is already computed, the IBM job already ran. This script finds
the QPU records in that directory and re-POSTs each to the SOLANGE backend, where
LEON re-verifies the seal and (idempotently) upserts the row by its run-id.

  * Zero QPU time is consumed — nothing is re-executed on hardware.
  * Safe to run even if the rows still exist: the backend upserts by id, so a
    re-submit of the same record does not create a duplicate.

USAGE:
  python scripts/laguna/resubmit_qpu.py                 # scan ./out, re-submit QPU rows
  python scripts/laguna/resubmit_qpu.py --out ./out --dry-run   # list, don't POST
  python scripts/laguna/resubmit_qpu.py --api https://qcaihpc-simulation-api.onrender.com
"""
import argparse
import json
import sys
import urllib.request
from pathlib import Path

DEFAULT_API = "https://qcaihpc-simulation-api.onrender.com"


def is_qpu_record(rec):
    """A QPU record is a dict whose phase begins with 3B-QPU (real or dry-run)."""
    return isinstance(rec, dict) and str(rec.get("phase", "")).startswith("3B-QPU")


def resubmit(api, rec):
    body = json.dumps({"provenance": rec, "jw": {}}, default=str).encode()
    req = urllib.request.Request(
        api.rstrip("/") + "/api/simulate/hpc/submit",
        data=body, method="POST", headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def main():
    ap = argparse.ArgumentParser(description="Re-submit completed QPU records to SOLANGE (no QPU time).")
    ap.add_argument("--out", default="./out", help="directory holding the saved QPU record JSONs")
    ap.add_argument("--api", default=DEFAULT_API, help="SOLANGE backend base URL")
    ap.add_argument("--dry-run", action="store_true", help="list what would be sent, POST nothing")
    args = ap.parse_args()

    out = Path(args.out)
    if not out.is_dir():
        print(f"ERROR: --out directory not found: {out.resolve()}", file=sys.stderr)
        sys.exit(1)

    found = []
    for f in sorted(out.glob("*.json")):
        try:
            rec = json.loads(f.read_text())
        except Exception:
            continue            # not a record file (or unreadable) — skip quietly
        if is_qpu_record(rec):
            found.append((f, rec))

    if not found:
        print(f"No QPU records (phase 3B-QPU*) found in {out.resolve()}.")
        print("If ./out was cleared, restore instead with:  solange_qpu.py --retrieve <JOB_ID> "
              "--key <KEY> --side <side> --submit   (also zero QPU time).")
        return

    print(f"Found {len(found)} QPU record(s) in {out.resolve()}:")
    for f, rec in found:
        tgt = rec.get("mutation_id") or rec.get("mutation_name") or "?"
        print(f"  · {f.name:42s} target={tgt:16s} phase={rec.get('phase')} "
              f"backend={rec.get('p3_backend')} p8={str(rec.get('p8_hash'))[:10]}…")

    if args.dry_run:
        print("\n--dry-run: nothing sent. Re-run without --dry-run to restore these rows.")
        return

    print()
    ok = 0
    for f, rec in found:
        try:
            resp = resubmit(args.api, rec)
            state = resp.get("db_status")
            print(f"  {f.name}: seal_ok={resp.get('seal_ok')} db={state} run_id={resp.get('run_id')}")
            if state in ("stored", "stored_no_payload"):
                ok += 1
        except Exception as e:
            print(f"  {f.name}: SUBMIT FAILED — {e}", file=sys.stderr)
    print(f"\nDone: {ok}/{len(found)} stored. Reload SOLANGE → Orchestration → Rung 4 (Phase 3B QPU).")


if __name__ == "__main__":
    main()
