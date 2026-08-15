#!/usr/bin/env python3
"""
resubmit_dmrg.py — restore already-completed DMRG classifier runs to SOLANGE
WITHOUT re-running them.

Same idea as resubmit_hpc.py (for solange_hpc.py/DMRG-SCF runs) and
resubmit_qpu.py (for real quantum-hardware runs), generalized to
solange_dmrg.py: any run NOT started with --submit still wrote a full,
sealed dmrg_class_<key>.json to its --out directory. That record is
complete and already sealed (dmrg_hash) — so getting it onto the platform
is a matter of POSTing the same JSON solange_dmrg.py's own --submit would
have sent, not recomputing anything.

This is exactly the situation tonight's TP53 mutation sweep (R282W, G245S,
R249S, ...) created: every run used --out ./out but none used --submit, so
the platform's live "measured" DMRG points (window.__dmrgRows on the
Classical Tractability Plane chart) never got them.

  * Nothing is re-executed — this only re-POSTs what is already on disk.
  * Safe to run even if the row already exists: the backend upserts by id.

USAGE:
  python scripts/laguna/resubmit_dmrg.py                 # scan ./out, re-submit DMRG rows
  python scripts/laguna/resubmit_dmrg.py --out ./out --dry-run   # list, don't POST
  python scripts/laguna/resubmit_dmrg.py --api https://qcaihpc-simulation-api.onrender.com
"""
import argparse
import json
import sys
import urllib.request
from pathlib import Path

DEFAULT_API = "https://qcaihpc-simulation-api.onrender.com"


def is_dmrg_record(rec):
    """Recognize a solange_dmrg.py output record: it always carries dmrg_hash
    (the seal) and bqp_class (the classifier's verdict) together, and
    dmrg_energies (the per-M ladder) — a combination no other record type in
    this pipeline produces."""
    return (isinstance(rec, dict) and rec.get("dmrg_hash") and rec.get("bqp_class")
            and "dmrg_energies" in rec)


def resubmit(api, rec):
    body = json.dumps(rec, default=str).encode()
    req = urllib.request.Request(
        api.rstrip("/") + "/api/simulate/hpc/dmrg/submit",
        data=body, method="POST", headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def main():
    ap = argparse.ArgumentParser(description="Re-submit completed DMRG classifier "
                                             "records to SOLANGE (nothing re-executed).")
    ap.add_argument("--out", default="./out", help="directory holding dmrg_class_*.json")
    ap.add_argument("--api", default=DEFAULT_API, help="SOLANGE backend base URL")
    ap.add_argument("--dry-run", action="store_true", help="list what would be sent, POST nothing")
    ap.add_argument("--filter", default=None,
                    help="only match filenames containing this substring — e.g. "
                         "--filter R282W to target one target's record")
    args = ap.parse_args()

    out = Path(args.out)
    if not out.is_dir():
        print(f"ERROR: --out directory not found: {out.resolve()}", file=sys.stderr)
        sys.exit(1)

    found = []
    for f in sorted(out.glob("dmrg_class_*.json")):
        if args.filter and args.filter not in f.name:
            continue
        try:
            rec = json.loads(f.read_text())
        except Exception:
            continue            # not a record file (or unreadable) — skip quietly
        if is_dmrg_record(rec):
            found.append((f, rec))

    if not found:
        scope = f" matching --filter {args.filter!r}" if args.filter else ""
        print(f"No DMRG classifier records{scope} found in {out.resolve()}.")
        return

    print(f"Found {len(found)} DMRG classifier record(s) in {out.resolve()}:")
    for f, rec in found:
        print(f"  {f.name}  ·  {rec.get('key')}/{rec.get('side')}  ·  "
              f"CAS({rec.get('nelecas')},{rec.get('ncas')})  ·  "
              f"S_max={rec.get('s_max')}  ·  class={rec.get('bqp_class')}  ·  "
              f"seal={str(rec.get('dmrg_hash') or '')[:12]}…")

    if args.dry_run:
        print("\n--dry-run: nothing submitted.")
        return

    ok, fail = 0, 0
    for f, rec in found:
        try:
            resp = resubmit(args.api, rec)
            status = resp.get("status") or resp.get("db_status") or "?"
            print(f"  SUBMITTED {f.name} -> status={status} run_id={resp.get('run_id')}")
            ok += 1
        except Exception as e:
            print(f"  FAILED {f.name}: {e}", file=sys.stderr)
            fail += 1

    print(f"\n{ok} submitted, {fail} failed.")


if __name__ == "__main__":
    main()
