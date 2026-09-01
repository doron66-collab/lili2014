#!/usr/bin/env python3
"""
promote_gate1_checks.py — end-of-day merge: gate1_checks (the working cache
for mutations targets.json has never seen) -> targets.json (the permanent,
single source of truth for target facts), then empties the cache.

Why this exists: an NGS report can surface any gene/mutation. Deciding
"does this residue even have a resolvable structure" from scratch every time
a name resurfaces is wasted work once a human has worked it out once — see
backend/routes/gate1.py's own module docstring for the full Gate 1 picture.
gate1_checks (POST /api/gate1/check) is where that one-time verdict lands
during the day; this script is what promotes it into the permanent record so
targets.json never has to be reconciled against a second, silently-drifting
source of the same facts.

Run manually, end of day, reviewing `git diff targets.json` before
committing — the same "verify, don't trust" posture applied to every other
write into targets.json in this project. This script never runs
automatically or on a schedule.

USAGE:
  python3 promote_gate1_checks.py --email you@x.com --password ...
  (or SOLANGE_EMAIL / SOLANGE_PASSWORD env vars, same convention as
  solange_hpc.py's --agent mode)

  --dry-run prints what would change without writing targets.json or
  deleting anything from gate1_checks.
"""
import argparse
import json
import os
import sys
import urllib.request
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_TARGETS_PATH = _REPO_ROOT / "targets.json"

# Same anon key + login flow as solange_hpc.py's --agent mode — this script
# only ever authenticates as a real user (never a service-role key), because
# recording who promoted a verdict matters the same way it matters for any
# other Gate 1/2 sign-off.
_SUPA_URL = "https://lzzuxtnubznrkxwxjaab.supabase.co"
_SUPA_ANON = ("eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imx6"
              "enV4dG51Ynpucmt4d3hqYWFiIiwicm9sZSI6ImFub24iLCJpYXQiOjE3Nzg5MTQ0MzgsImV4cCI6"
              "MjA5NDQ5MDQzOH0.fKApXc3ZPHXh4O008A5oFE5vbTNqJ168AI9NzIl4vHA")


def _supabase_login(email, password):
    req = urllib.request.Request(
        _SUPA_URL + "/auth/v1/token?grant_type=password",
        data=json.dumps({"email": email, "password": password}).encode(),
        method="POST",
        headers={"apikey": _SUPA_ANON, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())["access_token"]


def _api_get(base, path, token):
    req = urllib.request.Request(base + path, headers={"Authorization": "Bearer " + token})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def _api_delete(base, path, token):
    req = urllib.request.Request(base + path, method="DELETE",
                                  headers={"Authorization": "Bearer " + token})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", default=os.environ.get("SOLANGE_EMAIL"),
                     help="login email — or set SOLANGE_EMAIL")
    ap.add_argument("--password", default=os.environ.get("SOLANGE_PASSWORD"),
                     help="login password — or set SOLANGE_PASSWORD")
    ap.add_argument("--submit", default="https://qcaihpc-simulation-api.onrender.com",
                     help="SOLANGE backend base URL")
    ap.add_argument("--dry-run", action="store_true",
                     help="print what would change; write nothing, delete nothing")
    args = ap.parse_args()
    if not args.email or not args.password:
        ap.error("--email/--password required (or SOLANGE_EMAIL/SOLANGE_PASSWORD) — "
                  "same as solange_hpc.py's --agent mode.")

    token = _supabase_login(args.email, args.password)
    pending = _api_get(args.submit, "/api/gate1/list", token).get("records", [])
    if not pending:
        print("gate1_checks is empty — nothing to promote.")
        return

    print(f"{len(pending)} row(s) pending in gate1_checks.\n")
    targets = json.loads(_TARGETS_PATH.read_text())
    mutations = targets.setdefault("mutations", {})
    promoted = []
    for row in pending:
        target = row["target"]
        if mutations.get(target, {}).get("structure_caveat"):
            print(f"  SKIP {target}: targets.json already has a structure_caveat for "
                  f"this key — resolve the conflict by hand, not by overwriting.",
                  file=sys.stderr)
            continue
        verdict = "VERIFIED" if row.get("resolved") else "UNRESOLVED"
        caveat = (f"{verdict} — {row['reason']} (promoted from gate1_checks, "
                  f"checked by {row.get('checked_by', '?')} on {row.get('updated_at', '?')})")
        entry = mutations.setdefault(target, {})
        entry["structure_caveat"] = caveat
        if row.get("pdb_id") and not entry.get("pdb"):
            entry["pdb"] = row["pdb_id"]
        promoted.append(target)
        print(f"  {target}: {verdict} — will merge into targets.json")

    if not promoted:
        print("\nNothing to promote (every row skipped — see above).")
        return

    if args.dry_run:
        print(f"\n[dry run] would write {_TARGETS_PATH} and delete "
              f"{len(promoted)} row(s) from gate1_checks. No changes made.")
        return

    _TARGETS_PATH.write_text(json.dumps(targets, indent=2, ensure_ascii=False) + "\n")
    print(f"\nWrote {_TARGETS_PATH}.")

    for target in promoted:
        _api_delete(args.submit, f"/api/gate1/check/{target}", token)
    print(f"Cleared {len(promoted)} promoted row(s) from gate1_checks.")
    print("\nReview with `git diff targets.json`, then run "
          "`python scripts/laguna/verify_consistency.py`, then commit and push as usual.")


if __name__ == "__main__":
    main()
