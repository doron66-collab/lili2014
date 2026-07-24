#!/usr/bin/env python3
"""
merge_jw_hamiltonians.py — merge solange_hpc.py's generated jw_<KEY>_*.json files
into BOTH copies of jw_hamiltonians.json (root, used by solange_qpu.py by default;
backend/jw_hamiltonians.json, used by the backend's JW-download endpoint) so a
newly-generated target becomes runnable on real quantum hardware.

Each jw_<KEY>_<side>_...json file solange_hpc.py writes has the shape
{ "<KEY>": { "<side>": {...entry...} } } — this script merges each entry into
the two hamiltonian files by key/side (adds/overwrites just that key+side,
leaves every other entry untouched), then writes both files back
(indent=2, sorted keys) so the diff stays clean and reviewable.

SAFETY (learned the hard way — an earlier version of this script accepted a
whole directory and silently pulled in unrelated STALE files sitting in ./out
from other runs, overwriting a working, QPU-ready TP53_C275F Hamiltonian with
a leftover CAS(8,8) classical-HPC file that had ZERO Pauli terms):

  1. --keys is REQUIRED — an explicit allowlist of the top-level key(s) you
     intend to touch (e.g. --keys KEAP1_LOF,STK11_LKB1). Any source file whose
     top-level key is not in this list is skipped and reported, never merged —
     so a stray file for an unrelated target can never slip in even if it's
     sitting in the same directory.
  2. Any entry with zero/missing Pauli terms ("terms" empty or absent) is
     REFUSED — never written, even if its key is in the allowlist. A QPU
     Hamiltonian with no terms is not a smaller/different result, it's broken.

USAGE (from the repo root, after generating the source files with solange_hpc.py):
  python scripts/laguna/merge_jw_hamiltonians.py --keys KEAP1_LOF,STK11_LKB1 ./out
  python scripts/laguna/merge_jw_hamiltonians.py --keys KEAP1_LOF ./out/jw_KEAP1_LOF_native_*.json

Nothing here talks to SOLANGE or spends any compute/quantum time — it only
edits the two local JSON files. Review with `git diff` before committing.
"""
import argparse
import glob
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_TARGETS = [_REPO_ROOT / "jw_hamiltonians.json", _REPO_ROOT / "backend" / "jw_hamiltonians.json"]


def find_source_files(args):
    files = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            files += sorted(p.glob("jw_*.json"))
        else:
            files += sorted(Path(x) for x in glob.glob(a))
    return [f for f in files if f.exists()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="source file(s) or a directory of jw_*.json files")
    ap.add_argument("--keys", required=True,
                    help="comma-separated allowlist of top-level keys you intend to merge "
                         "(e.g. KEAP1_LOF,STK11_LKB1) — anything else found is skipped, never merged")
    args = ap.parse_args()
    allowlist = {k.strip() for k in args.keys.split(",") if k.strip()}

    sources = find_source_files(args.paths)
    if not sources:
        print("No jw_*.json source files found in the given path(s).", file=sys.stderr)
        sys.exit(1)

    merged = {}
    skipped_off_list, skipped_empty = [], []
    for f in sources:
        try:
            data = json.loads(f.read_text())
        except Exception as e:
            print(f"  skip {f.name}: could not parse ({e})", file=sys.stderr)
            continue
        for key, sides in data.items():
            if key not in allowlist:
                skipped_off_list.append(f"{f.name} ({key})")
                continue
            for side, entry in sides.items():
                terms = entry.get("terms")
                if not terms:
                    skipped_empty.append(f"{f.name} ({key}/{side})")
                    continue
                merged.setdefault(key, {})[side] = entry
                print(f"  read {f.name}: {key}/{side} — {len(terms)} Pauli terms")

    if skipped_off_list:
        print(f"\nSkipped (key not in --keys allowlist, left untouched): {', '.join(skipped_off_list)}")
    if skipped_empty:
        print(f"REFUSED (zero/missing Pauli terms — would corrupt the entry): {', '.join(skipped_empty)}",
              file=sys.stderr)

    if not merged:
        print("\nNothing valid to merge (see skip/refuse reasons above).", file=sys.stderr)
        sys.exit(1)

    for target in _TARGETS:
        existing = json.loads(target.read_text()) if target.exists() else {}
        added = []
        for key, sides in merged.items():
            existing.setdefault(key, {})
            for side, entry in sides.items():
                is_new = side not in existing[key]
                existing[key][side] = entry
                added.append(f"{key}/{side}" + ("" if is_new else " (overwrote existing)"))
        target.write_text(json.dumps(existing, indent=2, sort_keys=True) + "\n")
        print(f"\nWROTE {target.relative_to(_REPO_ROOT)} — merged: {', '.join(added)}")

    print("\nDone. Review with `git diff jw_hamiltonians.json backend/jw_hamiltonians.json`,")
    print("then commit and push as usual.")


if __name__ == "__main__":
    main()
