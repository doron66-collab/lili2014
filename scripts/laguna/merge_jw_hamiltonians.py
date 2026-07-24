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

USAGE (from the repo root, after generating the source files with solange_hpc.py):
  python scripts/laguna/merge_jw_hamiltonians.py ./out
  python scripts/laguna/merge_jw_hamiltonians.py ./out/jw_KEAP1_LOF_native_*.json  (glob also works)

Nothing here talks to SOLANGE or spends any compute/quantum time — it only
edits the two local JSON files. Review with `git diff` before committing.
"""
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
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    sources = find_source_files(sys.argv[1:])
    if not sources:
        print("No jw_*.json source files found in the given path(s).", file=sys.stderr)
        sys.exit(1)

    # Collect {key: {side: entry}} from every source file first, so each of the
    # two hamiltonian files is only opened/written once even with many sources.
    merged = {}
    for f in sources:
        try:
            data = json.loads(f.read_text())
        except Exception as e:
            print(f"  skip {f.name}: could not parse ({e})", file=sys.stderr)
            continue
        for key, sides in data.items():
            merged.setdefault(key, {}).update(sides)
        print(f"  read {f.name}: {list(data.keys())}")

    if not merged:
        print("Nothing valid to merge.", file=sys.stderr)
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
        print(f"WROTE {target.relative_to(_REPO_ROOT)} — merged: {', '.join(added)}")

    print("\nDone. Review with `git diff jw_hamiltonians.json backend/jw_hamiltonians.json`,")
    print("then commit and push as usual.")


if __name__ == "__main__":
    main()
