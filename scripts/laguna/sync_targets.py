#!/usr/bin/env python3
"""
sync_targets.py — push targets.json's values into the frontend tables, so nobody
hand-edits them.

WHY
targets.json is the single source of truth, and the backend reads from it directly.
The frontend cannot: it is a static file served by Netlify, and a runtime fetch would
mean that a network failure leaves the interface with no target data at all — an
unacceptable failure mode for a live demonstration. So the frontend keeps a copy, and
this script makes that copy DERIVED rather than authored.

The frontend was where the contradiction was actually visible: the interface displayed
Class A for the anchor target while the dissertation said B. Leaving those tables
hand-edited would leave the failure surface open while claiming it was closed.

HOW
Surgical, not generative. Only the contested field VALUES are rewritten in place
(bqp_class, active_electrons, full_qubits), matched by key. The surrounding structure,
formatting, and — importantly — the hand-written comments explaining each decision are
left untouched, because that reasoning is worth more than the convenience of emitting
the block wholesale.

    python scripts/laguna/sync_targets.py --check   # report drift, write nothing (CI)
    python scripts/laguna/sync_targets.py           # apply

Exit code is non-zero in --check mode if the frontend has drifted from the source.
"""
import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC  = ROOT / "targets.json"
UI   = ROOT / "public" / "Assignment10_Prototype.html"

FIELDS_NUM = ("active_electrons", "full_qubits")


def _block(text, opener):
    """Return (start, end) of a top-level JS object literal by its opening line."""
    i = text.find(opener)
    if i == -1:
        raise SystemExit(f"could not find '{opener}' in {UI}")
    j = text.find("\n};", i)
    return i, j


def sync(apply: bool):
    src = json.loads(SRC.read_text())
    html = UI.read_text(encoding="utf-8")
    changes = []

    # ── GENE_MAP: one line per gene ───────────────────────────────────────────
    gi, gj = _block(html, "const GENE_MAP = {")
    gene_blk = html[gi:gj]
    new_gene_blk = gene_blk

    for gene, t in src["genes"].items():
        m = re.search(rf"^(\s*{re.escape(gene)}:\s*\{{)([^}}]*)(\}})", new_gene_blk, re.M)
        if not m:
            continue
        body, orig = m.group(2), m.group(2)
        for field in FIELDS_NUM:
            want = t.get(field)
            if want is None:
                continue
            def repl(mm, w=want):
                return f"{mm.group(1)}{w}"
            body, n = re.subn(rf"({re.escape(field)}:\s*)(\d+)", repl, body)
            if n:
                got = re.search(rf"{re.escape(field)}:\s*(\d+)", orig)
                if got and int(got.group(1)) != want:
                    changes.append(f"GENE_MAP.{gene}.{field}: {got.group(1)} -> {want}")
        # bqp_class: present in source -> must be present here; absent -> must be absent
        want_cls = t.get("bqp_class")
        has = re.search(r"bqp_class:\s*'(\w)'", body)
        if want_cls and has and has.group(1) != want_cls:
            changes.append(f"GENE_MAP.{gene}.bqp_class: {has.group(1)} -> {want_cls}")
            body = re.sub(r"(bqp_class:\s*')\w(')", rf"\g<1>{want_cls}\g<2>", body)
        elif want_cls and not has:
            changes.append(f"GENE_MAP.{gene}.bqp_class: MISSING, source says {want_cls}")
        elif has and not want_cls:
            changes.append(f"GENE_MAP.{gene}.bqp_class: present ({has.group(1)}) but source "
                           f"has none — this gene classifies per mutation")
        if body != orig:
            new_gene_blk = new_gene_blk[:m.start(2)] + body + new_gene_blk[m.end(2):]

    # ── MUTATION_CLASS: multi-line entries ────────────────────────────────────
    mi, mj = _block(html, "const MUTATION_CLASS = {")
    mut_blk = html[mi:mj]
    new_mut_blk = mut_blk

    for key, t in src["mutations"].items():
        m = re.search(rf"^(  {re.escape(key)}:\s*\{{)(.*?)(\n  \}})", new_mut_blk, re.S | re.M)
        if not m:
            continue          # not every backend mutation has a frontend entry
        body, orig = m.group(2), m.group(2)
        pairs = [("bqp_class", t.get("bqp_class")),
                 ("active_electrons", t.get("full_electrons")),
                 ("full_qubits", t.get("full_qubits"))]
        for field, want in pairs:
            if want is None:
                continue
            if field == "bqp_class":
                got = re.search(r"bqp_class:\s*'(\w)'", body)
                if got and got.group(1) != want:
                    changes.append(f"MUTATION_CLASS.{key}.bqp_class: {got.group(1)} -> {want}")
                    body = re.sub(r"(bqp_class:\s*')\w(')", rf"\g<1>{want}\g<2>", body)
            else:
                got = re.search(rf"{re.escape(field)}:\s*(\d+)", body)
                if got and int(got.group(1)) != want:
                    changes.append(f"MUTATION_CLASS.{key}.{field}: {got.group(1)} -> {want}")
                    body = re.sub(rf"({re.escape(field)}:\s*)\d+", rf"\g<1>{want}", body)
        if body != orig:
            new_mut_blk = new_mut_blk[:m.start(2)] + body + new_mut_blk[m.end(2):]

    if new_gene_blk != gene_blk or new_mut_blk != mut_blk:
        html = html[:gi] + new_gene_blk + html[gj:]
        # recompute the second block's offsets — the first edit may have shifted them
        mi2, mj2 = _block(html, "const MUTATION_CLASS = {")
        html = html[:mi2] + new_mut_blk + html[mj2:]

    if not changes:
        print("frontend tables already match targets.json — nothing to do.")
        return 0

    print(f"{len(changes)} value(s) drifted from targets.json:")
    for c in changes:
        print(f"  · {c}")
    if apply:
        UI.write_text(html, encoding="utf-8")
        print("\napplied. Run scripts/laguna/verify_consistency.py to confirm.")
        return 0
    print("\n--check mode: nothing written. Run without --check to apply.")
    return 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report drift and exit non-zero; write nothing")
    args = ap.parse_args()
    return sync(apply=not args.check)


if __name__ == "__main__":
    sys.exit(main())
