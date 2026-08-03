#!/usr/bin/env python3
"""
verify_consistency.py — cross-check every factual claim the platform makes about
its targets against every other place that makes it.

WHY THIS EXISTS
A contradiction was found by accident: the dissertation states "C275F is Class B"
while GENE_MAP and simulate.py both say "A". Finding the next one by accident is
not a strategy. This script mechanically compares the four places that assert
facts about targets and reports every disagreement, so the question stops being
"has anyone noticed?" and becomes "does the check pass?".

SOURCES COMPARED
  1. public/Assignment10_Prototype.html  — GENE_MAP (per gene)
  2. backend/routes/simulate.py          — MUTATION_CONFIGS (per mutation)
  3. public/dissertation_revised.html    — the prose being defended
  4. jw_hamiltonians.json                — the stored Hamiltonians

The dissertation is treated as authoritative where it speaks: it is the artifact
under examination, and code that contradicts it is wrong by definition — a defence
cannot rest on "the implementation disagrees with the thesis".

WHAT IT CANNOT DO
Only mechanically comparable claims are checked — numbers and class letters. A
claim the dissertation makes in prose but nowhere in code, or a piece of reasoning
that is simply wrong, is out of reach here. Passing this is necessary, not sufficient.

    python scripts/laguna/verify_consistency.py

Exit code is non-zero if any contradiction is found, so it can gate CI.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
UI   = ROOT / "public" / "Assignment10_Prototype.html"
API  = ROOT / "backend" / "routes" / "simulate.py"
DISS = ROOT / "public" / "dissertation_revised.html"
JW   = ROOT / "jw_hamiltonians.json"


def gene_map():
    """GENE_MAP entries: gene -> {active_electrons, full_qubits, bqp_class}."""
    s = UI.read_text(encoding="utf-8")
    blk = s[s.find("const GENE_MAP = {"):]
    blk = blk[:blk.find("\n};")]
    out = {}
    for m in re.finditer(r"^\s*(\w+):\s*\{([^}]*)\}", blk, re.M):
        g, body = m.group(1), m.group(2)
        ae = re.search(r"active_electrons:\s*(\d+)", body)
        fq = re.search(r"full_qubits:\s*(\d+)", body)
        bc = re.search(r"bqp_class:\s*'(\w)'", body)
        out[g] = {"active_electrons": int(ae.group(1)) if ae else None,
                  "full_qubits": int(fq.group(1)) if fq else None,
                  "bqp_class": bc.group(1) if bc else None}
    return out


def mutation_class():
    """MUTATION_CLASS entries: key -> {bqp_class, active_electrons, full_qubits}.

    Some facts are per mutation, not per gene — TP53's C275F and Y220C differ in both
    class and active-space size. Where this table speaks it OVERRIDES GENE_MAP, and the
    checker must compare against it rather than flagging the gene-level default as a
    contradiction of a value it was never meant to describe."""
    s = UI.read_text(encoding="utf-8")
    blk = s[s.find("const MUTATION_CLASS = {"):]
    blk = blk[:blk.find("\n};")]
    out = {}
    for m in re.finditer(r"^\s*(\w+):\s*\{(.*?)\n  \}", blk, re.S | re.M):
        k, body = m.group(1), m.group(2)
        bc = re.search(r"bqp_class:\s*'(\w)'", body)
        ae = re.search(r"active_electrons:\s*(\d+)", body)
        fq = re.search(r"full_qubits:\s*(\d+)", body)
        out[k] = {"bqp_class": bc.group(1) if bc else None,
                  "active_electrons": int(ae.group(1)) if ae else None,
                  "full_qubits": int(fq.group(1)) if fq else None}
    return out


def mutation_configs():
    """MUTATION_CONFIGS entries: key -> {full_electrons, full_qubits, bqp_class, jw_source}."""
    s = API.read_text(encoding="utf-8")
    blk = s[s.find("MUTATION_CONFIGS"):]
    out = {}
    for m in re.finditer(r'"(\w+)":\s*\{(.*?)\n    \}', blk, re.S):
        k, body = m.group(1), m.group(2)
        if "bqp_class" not in body:
            continue
        fe = re.search(r'"full_electrons":\s*(\d+)', body)
        fq = re.search(r'"full_qubits":\s*(\d+)', body)
        bc = re.search(r'"bqp_class":\s*"(\w)"', body)
        js = re.search(r'"jw_source":\s*\("(\w+)",\s*"(\w+)"\)', body)
        out[k] = {"full_electrons": int(fe.group(1)) if fe else None,
                  "full_qubits": int(fq.group(1)) if fq else None,
                  "bqp_class": bc.group(1) if bc else None,
                  "jw_source": (js.group(1), js.group(2)) if js else None}
    return out


def dissertation_claims():
    """Explicit, machine-checkable assertions in the prose."""
    s = DISS.read_text(encoding="utf-8")
    claims = {"classes": {}, "electrons": {}}
    # "C275F is Class B"
    for m in re.finditer(r"\b([A-Za-z0-9]+) is Class ([ABC])\b", s):
        claims["classes"][m.group(1)] = m.group(2)
    # "the three KEAP1 target mutations (...), classified Class B"
    for m in re.finditer(r"\b([A-Z][A-Z0-9]+)(?:/\w+)? target mutations \([^)]*\), classified Class ([ABC])", s):
        claims["classes"].setdefault(m.group(1), m.group(2))
    # "~155 active electrons", "~152 active electrons"
    for m in re.finditer(r"\b([A-Z][A-Z0-9]{2,})\b[^<]{0,160}?~(\d+) active electrons", s):
        claims["electrons"].setdefault(m.group(1), int(m.group(2)))
    return claims


def main():
    gm, mc, dc = gene_map(), mutation_configs(), dissertation_claims()
    mcls = mutation_class()
    jw = json.loads(JW.read_text())
    problems, checked = [], 0

    print("=" * 78)
    print("SOLANGE — cross-source consistency check")
    print("=" * 78)

    # ── 1. Class letters: code vs the dissertation ────────────────────────────
    print("\n[1] BQP class — code vs dissertation")
    for name, cls in sorted(dc["classes"].items()):
        gene = re.match(r"[A-Z]+[0-9]*[A-Z]*", name).group(0)
        gene = gene if gene in gm else next((g for g in gm if name.startswith(g)), None)
        for src, val in (("GENE_MAP", gm.get(gene, {}).get("bqp_class") if gene else None),
                         ("simulate.py", next((v["bqp_class"] for k, v in mc.items()
                                               if name.upper() in k.upper()), None))):
            if val is None:
                continue
            checked += 1
            if val != cls:
                problems.append(f"{name}: dissertation says Class {cls}, {src} says Class {val}")
                print(f"    MISMATCH  {name:12} dissertation={cls}  {src}={val}")
            else:
                print(f"    ok        {name:12} dissertation={cls}  {src}={val}")

    # ── 2. Class letters: GENE_MAP vs simulate.py ─────────────────────────────
    print("\n[2] BQP class — GENE_MAP vs simulate.py")
    for key, cfg in sorted(mc.items()):
        gene = key.split("_")[0]
        # Mutation-level table wins where it speaks — that is the whole point of it.
        ui_cls, src = (mcls[key]["bqp_class"], "MUTATION_CLASS") if key in mcls \
                      else (gm.get(gene, {}).get("bqp_class"), f"GENE_MAP[{gene}]")
        if ui_cls is None or cfg["bqp_class"] is None:
            continue
        checked += 1
        if ui_cls != cfg["bqp_class"]:
            problems.append(f"{key}: {src}={ui_cls} vs simulate.py={cfg['bqp_class']}")
            print(f"    MISMATCH  {key:14} {src}={ui_cls}  simulate.py={cfg['bqp_class']}")
        else:
            print(f"    ok        {key:14} both {cfg['bqp_class']}  (via {src})")

    # ── 3. Electron counts: GENE_MAP vs simulate.py ───────────────────────────
    print("\n[3] Active-electron count — GENE_MAP vs simulate.py")
    for key, cfg in sorted(mc.items()):
        gene = key.split("_")[0]
        ui_e, src = (mcls[key]["active_electrons"], "MUTATION_CLASS") if key in mcls \
                    else (gm.get(gene, {}).get("active_electrons"), f"GENE_MAP[{gene}]")
        if ui_e is None or cfg["full_electrons"] is None:
            continue
        checked += 1
        if ui_e != cfg["full_electrons"]:
            problems.append(f"{key}: {src}={ui_e}e vs simulate.py full_electrons={cfg['full_electrons']}e")
            print(f"    MISMATCH  {key:14} {src}={ui_e}e  simulate.py={cfg['full_electrons']}e")
        else:
            print(f"    ok        {key:14} both {cfg['full_electrons']}e  (via {src})")

    # ── 4. Jordan-Wigner arithmetic: qubits must be 2x electrons ──────────────
    print("\n[4] Jordan-Wigner arithmetic — qubits == 2 x electrons")
    for label, src in (("GENE_MAP", [(g, v["active_electrons"], v["full_qubits"]) for g, v in gm.items()]),
                       ("simulate.py", [(k, v["full_electrons"], v["full_qubits"]) for k, v in mc.items()])):
        for name, e, q in src:
            if e is None or q is None:
                continue
            checked += 1
            if q != 2 * e:
                problems.append(f"{label} {name}: {e}e should map to {2*e}q under JW, but says {q}q")
                print(f"    MISMATCH  {label} {name:14} {e}e -> stated {q}q, expected {2*e}q")

    # ── 5. jw_source keys must exist in the Hamiltonian file ──────────────────
    print("\n[5] simulate.py jw_source -> jw_hamiltonians.json")
    for key, cfg in sorted(mc.items()):
        if not cfg["jw_source"]:
            continue
        k, side = cfg["jw_source"]
        checked += 1
        if k not in jw or side not in (jw.get(k) or {}):
            problems.append(f"{key}: jw_source ({k}, {side}) missing from jw_hamiltonians.json")
            print(f"    MISSING   {key:14} -> {k}/{side}")
        else:
            print(f"    ok        {key:14} -> {k}/{side} ({jw[k][side].get('compound')})")

    print("\n" + "=" * 78)
    print(f"{checked} comparisons made · {len(problems)} contradiction(s) found")
    print("=" * 78)
    for p in problems:
        print(f"  ✗ {p}")
    if problems:
        print("\nThe dissertation is authoritative where it speaks: it is the artifact under\n"
              "examination, and code contradicting it is wrong by definition.")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
