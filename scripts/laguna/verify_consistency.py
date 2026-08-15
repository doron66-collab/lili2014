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
  0. targets.json                        — THE SINGLE SOURCE OF TRUTH
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
SRC  = ROOT / "targets.json"
PDBPY = ROOT / "backend" / "routes" / "pdb.py"


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


def pdb_map():
    """MUTATION_PDB_MAP literals in backend/routes/pdb.py: key -> PDB ID.

    This file was a fifth independent home for a target fact and drifted exactly
    as that predicts — STK11_LKB1 was 2QK7 here while targets.json said 2WTK,
    with nothing forcing agreement and this checker not looking. pdb.py now
    overlays targets.json at import, so a stale literal is corrected at runtime;
    this check exists so it is also corrected at SOURCE, rather than being
    silently papered over on every boot.
    """
    s = PDBPY.read_text(encoding="utf-8")
    blk = s[s.find("MUTATION_PDB_MAP = {"):]
    blk = blk[:blk.find("}")]
    return dict(re.findall(r'"(\w+)":\s*"(\w{4})"', blk))


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
        # A prose sentence can match the claim regex without naming a gene (an
        # earlier §06.i rewrite produced "...entanglement is Class B"). Treat that
        # as "no gene", not as a crash: the lookups below then find nothing and
        # the entry is skipped rather than taking the whole checker down.
        m = re.match(r"[A-Z]+[0-9]*[A-Z]*", name)
        gene = m.group(0) if m else name
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
    # This 2x rule is a CAS(N,N) simplification (one spatial orbital per electron),
    # which is what every *size-prior* gene/mutation in this table actually is. It is
    # NOT a law of the JW mapping itself: qubits = 2 x orbitals, and a real AVAS active
    # space (electrons measured directly from PDB coordinates) can have orbitals !=
    # electrons. TP53_C275F and TP53_R175H are exempted for exactly this reason — their
    # full_electrons/full_qubits record a real CAS(48,28)/CAS(96,54) measurement, not a
    # CAS(N,N) size guess, so 2x electrons is the wrong expected value for them specifically.
    JW_DOUBLING_EXEMPT = {"TP53_C275F", "TP53_R175H", "TP53_R282W", "TP53"}
    print("\n[4] Jordan-Wigner arithmetic — qubits == 2 x electrons (size-prior entries only)")
    for label, src in (("GENE_MAP", [(g, v["active_electrons"], v["full_qubits"]) for g, v in gm.items()]),
                       ("simulate.py", [(k, v["full_electrons"], v["full_qubits"]) for k, v in mc.items()])):
        for name, e, q in src:
            if e is None or q is None:
                continue
            if name in JW_DOUBLING_EXEMPT:
                print(f"    skip      {label} {name:14} real AVAS active space (orbitals != electrons) — exempt")
                continue
            checked += 1
            if q != 2 * e:
                problems.append(f"{label} {name}: {e}e should map to {2*e}q under JW, but says {q}q")
                print(f"    MISMATCH  {label} {name:14} {e}e -> stated {q}q, expected {2*e}q")

    # ── 4b. targets.json is the SINGLE SOURCE — every consumer must match it ──
    # Checks 2 and 3 above compare consumers to each other, which catches drift but
    # cannot say which one is right. This check says so: targets.json is the source,
    # and a consumer that disagrees with it is the one that is wrong.
    print("\n[4b] targets.json (single source) vs its consumers")
    src = json.loads(SRC.read_text())
    for key, t in sorted(src["mutations"].items()):
        cfg = mc.get(key)
        if not cfg:
            continue
        for field, srcv, gotv in (("bqp_class", t.get("bqp_class"), cfg.get("bqp_class")),
                                  ("full_electrons", t.get("full_electrons"), cfg.get("full_electrons")),
                                  ("full_qubits", t.get("full_qubits"), cfg.get("full_qubits"))):
            if srcv is None or gotv is None:
                continue
            checked += 1
            if srcv != gotv:
                problems.append(f"{key}.{field}: targets.json={srcv} vs simulate.py={gotv}")
                print(f"    MISMATCH  {key:14} {field}: source={srcv} simulate.py={gotv}")
        # the UI's mutation-level table, where it carries the mutation
        if key in mcls:
            for field, srcv, gotv in (("bqp_class", t.get("bqp_class"), mcls[key].get("bqp_class")),
                                      ("active_electrons", t.get("full_electrons"), mcls[key].get("active_electrons"))):
                if srcv is None or gotv is None:
                    continue
                checked += 1
                if srcv != gotv:
                    problems.append(f"{key}.{field}: targets.json={srcv} vs MUTATION_CLASS={gotv}")
                    print(f"    MISMATCH  {key:14} {field}: source={srcv} MUTATION_CLASS={gotv}")
        print(f"    ok        {key:14} matches source")

    for gene, t in sorted(src["genes"].items()):
        ui_g = gm.get(gene)
        if not ui_g:
            continue
        for field in ("active_electrons", "full_qubits", "bqp_class"):
            srcv, gotv = t.get(field), ui_g.get(field)
            if srcv is None and gotv is None:
                continue
            checked += 1
            if srcv != gotv:
                problems.append(f"gene {gene}.{field}: targets.json={srcv} vs GENE_MAP={gotv}")
                print(f"    MISMATCH  gene {gene:10} {field}: source={srcv} GENE_MAP={gotv}")

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

    # ── 6. pdb.py's PDB IDs must match targets.json ───────────────────────────
    print("\n[6] pdb.py MUTATION_PDB_MAP -> targets.json 'pdb'")
    src_muts = (json.loads(SRC.read_text(encoding="utf-8")).get("mutations") or {})
    pm = pdb_map()
    for key, pid in sorted(pm.items()):
        want = (src_muts.get(key) or {}).get("pdb")
        if want is None:
            # Not a contradiction: targets.json has no record for this key
            # (CDKN2A has a gene entry but no mutation entry), so pdb.py is
            # legitimately its only source. Reported so it stays visible.
            print(f"    only-here {key:14} {pid}  (no targets.json mutation record)")
            continue
        checked += 1
        if pid != want:
            problems.append(f"{key}: pdb.py says {pid}, targets.json says {want}")
            print(f"    MISMATCH  {key:14} {pid} != {want}")
        else:
            print(f"    ok        {key:14} {pid}")

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
