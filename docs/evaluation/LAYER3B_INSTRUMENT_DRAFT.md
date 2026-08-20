# Layer 3B — Output-Level Evaluation Instrument (DRAFT)

**Status: DRAFT, not pre-registered, not in use.** Nothing here is committed to. This
document exists so that when Phase 4 evaluation begins (Jan–Apr 2027, §05), there is a
starting point rather than a blank page. It will change before it is frozen.

**Do not pre-register this as-is.** Freeze on OSF only immediately before running the
panel, and only after the case set has been re-verified against the codebase at that time
(see "Re-verify before freezing" at the end).

---

## 1 · What this instrument is for

Phase 4's existing layers evaluate the *research* (Hevner G1–G7, Layer 3) and the
*system's integrity properties* (audit trail, reproducibility, Part 11 — Layer 1). Neither
shows an expert a **specific system output** and asks them to judge it.

This layer does that, for one design principle:

> **DP5** — *a result's standing is bounded by its weakest input, and the bound travels
> with the result.*

DP5 is the principle §06.iv describes as derived from failures observed in this artifact
rather than proposed in the abstract. It is implemented and enforced in running code
(`public/Assignment10_Prototype.html`, `STRUCTURALLY_UNRESOLVED`, ~line 4440), which
blocks dispatch rather than only footnoting the limitation. What has never been tested is
whether that enforcement is **legible to someone who did not build it**.

That distinction matters and should be stated plainly in the write-up: this is not an
evaluation of whether the mechanism exists — it demonstrably does — but of whether its
output communicates the bound to a qualified reader unaided.

---

## 2 · Case set

All cases are drawn from live system state, not constructed for the evaluation. The five
refusal cases are exactly the entries of `STRUCTURALLY_UNRESOLVED`; each is blocked from
dispatch in the running platform and carries a written reason shown to the user.

The set covers both independent grounds for refusal:

- **Ground A — a structure exists but does not resolve the mutated residue** (cases R1, R2)
- **Ground B — there is no single mutated residue to build an active space around** (R3, R4, R5)

| ID | Case | Type | Reason shown by the system |
|----|------|------|-----------------------------|
| **P1** | TP53 C275F | positive control | *(not blocked)* measured, structure verified, Class B |
| **R1** | KEAP1 R320Q | refusal — ground A | residue 320 falls 2 residues before the Kelch domain's resolved boundary (322–609, PDB 1U6D/2FLU) — inside a documented flexible IVR–Kelch linker, not the resolved propeller |
| **R2** | STK11 F354L | refusal — ground A | residue 354 falls outside LKB1's resolved kinase domain (43–347, PDB 2WTK) — in a C-terminal region with no identifiable folded domain |
| **R3** | KEAP1 LOF | refusal — ground B | frameshift at residue 196 (pre-Kelch domain) — no single mutant residue to anchor an active space to; the modeled Hamiltonian is borrowed from the unrelated point mutation G333C |
| **R4** | CDKN2A LOF | refusal — ground B | homozygous deletion, not a point mutation — no residue to anchor an active space to |
| **R5** | ARID2 LOF | refusal — ground B | nonsense mutation (p.Gln1118\*) — no single mutant residue to anchor; the cited structure (AlphaFold-Q68CP9) is a computational prediction, not an experimental structure, and the modeled "mutant" compound is backbone-only |
| **D1** | STK11 LKB1 **vs** STK11 D194N | discrimination pair | *(neither blocked)* same gene, same PDB (2WTK): LKB1 is a **size prior** (Class A, 152e/304q, no measurement); D194N is **measured** (Class B, 54e/62q, DMRG CAS(54,31), ΔE=0.001 mHa) |
| **X1** | Formamide shown under NF1, RB1, TP53-LOF, CDKN2A | proxy-disclosure trap | *(not blocked)* one CAS(2e,2o) calculation on a shared model compound, displayed under four targets; the compound name is rendered alongside every number |

**Why D1 and X1 are strong cases.** D1 is the sharpest pair in the repository: two entries
under the same gene, citing the same structure, where one is a measurement and the other
is a size estimate — and the numbers themselves reveal it (304 = 2×152, the Jordan-Wigner
doubling signature of an unmeasured prior; 62 ≠ 2×54, the signature of a real AVAS
selection). X1 reproduces the exact historical failure that produced DP5.

**Note on P1's honest framing.** C275F is Class B *on the current generic active-space
setup*. A chemist-defined setup — orbital selection, and whether a metal ion belongs in the
active space — could move it. Evaluators should be told this; the positive control tests
traceability, not the finality of the class letter.

---

## 3 · Constructs

Six constructs. Two are binary; four are Likert 1–7, parallel to the existing Layer 3
instrument.

| # | Construct | Scale | Applied to |
|---|-----------|-------|-----------|
| C1 | **Refusal correctness** — was withholding the right call? | binary | R1–R5 |
| C2 | **Refusal legibility** — was the reason understandable unaided? | 1–7 | R1–R5 |
| C3 | **Ground discrimination** — is it clear *which kind* of problem this is (no coverage vs. no residue)? | 1–7 | R1–R5 |
| C4 | **Prior-vs-measurement distinguishability** — visible from the interface alone? | 1–7 | D1 |
| C5 | **Proxy disclosure** — recognised as one calculation, not four agreeing results? | binary | X1 |
| C6 | **Evidence traceability** — can you reach the measurement behind the class letter (`bqp_class_source`)? | 1–7 | P1, D1 |

C3 is not in the original proposal. It is included because the two refusal grounds are
genuinely different — one is a coverage failure, the other is a category failure — and a
system that blocks both with reasons the reader cannot tell apart has communicated less
than it appears to.

---

## 4 · Acceptance thresholds (to be fixed at pre-registration)

| Construct | Threshold | Rationale |
|-----------|-----------|-----------|
| **C1 Refusal correctness** | **100%** — every evaluator, every case, binary, no safety net | This is the criterion designed to be *able* to fail. If the system dispatches a case it should have declined, DP5 fails and is reported as failed. Do not soften this after seeing results. |
| C2 Refusal legibility | mean ≥ 5.0, no individual < 4 | matches existing Layer 3 rule |
| C3 Ground discrimination | mean ≥ 5.0, no individual < 4 | same rule; a low score here is a wording finding, not a mechanism failure |
| C4 Prior-vs-measurement | mean ≥ 5.0, no individual < 4 | same rule |
| **C5 Proxy disclosure** | **≥ 75% of evaluators** | Deliberately below the others. This is the hardest trap in the set; an evaluator who falls for it is a finding worth reporting, not a system failure. |
| C6 Evidence traceability | mean ≥ 5.0, no individual < 4 | same rule |

**Inter-rater reliability.** Report Krippendorff's α for C1, C3, C4 — the constructs that
should be objective. α < 0.67 flags a construct needing refinement rather than a system
defect, and should be reported as such.

---

## 5 · Panel and administration

Uses the existing Phase 4 panel (n = 4–6, §05); no additional recruitment.

| Rater | Cases |
|-------|-------|
| Computational chemist | R1, R2, D1, X1 — structural coverage and measurement-vs-prior judgements |
| Regulatory-informed reviewer | R3, R4, R5, P1 — disclosure and evidence-trail judgements |
| All raters | C1 on all five refusal cases (the binary criterion needs full coverage) |

Each case is shown as the platform renders it — the blocked row and its reason, not a
prose summary written for the evaluation. Evaluators are given no explanation beyond what
the interface itself shows; that is precisely what is being measured.

---

## 6 · Re-verify before freezing

This draft was written on 2026-08-20 against the repository as it then stood. Before any
OSF pre-registration, re-verify every case, because a case that no longer matches the
system invalidates the instrument:

1. `STRUCTURALLY_UNRESOLVED` still contains R1–R5, with these reasons
   (`public/Assignment10_Prototype.html`)
2. `targets.json` still shows STK11_LKB1 as a size prior (Class A, 152e/304q) and
   STK11_D194N as measured (Class B, 54e/62q)
3. C275F's class and active-space range are current — they moved during this project's
   history and may move again under a chemist-defined setup
4. The formamide proxy is still shared across the four named genes, and the compound name
   is still rendered alongside the number
5. `python scripts/laguna/verify_consistency.py` returns 0 contradictions

**A known data-layer discrepancy, unresolved at time of writing:** `targets.json` records
`bqp_class: "B"` for ARID2_LOF (and for SMARCA4_R1192C and ARID1A_R1020S) while their own
`bqp_class_source` fields state these are size priors, not measurements, and their
`structure_caveat` fields state the structures are unresolved or invalid. The user-facing
behaviour is correct — ARID2_LOF is blocked in the interface, and the other two are not
surfaced at all — so no evaluator would encounter the discrepancy through the UI. It is
recorded here because an evaluator given direct access to `targets.json` would, and because
deciding whether those entries should return null (as U2AF1_S34F correctly does) is an open
design question, not a settled one.
