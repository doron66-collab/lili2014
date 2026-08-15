# SOLANGE™ — Project Memory

## Communication
- Doron communicates in Hebrew. Write replies in plain right-to-left Hebrew; when an
  English term is needed, place it in parentheses immediately after its Hebrew
  equivalent (e.g. "מוטציה (mutation)", "תחום הליבה (core domain)") — never as a bare
  Latin token embedded mid-sentence, which breaks bidi (RTL/LTR) rendering and garbles
  the sentence on his screen. This applies to every conversation, not just one session.

## Acronym
**SOLANGE™** = Scientific Oncology Legacy Advancing Non-druggable Ground-state Energetics

## Platform
- Full name: SOLANGE™ 3D Quantum Simulation Platform
- URL: https://solange-platform.bio
- No guest login. Retired 2026-08-03 — a shared, publicly-known password (it had
  been repeated in this file and in chat many times) could inherit whatever NGS
  report was last loaded in the same browser via localStorage, and could reach
  the destructive Reset control. Every account is now individually provisioned
  via Supabase Auth ("Add user"), so there is always a real user_id behind every
  session. Three tiers only: admin, executive, full user — see showPlatform() in
  Assignment10_Prototype.html.

## Dissertation
- Candidate: Doron Cohen
- Program: PhD/DTech
- University: Claremont Graduate University (CGU), California
- Course: IST 697
- Supervisor: Prof. Itamar Shabtai
- Year: CGU 2026

## Science
- Domain: Non-druggable NSCLC tumor-suppressor mutations
- Targets: four — TP53, KEAP1, CDKN2A, STK11 (CDKN2A confirmed as the fourth written target, 2026-07-26; the earlier "exclude CDKN2A from written text" rule is retired — the dissertation names all four consistently)
- Key mutation: C275F — 18–48 active electrons / 24–56 qubits, AVAS-tested range (criterion-
  dependent, not a single value); DMRG-measured Class B at 36e (S_max=0.25) and 48e
  (S_max=0.35), both real, converged, FINAL results. The "44e/88q" figure used throughout
  the dissertation and codebase until 2026-08-14/15 was an architectural estimate that
  traced back to no verifiable derivation — replaced everywhere with the real measured
  range. Five real DMRG classifications exist now (all Class B): C275F, R282W, G245S,
  R249S, SETD2 R1625C. R175H is a provisional Class A (real, reproduced non-convergence
  at M=250→500 across three independent runs), pending a completed bond-dimension ladder.
- Classical limit: CCSD(T) breaks down past ~18e
- Phase 3A: Live classical proxy — VQE ground-state energies, 4-qubit JW Hamiltonian, PySCF CAS(2e,2o)/STO-3G
- Phase 3B: Full quantum hardware — IBM Heron r3, sqDRIFT sample-based quantum diagonalization, 94+ qubits

## IBM
- Algorithm: sqDRIFT (sample-based quantum diagonalization; quantum sampling + classical diagonalization, NOT noise characterization)
- Hardware: IBM Heron r3
- Contact: Michal Rosen-Zvi (ROSEN@il.ibm.com) — leaving IBM, referred to Ella (quantum partnerships)
- Ella: leads quantum partnerships at IBM Research Israel

## Heron r3 access — strategy (in negotiation, 2026-07; leave until access lands)
- STATUS: Doron negotiating with the university to fund r3 runs. Not yet available.
- The anchor target (TP53 C275F) is now DMRG-measured Class B (S_max=0.25 at 36e, 0.35 at
  48e, both real and converged) — classically tractable at its measured range, so it is
  quantum-advantaged, not quantum-necessary. It does NOT carry the quantum-necessity case
  for r3 access on its own anymore. That case now rests on TP53 R175H: real, reproduced
  DMRG non-convergence at M=250→500 (three independent runs) up to 108 qubits — provisional
  Class A, the actual candidate for "why real hardware is needed," not C275F.
  108q still fits Heron r3 (~156q). NOT via naive VQE (NISQ noise kills a deep circuit this
  size), but via **SqDRIFT / sample-based quantum diagonalization** — QPU samples dominant
  configs, classical diagonalizes the sampled subspace (noise-robust). IBM already demoed
  32e/100q (C13Cl2, arXiv:2603.08696); R175H's range is comparable in scale → plausibly in
  reach, NOT guaranteed, and hard to validate (no exact classical ground truth at this scale
  — compare vs DMRG at higher bond dimension once affordable).
- FRAMING for the dissertation: r3 is a STRENGTHENING of Phase 3B, NOT a prerequisite. The
  contribution is the ARCHITECTURE (already substantiated). Three value tiers, in order:
  (1) operational scale-up on real hardware [safest — proves the governed pipeline holds at
  larger circuit width; directly answers RQ II]; (2) real P3/P4 provenance/telemetry at scale
  [strengthens the Part-11 claim §06.iii]; (3) the R175H quantum-necessity attempt via SqDRIFT
  [upside/risky — DO NOT bet the thesis on tier 3]. Also: separate publication + BLAIS grant
  strength.

## LEON
- **LEON** = Lineage-Evidence Orchestration & Notarization — the single notarization authority
- Canonical module: `backend/routes/leon.py` (owns the P8 seal: build_p8_payload / build_p8_seal / notarize / reverify)
- `simulate.py` (/hpc/submit) and `provenance.py` (/runs/{id}/verify) BOTH delegate to leon — one source of truth, seal never drifts
- Principle: "verify, don't trust" (DP1, §06.iv). Named in memory of Doron's father, Leon
- UI: results show "notarized by LEON"; dissertation has dedication + §06.ii LEON paragraph + §06.iv + glossary entry

## TODO — session vs. history separation (design, not started)
- Proposed 2026-08-15: separate what the end-user dynamic list shows (only targets that came
  up in the currently-loaded NGS report; cleared on Reset) from a permanent, never-cleared
  history store (every run ever submitted to the platform, across its entire lifetime — not
  the admin tab). On Reset, the current-session results page is cleared but its records still
  get archived into the history view, not deleted. Motivated by a real gap found this session:
  DMRG classifier submissions (e.g. SETD2_R1625C, PDB 5JJY) land in the live "DMRG A/B/C
  Classifications" table regardless of any NGS session, while the NGS-driven dynamic gene list
  maps genes through a completely separate, generic GENE_MAP entry (e.g. SETD2 → PDB 5JLB, a
  different structure/compound) that has no link to that specific mutation record — two
  unconnected sources of truth, the same shape of problem as the 44e/88q drift this session
  spent hours fixing elsewhere. Needs a real session_id/is_current distinction at the data
  layer, not a UI-only fix — scope properly before implementing, not mid-session.

## Standing Tasks
- **Dissertation sync (ALWAYS):** every change to the SOLANGE system must be reflected in the dissertation (§06.ii LEON, §06.iii audit-in-running-code, §06.iv DP1–DP4). Treat this as a permanent, non-optional step of any feature/UI change. **Single file only:** `public/dissertation_revised.html` is the one real file (it's what Netlify actually serves); root `dissertation_revised.html` is a symlink to it, kept only for convenience — never recreate it as a second real file (two copies drifted apart once already, from a parallel session editing only one side; the symlink makes that impossible now).
- **UI ladder order:** Orchestration tab blocks are ordered small→large→quantum (laptop/in-browser VQE at top → HPC classical → DMRG → QPU at bottom); keep new blocks in that ladder.
- **Consistency check (ALWAYS, after any change to target facts):** run
  `python scripts/laguna/verify_consistency.py`. Two seconds, no environment needed —
  it only reads four repo files. Expect `0 contradiction(s) found`; non-zero exit means
  a real contradiction. This exists because target facts are hardcoded in **four
  independent places** (GENE_MAP ~45 genes, MUTATION_CLASS, backend MUTATION_CONFIGS,
  and the dissertation prose), so nothing structurally forces them to agree. That is how
  "C275F is Class B" in the dissertation sat next to `bqp_class:'A'` in both the frontend
  and the backend — the demo overclaiming, on the anchor target, relative to the document
  it demonstrates. Found by accident while reading; the checker exists so the next one
  is not. **Necessary, not sufficient:** it compares numbers and class letters only. A
  claim made in one place only, or reasoning that is simply wrong, is out of its reach.
- **Which table owns a class:** genes the dissertation records as having *"no
  single-residue anchor"* (KEAP1, CDKN2A, STK11) carry `bqp_class` on GENE_MAP — for them
  the class is a property of the domain, so gene-level is correct. Genes whose mutations
  genuinely differ (TP53: C275F = B at 44e, Y220C = C at 38e) carry **no** gene-level
  class and live in `MUTATION_CLASS`. Never add a gene-level class to a gene that
  classifies per mutation: an unlisted variant must come back *not classified*, never a
  size proxy. Every entry needs its dissertation quote in `source`.
- **Benchmarks:** `python scripts/laguna/verify_benchmark.py` re-derives every reference
  energy from its own stated specification. A reference nobody can regenerate is an
  assertion, not a reference.
- **Single source of truth — `targets.json`.** Target facts (bqp_class, electron and
  qubit counts) are edited **there and nowhere else**. The backend overlays it onto
  `MUTATION_CONFIGS` at import; the frontend's `GENE_MAP`/`MUTATION_CLASS` are **derived**,
  kept in sync by `python scripts/laguna/sync_targets.py` (add `--check` for a read-only
  CI gate). Never hand-edit those tables — a hand-edit is exactly how the frontend came to
  display Class A while the dissertation said B. After changing `targets.json`: run
  `sync_targets.py`, then `verify_consistency.py`. Every class needs its dissertation
  quote in `bqp_class_source`.

## Evidence discipline (non-negotiable)
- **Never present a proxy value as a target property.** Platform energies are computed
  over CAS(2e,2o)/STO-3G **model compounds** — toluene for a Phe sidechain, methanethiol
  for Cys — and compounds are **shared across genes**, so formamide's number appears
  identically under NF1, RB1, TP53-LOF and CDKN2A. Always name the compound alongside the
  number, or four views of one calculation read as four independent results agreeing.
- **Scope gates everything.** A measurement is authoritative only for the active space it
  covered: `solangeDmrgCoversSite` (UI) and `_covers_site` (gateway) both require ≥90% of
  the target's site. A CAS(6,4) run does not classify a 44e site — it annotates, it never
  overrides.
- **Confidence is capped by provenance, not by presence** (`backend/routes/gateway.py`).
  Unstated provenance is treated as unverified, never assumed good.
- **Grades are earned from inputs, never declared** (`backend/routes/evidence.py`). No
  target currently supports a chemical or biological finding; `POC_DISCLAIMER` is the one
  canonical wording — mirror it, don't rewrite it.
- This is an **Information Systems** dissertation. The architecture is the artifact; the
  chemistry is not. Demonstrating the pipeline is the claim — never that a computed value
  is chemically correct.

## Compliance
- Standard: FDA 21 CFR §11.10(e)
- Provenance: P1–P9 records
- Storage: Supabase (service_role key required for inserts)
- Audit hash: P8, truncated to 16 chars with copy pill

## Tech Stack
- Frontend: React + TypeScript (Netlify, deploys from main)
- Backend: FastAPI Python (Render, deploys from main)
- DB: Supabase
- Email: Gmail SMTP via smtplib SSL, GMAIL_APP_PASSWORD env var
- 3D intro: IntroScreen.tsx

## Git
- ALWAYS push to BOTH branches:
  - git push origin <branch>
  - git push origin <branch>:main
- Feature branch: claude/code-access-clarification-ab1W8
- Netlify and Render both watch: main

## BLAIS 2026
- Award range: $10,000–$25,000
- Deadline: August 15, 2026
- Submit to: Eusebio.Alvaro@cgu.edu
- 100-word summary: locked
- 1,500-word section: in progress
- CGU strategic priorities addressed: Human Health & Flourishing, Data Analysis & Computational Mathematics
- Undergraduate partner faculty: TBD (Prof. Shabtai working on it)

## CSS Variables (Assignment10_Prototype.html)
- --white: #f1f5f9
- --gray: #cbd5e1
- --navy: #0B1E3D
- --teal: #06b6d4
- --green: #22c55e
- Dark text fix: always use var(--white) or var(--gray), never #223244 or #475569
