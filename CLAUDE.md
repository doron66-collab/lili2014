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
  range. Eight real DMRG classifications exist now, ALL Class B: C275F, R282W, G245S,
  R249S, SETD2 R1625C, KEAP1 G333C, STK11 D194N, R175H. No Class A candidate currently
  exists among measured targets (see R175H finding below — its earlier provisional Class A
  reading did not survive a warm-started re-run).
- R175H methodological finding (2026-08-15): three independent cold-start DMRG runs showed
  large non-convergence (ΔE≈1.4–1.8 Ha at M=250→500) — correctly read, at the time, as
  provisional Class A per the classifier's own non-convergence rule. A later warm-started
  re-run (loading MPS state accumulated across those earlier HPC-ticket-interrupted attempts,
  via the --scratch resume fix in solange_dmrg.py) reached a materially better energy already
  at its nominal M=250 step, and converged cleanly (ΔE=0.006 mHa, then 0.011 mHa under the
  identical M=250,500,1000 schedule the cold-start runs used, ruling out a too-close-M-values
  artifact) — both give S_max≈0.78. Conclusion: the earlier large ΔE was a cold-start
  optimization artifact (10 sweeps from a random state, insufficient to reach the true
  minimum at low M for a 54-orbital active space), not genuine high entanglement. Final:
  Class B, S_max=0.78. Worth remembering for any future large-active-space DMRG run: a
  cold-start bond-dimension convergence test can give a false non-convergence signal; warm-
  starting (or more sweeps per M) is needed to tell real entanglement growth apart from
  optimization failure at low M.
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
- As of 2026-08-15, all eight real DMRG classifications (C275F, R282W, G245S, R249S,
  SETD2 R1625C, KEAP1 G333C, STK11 D194N, R175H) are Class B — classically tractable at
  their measured ranges. None currently carries a quantum-necessity case; every target
  tested so far is quantum-advantaged, not quantum-necessary. This does NOT block r3 access
  — see the three-tier value framing below, where operational scale-up (tier 1) and real
  provenance-at-scale (tier 2) both stand on their own regardless of any target's class —
  but it does mean there is presently no confirmed candidate for tier 3 (a quantum-necessity
  demonstration). A future chemist-defined active space, or a not-yet-tested gene (KEAP1's
  gene-level generic case, ARID1A pending a valid structure), could still surface one; none
  has yet.
- FRAMING for the dissertation: r3 is a STRENGTHENING of Phase 3B, NOT a prerequisite. The
  contribution is the ARCHITECTURE (already substantiated). Three value tiers, in order:
  (1) operational scale-up on real hardware [safest — proves the governed pipeline holds at
  larger circuit width; directly answers RQ II]; (2) real P3/P4 provenance/telemetry at scale
  [strengthens the Part-11 claim §06.iii]; (3) a quantum-necessity demonstration via SqDRIFT,
  once a Class A candidate exists [upside/risky, no confirmed candidate as of 2026-08-15 —
  DO NOT bet the thesis on tier 3]. Also: separate publication + BLAIS grant
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
