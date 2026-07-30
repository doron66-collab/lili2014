# SOLANGE™ — Project Memory

## Acronym
**SOLANGE™** = Scientific Oncology Legacy Advancing Non-druggable Ground-state Energetics

## Platform
- Full name: SOLANGE™ 3D Quantum Simulation Platform
- URL: https://solange-platform.bio
- Guest login: guest@solange.bio / Solange2026

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
- Key mutation: C275F (44 electrons / 88 qubits)
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
- Feasibility of the 44e/88q anchor (TP53 C275F): the blocker was ACCESS, not feasibility.
  88q fits in Heron r3 (~156q). NOT via naive VQE (NISQ noise kills a deep 88q circuit),
  but via **SqDRIFT / sample-based quantum diagonalization** — QPU samples dominant configs,
  classical diagonalizes the sampled subspace (noise-robust). IBM already demoed 32e/100q
  (C13Cl2, arXiv:2603.08696); 44e/88q is comparable → plausibly in reach, NOT guaranteed,
  and hard to validate (no exact classical ground truth at 44e → compare vs DMRG/CCSD(T)).
- FRAMING for the dissertation: r3 is a STRENGTHENING of Phase 3B, NOT a prerequisite. The
  contribution is the ARCHITECTURE (already substantiated). Three value tiers, in order:
  (1) operational scale-up on real hardware [safest — proves the governed pipeline holds at
  larger circuit width; directly answers RQ II]; (2) real P3/P4 provenance/telemetry at scale
  [strengthens the Part-11 claim §06.iii]; (3) the 44e anchor attempt via SqDRIFT [upside/risky].
  DO NOT bet the thesis on tier 3. Also: separate publication + BLAIS grant strength.

## LEON
- **LEON** = Lineage-Evidence Orchestration & Notarization — the single notarization authority
- Canonical module: `backend/routes/leon.py` (owns the P8 seal: build_p8_payload / build_p8_seal / notarize / reverify)
- `simulate.py` (/hpc/submit) and `provenance.py` (/runs/{id}/verify) BOTH delegate to leon — one source of truth, seal never drifts
- Principle: "verify, don't trust" (DP1, §06.iv). Named in memory of Doron's father, Leon
- UI: results show "notarized by LEON"; dissertation has dedication + §06.ii LEON paragraph + §06.iv + glossary entry

## Standing Tasks
- **Dissertation sync (ALWAYS):** every change to the SOLANGE system must be reflected in the dissertation (§06.ii LEON, §06.iii audit-in-running-code, §06.iv DP1–DP4). Treat this as a permanent, non-optional step of any feature/UI change. **Single file only:** `public/dissertation_revised.html` is the one real file (it's what Netlify actually serves); root `dissertation_revised.html` is a symlink to it, kept only for convenience — never recreate it as a second real file (two copies drifted apart once already, from a parallel session editing only one side; the symlink makes that impossible now).
- **UI ladder order:** Orchestration tab blocks are ordered small→large→quantum (laptop/in-browser VQE at top → HPC classical → DMRG → QPU at bottom); keep new blocks in that ladder.

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
