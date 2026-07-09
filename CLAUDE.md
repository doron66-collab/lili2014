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
- Targets: TP53, STK11, KEAP1 (NOT CDKN2A in written text)
- Key mutation: C275F (44 electrons / 88 qubits)
- Classical limit: CCSD(T) breaks down past ~18e
- Phase 3A: Live classical proxy — VQE ground-state energies, 4-qubit JW Hamiltonian, PySCF CAS(2e,2o)/STO-3G
- Phase 3B: Full quantum hardware — IBM Heron r3, sqDRIFT sample-based quantum diagonalization, 94+ qubits

## IBM
- Algorithm: sqDRIFT (sample-based quantum diagonalization; quantum sampling + classical diagonalization, NOT noise characterization)
- Hardware: IBM Heron r3
- Contact: Michal Rosen-Zvi (ROSEN@il.ibm.com) — leaving IBM, referred to Ella (quantum partnerships)
- Ella: leads quantum partnerships at IBM Research Israel

## LEON
- **LEON** = Lineage-Evidence Orchestration & Notarization — the single notarization authority
- Canonical module: `backend/routes/leon.py` (owns the P8 seal: build_p8_payload / build_p8_seal / notarize / reverify)
- `simulate.py` (/hpc/submit) and `provenance.py` (/runs/{id}/verify) BOTH delegate to leon — one source of truth, seal never drifts
- Principle: "verify, don't trust" (DP1, §06.iv). Named in memory of Doron's father, Leon
- UI: results show "notarized by LEON"; dissertation has dedication + §06.ii LEON paragraph + §06.iv + glossary entry

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
