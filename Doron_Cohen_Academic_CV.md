DORON COHEN
[Street Address, City, State, ZIP]
doron66@gmail.com | [Phone Number] | [LinkedIn / personal site URL]

---

## EDUCATION

**Claremont Graduate University (CGU)**, Claremont, California
PhD / DTech Candidate — *[exact degree title, e.g., "Doctor of Philosophy in Information Systems and Technology"]*
Expected graduation: 2026
Relevant coursework: IST 697
Dissertation Supervisor: Prof. Itamar Shabtai
Dissertation: *[Working title — e.g., "Quantum Simulation Approaches to Non-Druggable Tumor-Suppressor Mutations in Non-Small-Cell Lung Cancer"]*

**[Prior institution]**, [City, State]
[Degree], [Field] — [Year]

---

## RESEARCH FOCUS

Computational and quantum-simulation methods for non-druggable tumor-suppressor mutations in non-small-cell lung cancer (NSCLC), with emphasis on TP53, STK11, and KEAP1. Work bridges classical quantum chemistry (CASSCF/CCSD(T)), variational quantum eigensolver (VQE) methods, and near-term quantum hardware (superconducting qubit devices), with attention to regulatory-grade data provenance and reproducibility.

---

## DISSERTATION PROJECT — SOLANGE™ 3D Quantum Simulation Platform

*Scientific Oncology Legacy Advancing Non-druggable Ground-state Energetics*

Designed and built a full-stack research platform (React/TypeScript frontend, FastAPI/Python backend, Supabase persistence) implementing a two-phase quantum simulation pipeline for oncology target validation:

- **Phase 3A — Classical proxy:** Live VQE ground-state energy calculations using 4-qubit Jordan-Wigner Hamiltonians derived from PySCF CASSCF(2e,2o)/STO-3G computations, addressing the breakdown of classical CCSD(T) methods beyond ~18 active electrons.
- **Phase 3B — Quantum hardware roadmap:** Architecture for full quantum-hardware execution on IBM Heron r3 devices (94+ qubits), incorporating sqDRIFT (sample-based quantum diagonalization), scaling toward the C275F mutation model (44 electrons / 88 qubits).
- **Expansion gene framework:** Extended the simulation engine from a single hardcoded gene to 23 configurable tumor-suppressor/oncogene loss-of-function targets, each independently parameterized (electron count, qubit count, structural source, BQP complexity class, target hardware era) and backed by 29 distinct, independently computed CASSCF(2,2)/STO-3G Jordan-Wigner Hamiltonians built from gene-specific functional-residue model compounds.
- **Data integrity safeguards:** Identified and corrected a placeholder-data defect in which multiple expansion genes shared a single gene's Hamiltonian, and a CASSCF orbital-convergence bug arising from sequential same-process computation; resolved via per-gene subprocess isolation to guarantee independent, scientifically valid results for every gene.
- **Regulatory compliance:** Implemented FDA 21 CFR §11.10(e)-aligned provenance tracking (P1–P9 audit records) with cryptographic audit-hash verification for computational result traceability.

External collaboration: IBM Research (quantum hardware partnership discussions, sqDRIFT sample-based quantum diagonalization algorithm).

---

## GRANTS / AWARDS

- BLAIS Award application, 2026 (CGU) — $10,000–$25,000 — *[Submitted / In Progress]*, addressing CGU strategic priorities in Human Health & Flourishing and Data Analysis & Computational Mathematics.

---

## TECHNICAL SKILLS

- **Quantum / computational chemistry:** Variational Quantum Eigensolver (VQE), Jordan-Wigner transformation, CASSCF, CCSD(T), PySCF
- **Quantum hardware:** IBM Qiskit ecosystem, IBM Heron r3, sample-based quantum diagonalization (sqDRIFT)
- **Software:** Python, TypeScript/React, FastAPI, Supabase/PostgreSQL
- **Deployment:** Netlify (frontend CI/CD), Render (backend CI/CD)
- **Regulatory/compliance:** FDA 21 CFR Part 11 electronic records and audit-trail design

---

## PUBLICATIONS / MANUSCRIPTS

- [Add dissertation chapters, preprints, or papers in progress]

---

## PRESENTATIONS

- [Add conference talks, poster sessions, or guest lectures]

---

## PROFESSIONAL EXPERIENCE

**[Job Title]**, [Organization] — [Dates]
- [Responsibility / accomplishment]

---

## REFERENCES

Available upon request.

---

*Note: This draft was generated from verified facts already documented in the SOLANGE™ project records (CLAUDE.md project memory) — institution, supervisor, course, platform architecture, and technical accomplishments are accurate as of this session. Bracketed fields are placeholders for biographical information (address, phone, prior degrees, employment history, publications, presentations) that should be filled in directly, since I don't have reliable source data for them.*
