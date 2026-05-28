# SOLANGE™ — Artifacts Registry

## Status as of 2026-05-28

| # | Artifact | File | Status | Last Updated | Notes |
|---|---|---|---|---|---|
| 1 | **Live Platform** | `public/Assignment10_Prototype.html` | ✅ Live | 2026-05-28 | Deployed on Netlify from `main`. Multi-user, FDA-compliant, 4-tab interface. |
| 2 | **Backend API** | `backend/routes/simulate.py` | ✅ Live | 2026-05-28 | FastAPI on Render Starter. PennyLane 0.38.0, live VQE circuit execution. |
| 3 | **PennyLane VQE** | `backend/routes/simulate.py` | ✅ Live | 2026-05-28 | 4-qubit AllSinglesDoubles UCCSD, 80 Adam steps, CAS(2e,2o)/STO-3G, converges to 1e-05 Ha of exact. |
| 4 | **JW Hamiltonians** | `backend/jw_hamiltonians.json` | ✅ Complete | 2026-05-15 | Real PySCF CASSCF(2,2) derived. 27 Pauli terms per mutation. 6 mutations × 2 sides. |
| 5 | **NGS Pipeline** | `public/Assignment10_Prototype.html` | ✅ Live | 2026-05-28 | Real patient VCF (MI25-0349). Mutation identification → Hamiltonian selection → VQE. |
| 6 | **FDA Compliance Layer** | `public/Assignment10_Prototype.html` + `backend/routes/simulate.py` | ✅ Live | 2026-05-28 | Full P1–P9 provenance per run. SHA-256 P8 seal. Supabase storage. 21 CFR §11.10(e). |
| 7 | **3D Molecular Visualizer** | `src/NSCLCViewer.tsx` + `src/TP53LoopsViewer.tsx` | ✅ Live | 2026-05-15 | Three.js. Real-time 3D rendering of TP53 C275F active site, loop-sheet-helix, sqDRIFT wavefunction rings. |
| 8 | **VQE Convergence Chart** | `public/Assignment10_Prototype.html` | ✅ Live | 2026-05-28 | Real PennyLane data. SSE streaming — curve builds live during 6-7s computation. HF→E₀ Y-axis. |
| 9 | **Login Notification** | `backend/routes/notify.py` | ✅ Live | 2026-05-26 | Gmail SMTP via smtplib SSL. Notifies doron66@gmail.com on each login. |
| 10 | **Guest Account** | Supabase Auth | ✅ Live | 2026-05-26 | guest@solange.bio / Solange2026. Password change locked. Demo access for IBM/reviewers. |
| 11 | **arXiv Preprint** | `arxiv/main.tex` / `arxiv/main.pdf` | 🔄 Draft | 2026-05-28 | Structured, ready to submit. Phase 3A updated to reflect live PennyLane execution. |
| 12 | **Poster** | `arxiv/IST697_Final_Poster.pptx` | 🔄 Draft | 2026-05-28 | Phase 3A updated: live circuit, 4-qubit, PySCF JW Hamiltonian. PDF regeneration needed. |
| 13 | **Architecture Diagrams** | `public/QC_AI_HPC_Architecture.html` + `public/P1_P9_Architecture_Diagram.html` | ✅ Complete | 2026-05-15 | Three-layer QC-AI-HPC architecture. P1–P9 provenance flow. |
| 14 | **BLAIS 2026 Proposal** | — | 🔄 In Progress | 2026-05-28 | 100-word summary locked. 1,500-word section in progress. Deadline: Aug 15, 2026. |

## Phase Status

| Phase | Description | Status | Hardware |
|---|---|---|---|
| 3A minimal | 4-qubit live VQE, CAS(2e,2o), PennyLane default.qubit | ✅ **Live** | Render server |
| 3A intermediate | 24e/48q local site, CASSCF/DMRG | ⏳ Pending | NSF ACCESS HPC |
| 3B | 44e/88q full active site, sqDRIFT, IBM Heron r3 | ⏳ Pending | IBM Heron r3 |

## Key Facts (for context integrity)

- **SOLANGE™** = Scientific Oncology Legacy Advancing Non-druggable Ground-state Energetics
- **Targets:** TP53, STK11, KEAP1 (NOT CDKN2A in written text)
- **Key mutation:** C275F — 44e / 88q — only target within IBM Heron r3 demonstrated ceiling (Merz et al. 2026)
- **Classical wall:** CCSD(T) breaks down past ~18e
- **Phase 3A convergence:** 1.09e-05 Ha from CASSCF-exact in ~7s on Render Starter
- **IBM contact:** Ella (quantum partnerships, IBM Research Israel) — introduced via Michal Rosen-Zvi
