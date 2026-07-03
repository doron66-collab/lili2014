# SOLANGE™ — Artifacts Registry

**Status as of 2026-07-03.** Complete inventory of all project artifacts by category.
Canonical source files are tracked in git (repo `doron66-collab/lili2014`, feature
branch `claude/code-access-clarification-ab1W8`, deploys from `main`).

---

## A. Dissertation & Academic Writing

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| A1 | **Dissertation (main defense doc)** | `dissertation_revised.html` (+ `public/` copy) | HTML | ✅ Current | IST 697 DTech proposal. All 8 committee fixes applied (SqDRIFT, VQE↔SqDRIFT, BQP/QMA, consistency, venue, Part-11 stochasticity, novelty anchors, links). A/B/C taxonomy aligned to live system. |
| A2 | **arXiv preprint** | `arxiv/main.tex`, `arxiv/main.pdf` | LaTeX/PDF | 🔄 Draft | SqDRIFT + venue + BQP/QMA fixes applied. Needs `pdflatex` recompile before submission. |
| A3 | **Conference poster** | `arxiv/IST697_Final_Poster.pptx` | PowerPoint | 🔄 Draft | SqDRIFT + venue corrected. PDF export still needed. |
| A4 | **Theoretical Foundations** (lit review, TOGAF lens) | `theoretical_foundations.md`, `.html` (+ `public/`) | MD/HTML | ✅ Current | Venn synthesis, ADM cycle, §2.3 basis-set scope. SqDRIFT + evolution framing applied. |
| A5 | **Verified references** | `references_verified.md` | MD | ✅ Current | Reference list with verification notes. |
| A6 | **Academic CV** | `Doron_Cohen_Academic_CV.md` | MD | ✅ Current | SqDRIFT wording corrected. |

## B. Positioning / Advisor Documents (companion — not dissertation body)

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| B1 | **Value Proposition Memo** (to Prof. Shabtai → Molly) | `SOLANGE_Value_Proposition_Memo.md`, `.docx` | MD/Word | ✅ Current | 6 sections. §6.3 revised per Itamar's outcome-based feedback. SqDRIFT wording corrected. |
| B2 | **Value Proposition Notes** (working talking-points) | `SOLANGE_value_proposition_notes.md` | MD | ✅ Current | Internal notes behind B1. |
| B3 | **Quantum Hardware Approaches** (plain-language guide for Molly) | `Quantum_Computing_Hardware_Approaches_updated.docx` | Word | ✅ Current | 5-modality guide. Physical:logical ratios updated to late-2025 (Helios ~2:1, QuEra ~4.7:1); Sources section; Infleqtion added; memory-vs-computation caveat. |

## C. Live Platform — Frontend

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| C1 | **Main platform UI** | `public/Assignment10_Prototype.html` | HTML/JS | ✅ Live | Netlify (deploys from `main`). 4-tab interface, guest login, SSE VQE convergence chart, P1–P9 provenance display. Hosts the Non-Druggable Target Registry (A/B/C queue). |
| C2 | **React 3D app** | `src/App.tsx`, `NSCLCViewer.tsx`, `TP53LoopsViewer.tsx`, `PDBMolViewer.tsx`, `IntroScreen.tsx`, `DataPanel.tsx`, `main.tsx` | TypeScript/React | ✅ Live | Three.js + NGL. 3D active-site rendering; reads uploaded NGS variants via localStorage bridge; PDB/AlphaFold fallback. |

## D. Live Platform — Backend

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| D1 | **FastAPI backend** | `backend/main.py`, `backend/routes/` (`simulate`, `admin`, `notify`, `pdb`, `provenance`, `__init__`) | Python | ✅ Live | Render (Starter). CORS, 5 route modules. |
| D2 | **VQE engine + P1–P9 provenance** | `backend/routes/simulate.py` | Python | ✅ Live | PennyLane 0.38 `default.qubit`, 4-qubit UCCSD, single-run SSE (consolidated), user_id attribution, Phase 3A/3B boundary enforcement. |
| D3 | **User provisioning** | `backend/create_users.py` | Python | ✅ Live | Supabase Auth account creation (researcher/admin roles). |
| D4 | **DB schema + migration** | `backend/supabase_schema.sql`, `backend/migrations/2026-06-05_add_p2_p5_provenance_columns.sql` | SQL | ✅ Live | `simulation_runs`, `provenance_audit`, `users_profile`; RLS policies. |
| D5 | **Backend config** | `backend/render.yaml`, `requirements.txt`, `.env.example` | Config | ✅ Live | Deploy + dependency + env template. |

## E. Scientific Data & Compute Scripts

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| E1 | **JW Hamiltonians** | `backend/jw_hamiltonians.json`, `jw_hamiltonians.json` (root) | JSON | ✅ Current | Real PySCF CASSCF(2,2) JW Pauli terms. All 48 entries reconciled for CAS(2,2) integrity. |
| E2 | **CASSCF datasets** | `all_mutations_casscf.json`, `c275f_casscf_phase3a.json`, `pauli_coefficients.json` | JSON | ✅ Complete | Precomputed classical reference data. |
| E3 | **PySCF / CASSCF scripts** | `pyscf_c275f.py`, `pyscf_mutations.py`, `casscf_analysis.py`, `generate_expansion_jw.py`, `process_glossary.py` | Python | ✅ Complete | Hamiltonian generation + analysis pipeline. |

## F. Architecture Diagrams & Supporting HTML

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| F1 | **QC·AI·HPC Architecture** | `public/QC_AI_HPC_Architecture.html` | HTML | ✅ Complete | Three-layer architecture diagram. |
| F2 | **P1–P9 Architecture Diagram** | `public/P1_P9_Architecture_Diagram.html` | HTML | ✅ Complete | Provenance flow. |
| F3 | **Three-layer architecture** | `public/three_layer_architecture.html` | HTML | ✅ Current | Venue corrected (IBM Research 2026). |
| F4 | **Unified architecture** | `public/unified_architecture.html` | HTML | ✅ Complete | — |
| F5 | **Platform architecture provenance** | `public/platform_architecture_provenance.html` | HTML | ✅ Complete | — |
| F6 | **System Architecture Report** | `public/System_Architecture_Report.html` | HTML | ✅ Complete | — |
| F7 | **Figures / logos** | `public/figures/all_mutations_combined.svg`, `tp53_c275f.svg`, `public/logos/cgu-flame.png` | SVG/PNG | ✅ Complete | — |

## G. Build / Project Config

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| G1 | **Vite/React build config** | `package.json`, `tsconfig.json`, `netlify.toml`, `index.html` | Config | ✅ Live | Frontend build (deploys `public/` + `dist/`). |

## H. Project Memory / Registry

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| H1 | **Project memory** | `CLAUDE.md` | MD | ✅ Current | SOLANGE facts, git rules, SqDRIFT definition corrected. |
| H2 | **This registry** | `ARTIFACTS.md` | MD | ✅ Current | Artifact inventory (this file). |

## I. External / Not in Repo

| # | Artifact | Location | Type | Notes |
|---|---|---|---|---|
| I1 | **Committee review** | uploaded (not committed) | Word | Prof. Shabtai committee-review doc; drove the 8 dissertation fixes. |
| I2 | **BLAIS 2026 proposal** | — | — | 100-word summary locked; 1,500-word section in progress. Deadline Aug 15, 2026. |
| I3 | **Deployed site** | `https://solange-platform.bio` | Live URL | Netlify. Guest: guest@solange.bio / Solange2026. |

---

## ⚠ Maintenance Notes

- **Root ↔ `public/` HTML duplication.** Several files exist in BOTH the repo root
  and `public/` (`dissertation_revised.html`, `theoretical_foundations.html`,
  `Assignment10_Prototype.html`, `QC_AI_HPC_Architecture.html`,
  `P1_P9_Architecture_Diagram.html`). **`public/` is what Netlify serves** — always
  edit and keep both copies in sync (dissertation + theoretical_foundations are
  kept synced by the current workflow). Root `index.html` and
  `preview_apple_design.html` are entry/preview pages.
- **`dist/`** is a gitignored build output — stale copies live there; ignore, it
  rebuilds.
- **arXiv PDF (A2)** and **poster PDF (A3)** need regeneration from their sources
  before external use.
