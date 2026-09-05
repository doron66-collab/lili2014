# SOLANGE™ — Artifacts Registry

**Status as of 2026-09-04.** Complete inventory of all project artifacts by category.
Canonical source files are tracked in git (repo `doron66-collab/lili2014`, feature
branch `claude/code-access-clarification-ab1W8`, deploys from `main`).

---

## A. Dissertation & Academic Writing

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| A1 | **Dissertation (main defense doc)** | `public/dissertation_revised.html` | HTML | ✅ Current | The one real file — root `dissertation_revised.html` is a **symlink** to it, kept only for convenience (never a second real copy — a prior parallel-edit drift is exactly why it's a symlink now). Full line-by-line accuracy pass completed 2026-09-04: Gate 2 dispatch, SHCI mandatory two-method policy + real Class A cross-validation finding, the Rung 3→4 QPU escalation path, JW encoder generalization, a `targets.json` structure-caveat drift, a fabricated citation, and several TOC/reference/glossary mechanical errors all found and fixed. |
| A2 | **arXiv preprint** | `arxiv/main.tex`, `arxiv/main.pdf` | LaTeX/PDF | 🔄 Draft | Not re-verified against the dissertation's 2026-09-04 pass — check for drift before submission. |
| A3 | **Conference poster** | `arxiv/IST697_Final_Poster.pptx` | PowerPoint | ✅ Current | Already reflected the SHCI cross-check and R175H's corrected Class B verdict as of 2026-09-04; one fabricated citation ("Merz et al. 2026", shared with a since-fixed dissertation error) removed same day. PDF export still needed for external use. |
| A4 | **Theoretical Foundations** (lit review, TOGAF lens) | `theoretical_foundations.md`, `.html` (+ `public/` copy) | MD/HTML | ✅ Current | Root and `public/` copies confirmed byte-identical 2026-09-04 (manually kept in sync, not a symlink — unlike A1). |
| A5 | **Verified references** | `references_verified.md` | MD | ⚠ Not reviewed | Not touched in the 2026-09-04 accuracy pass — the dissertation's own reference list was checked directly instead. |
| A6 | **Academic CV** | `Doron_Cohen_Academic_CV.md` | MD | ✅ Current | Stale SOLANGE acronym ("...Advancing Non-druggable Ground-state Energetics") found and corrected 2026-09-04. |

## B. Positioning / Advisor Documents (companion — not dissertation body)

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| B1 | **Value Proposition Memo** (to Prof. Shabtai → Molly) | `SOLANGE_Value_Proposition_Memo.md`, `.docx` | MD/Word | ⚠ Not reviewed since 2026-07-03 | Not re-checked in the 2026-09-04 pass. |
| B2 | **Value Proposition Notes** (working talking-points) | `SOLANGE_value_proposition_notes.md` | MD | ⚠ Not reviewed since 2026-07-03 | Internal notes behind B1. |
| B3 | **Quantum Hardware Approaches** (plain-language guide) | `Quantum_Computing_Hardware_Approaches_updated.docx` | Word | ⚠ Not reviewed since 2026-07-03 | 5-modality guide. |
| B4 | **Production Deployment Model** | `PRODUCTION_DEPLOYMENT_MODEL.md` | MD | ⚠ Not reviewed | New since the last registry snapshot — not yet cross-checked against §08 of the dissertation, which covers the same ground. |

## C. Live Platform — Frontend

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| C1 | **Main platform UI** | `public/Assignment10_Prototype.html` | HTML/JS | ✅ Live | Netlify (deploys from `main`). The single largest live file in the repo — now hosts the full Rung 1–4 escalation ladder (laptop N2 self-test → HPC classical → DMRG/SHCI dual-method classifier → real QPU dispatch), Gate 1 (structural resolvability) and Gate 2 (mechanism-category sign-off + a real pre-dispatch block), the NGS upload/parse pipeline (PDF via pdf.js + OCR fallback via Tesseract.js, VCF, JSON — all routed through a human-review form before loading), the HPC dispatch queue view, and the Classical Tractability Plane chart. **Guest login was retired 2026-08-03** (a shared password could inherit a prior session's loaded NGS report via localStorage and reach the destructive Reset control) — every account is now individually provisioned via Supabase Auth; do not publish shared credentials for this platform anywhere. |
| C2 | **React 3D app** | `src/App.tsx`, `NSCLCViewer.tsx`, `TP53LoopsViewer.tsx`, `PDBMolViewer.tsx`, `IntroScreen.tsx`, `DataPanel.tsx`, `main.tsx` | TypeScript/React | ⚠ Not reviewed since 2026-07-03 | Three.js + NGL. Not touched in recent sessions' work; verify it still reads the current NGS-variant bridge format before relying on it in a demo. |

## D. Live Platform — Backend

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| D1 | **FastAPI backend** | `backend/main.py`, `backend/routes/` | Python | ✅ Live | Render. Now **11** route modules (up from 5): `simulate`, `admin`, `notify`, `pdb`, `provenance`, `gate1`, `gate2`, `gateway`, `leon`, `evidence`, `security_log`. |
| D2 | **VQE / HPC dispatch + P1–P9 provenance** | `backend/routes/simulate.py` | Python | ✅ Live | No longer a single 4-qubit demo path: `/hpc/dispatch` now branches on `job_type` (`dmrg`, `shci`, `qpu`, `screen_classify`) and on whether a real `geometry`/`avas` payload is present (Gate-2-cleared custom compounds, or a real classified target's own active space carried forward to Rung 4). SHCI cross-validation (`/hpc/shci/submit`, `/hpc/shci/list`) computes `class_agreement` (same A/B/C verdict) as a field distinct from `agreement` (chemical-accuracy energy match) — see A1's dissertation note on why conflating them misreads a real corroboration case as a contradiction. |
| D3 | **Structural resolvability gate** | `backend/routes/gate1.py` | Python | ✅ Live | `STRUCTURALLY_UNRESOLVED` refusal pattern — blocks a target with no residue to anchor an active space to. |
| D4 | **Mechanism-category gate** | `backend/routes/gate2.py` | Python | ✅ Live | `dispatch_custom_compound` — a real pre-dispatch gate (not only a taxonomy): saves geometry/category/sign-offs, computes `missing_requirements`, and refuses dispatch until every sign-off the target's mechanism category requires is on record. The category-specific *consumer* (a calculator that branches its computed quantity on the recorded category) remains unbuilt — every dispatch that clears this gate still runs the same ground-state calculation regardless of category. |
| D5 | **Computational Decision Gateway** | `backend/routes/gateway.py` | Python | ✅ Live | Evidence-based routing recommendation; `SITE_COVERAGE_MIN = 0.9` scope gate (`_covers_site`) — a measurement only counts if it covers ≥90% of a target's active site. |
| D6 | **LEON (notarization authority)** | `backend/routes/leon.py` | Python | ✅ Live | Single source of truth for the P8 seal (`build_p8_payload` / `build_p8_seal` / `notarize` / `reverify`); `simulate.py` and `provenance.py` both delegate to it rather than each computing their own seal. |
| D7 | **Evidence grading** | `backend/routes/evidence.py` | Python | ✅ Live | No target currently supports a chemical/biological finding above the platform's own disclaimer grade; grades are earned from inputs, never declared. |
| D8 | **User provisioning** | `backend/create_users.py` | Python | ✅ Live | Supabase Auth account creation — now the *only* way to get platform access, per C1's guest-login retirement. |
| D9 | **DB schema + migrations** | `backend/supabase_schema.sql`, `backend/migrations/` | SQL | ✅ Live | Has grown well beyond the original `simulation_runs`/`provenance_audit`/`users_profile` set — now also `hpc_dispatch`, `dmrg_classifications`, `shci_crossvalidations` (with `class_agreement`), `gate1_checks`, `gate2_records`, `custom_compounds`, `leon_audit`. No single migration file enumerates all of them; check each route module's own inline migration-SQL comments. |
| D10 | **Backend config** | `backend/render.yaml`, `requirements.txt`, `.env.example` | Config | ✅ Live | |

## E. Scientific Data & Compute Scripts

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| E1 | **Legacy fixed-demo JW Hamiltonians** | `backend/jw_hamiltonians.json`, `jw_hamiltonians.json` (root) | JSON | ✅ Current, superseded for new work | Real PySCF CASSCF(2,2) JW Pauli terms for the original 48-entry fixed compound library (Rung 1 / Phase 3A demo, and the 6 real QPU runs of Table 7.2). Still exercised by the N2 pipeline self-test; no longer the only path to a JW Hamiltonian — see E3. |
| E2 | **CASSCF datasets** | `all_mutations_casscf.json`, `c275f_casscf_phase3a.json`, `pauli_coefficients.json` | JSON | ✅ Complete | Precomputed classical reference data for the fixed compound library. |
| E3 | **JW term generator (generalized)** | `generate_expansion_jw.py` | Python | ✅ Current | `build_jw_terms()` generalized from a hardcoded `n_orb=2` to an arbitrary active-space size (verified against PySCF's own FCI solver to <1e-6 Ha at 4 orbitals; exercised without error at 8). This is what makes `geometry_target()` (E5) possible — a real classified target's own active space, not just the fixed library, can now become a QPU-ready Hamiltonian. |
| E4 | **PySCF / CASSCF scripts (legacy)** | `pyscf_c275f.py`, `pyscf_mutations.py`, `casscf_analysis.py`, `process_glossary.py` | Python | ✅ Complete | Original Hamiltonian generation + analysis pipeline for the fixed compound library. |
| E5 | **Classical rung — DMRG classifier** | `scripts/laguna/solange_dmrg.py` | Python | ✅ Live | Both `--compound` (fixed library) and `--geometry`+`--avas` (real PDB-carved cluster) modes. `--scratch` warm-start/resume. Orbital optimization now performed by the DMRG sweep itself (tensor-network solver replacing the CASSCF full-CI step), exercised at CAS(20,20)/40 qubits. For a real `--geometry` run, stores the raw `.xyz` **content** in its record (not the local file path — fixed 2026-09-02 after tracing why a downstream QPU-escalation attempt would have written a garbage geometry file). |
| E6 | **Classical rung — SHCI cross-validator** | `scripts/laguna/solange_shci.py` | Python | ✅ Live | Independent classification via Selected CI (Dice solver) — mandatory pairing with DMRG as of 2026-09-02 (not merely an optional second opinion). SOSCF (Newton) RHF fallback for borderline SCF convergence. Now also captures `geometry`/`avas`/`charge`/`spin` in its submitted record (added 2026-09-02 — previously never captured at all, which is why the SHCI table's own QPU-escalation button never appeared). |
| E7 | **Quantum rung — QPU dispatch agent** | `scripts/laguna/solange_qpu.py` | Python | ✅ Live | `geometry_target()` (new 2026-09-02): builds a QPU-ready Hamiltonian directly from a classified target's own geometry/AVAS/charge/spin, via E3's generalized JW encoder — closing the Rung 3→4 data path. **Architectural only, not yet exercised against real hardware**: no target measured to date carries a Class A verdict on a real (non-proxy) biological site, and `run_agent`'s `measure(..., hardware=True)` has no automatic dry-run gate. |
| E8 | **HPC orchestration agent** | `scripts/laguna/solange_hpc.py` | Python | ✅ Live | Pull-based dispatch loop (`--agent`) claiming `hpc_dispatch` rows over an outbound-only channel; branches on `job_type` and on whether `job.geometry` is present. Live per-stage progress reporting (RHF → CASSCF → VQE step, or DMRG bond dimension, or SHCI epsilon) back to the dispatch row while a job runs. |
| E9 | **HPC/PDB/classification helper scripts** | `avas_probe.py`, `build_qm_cluster.py`, `calibrate_bond_dimension.py`, `dmrgscf_block2.py`, `protonate.py`, `promote_gate1_checks.py`, `merge_jw_hamiltonians.py`, `resubmit_dmrg.py`, `resubmit_hpc.py`, `resubmit_qpu.py`, `solange_screen_and_classify.py`, `experiment_qc_vs_hpc.py` | Python | ✅ Live (varies) | Not individually re-verified in the 2026-09-04 pass; `protonate.py` (wraps `pdb2pqr`) and `promote_gate1_checks.py` are referenced directly by name in the dissertation (§01, §06.iv) and confirmed current there. |
| E10 | **Consistency / benchmark verification** | `scripts/laguna/verify_consistency.py`, `verify_benchmark.py`, `sync_targets.py` | Python | ✅ Live | Run after any `targets.json` change (per CLAUDE.md's standing rule). Both clean (0 contradictions; all benchmarks reproduce) as of 2026-09-04. **Necessary, not sufficient** — checks numeric fields and class letters only, not prose (a real drift in `targets.json`'s STK11 caveat prose vs. the dissertation's own more precise finding was found by a manual audit, not by this checker, 2026-09-04). |
| E11 | **Single source of truth for target facts** | `targets.json` (repo root) | JSON | ✅ Current | `bqp_class`, electron/qubit counts, `structure_caveat`. Frontend's `GENE_MAP`/`MUTATION_CLASS` and the backend's `MUTATION_CONFIGS` overlay are *derived* from this file via `sync_targets.py` — never hand-edited directly. |

## F. Architecture Diagrams & Supporting HTML

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| F1 | **Three-layer architecture** | `public/three_layer_architecture.html` | HTML | ⚠ Not reviewed since 2026-07-03 | |
| F2 | **Platform architecture provenance** | `public/platform_architecture_provenance.html` | HTML | ⚠ Not reviewed since 2026-07-03 | |
| F3 | **System Architecture Report** | `public/System_Architecture_Report.html` | HTML | ⚠ Not reviewed since 2026-07-03 | |
| F4 | **Figures / logos** | `public/figures/*.svg`, `public/logos/cgu-flame.png` | SVG/PNG | ⚠ Not reviewed since 2026-07-03 | |
| — | ~~QC·AI·HPC Architecture~~, ~~P1–P9 Architecture Diagram~~, ~~Unified architecture~~ | — | — | ❌ **Removed** | Previously listed here (`QC_AI_HPC_Architecture.html`, `P1_P9_Architecture_Diagram.html`, `unified_architecture.html`) — confirmed absent from the repo entirely as of 2026-09-04. Either deleted deliberately or lost without a corresponding registry update; not investigated further. |

## G. Build / Project Config

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| G1 | **Vite/React build config** | `package.json`, `tsconfig.json`, `netlify.toml`, `index.html` | Config | ✅ Live | |

## H. Regulatory & Evaluation Documentation

*(New category — none of this existed at the 2026-07-03 snapshot.)*

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| H1 | **Validation protocols** | `docs/regulatory/URS_2026-001.docx`, `IQ_Protocol_v1.0.docx`, `OQ_Protocol_v1.0.docx`, `RTM_2026-001.docx` | Word | ⚠ Not reviewed | User Requirements Spec, Installation/Operational Qualification, Requirements Traceability Matrix — a computer-system-validation-style documentation set. Not cross-checked against the current system in this pass. |
| H2 | **Benchmark protocol** | `docs/benchmark/Benchmark_Protocol_v1.0.docx` | Word | ⚠ Not reviewed | |
| H3 | **Phase 4 DSR evaluation instrument (draft)** | `docs/evaluation/LAYER3B_INSTRUMENT_DRAFT.md` | MD | ✅ Current, explicitly draft | An output-level rubric for DP5, the one design principle no other Phase 4 instrument tests. **Not pre-registered, not committed to** — freeze on OSF only immediately before the actual panel (Jan–Apr 2027), and re-verify every case against the codebase first (per its own §6 checklist) since a case can drift out of sync with the system, exactly like the artifacts in this registry can. |

## I. Project Memory / Registry

| # | Artifact | File(s) | Type | Status | Notes |
|---|---|---|---|---|---|
| I1 | **Project memory** | `CLAUDE.md` | MD | ✅ Current | The single most current source for project facts — check here before trusting a stale note in any other artifact, including this one. |
| I2 | **This registry** | `ARTIFACTS.md` | MD | ✅ Current | Rewritten 2026-09-04 after the prior 2026-07-03 snapshot was found to describe a retired feature (guest login) and omit ~2 months of the most significant architectural work (Gate 1/2, the DMRG+SHCI classifier, LEON, the Rung 3→4 QPU path) entirely. |

## J. External / Not in Repo

| # | Artifact | Location | Type | Notes |
|---|---|---|---|---|
| J1 | **Committee review** | uploaded (not committed) | Word | Prof. Shabtai committee-review doc; drove the original 8 dissertation fixes (predates this session's further work). |
| J2 | **BLAIS 2026 proposal** | `docs/blais_2026/` | PDF ×3 | Submitted 2026-08-15, 19:58 PDT (the deadline day). Proposal and budget both still carry the old, incorrect SOLANGE acronym expansion. The budget's core justification (R175H as the sole Class A hardware-time candidate) was superseded by a warm-started re-run that reclassified R175H to Class B roughly 10 hours *before* submission (commit `02a6326`, 2026-08-15 17:01:24 UTC) — known, discussed, and deliberately left as-is pending any committee response, on the reasoning that a fast-moving project's preliminary evidence is expected to keep evolving. |
| J3 | **Deployed site** | `https://solange-platform.bio` | Live URL | Netlify. **No shared/guest credentials exist any more** — every account is individually provisioned via Supabase Auth (see C1). Do not record a shared login here or anywhere else. |

---

## ⚠ Maintenance Notes

- **Symlink, not duplication, for the dissertation.** `dissertation_revised.html` at repo root is a **symlink** to `public/dissertation_revised.html` — there is exactly one real file. This replaced an earlier manual-sync arrangement after two copies drifted apart once already; never recreate it as a second real file.
- **`theoretical_foundations.html` is still manually synced** (root and `public/` are two real files, confirmed byte-identical 2026-09-04) — unlike the dissertation, this one has no symlink protection, so a future edit to only one side is a real risk. Consider symlinking it the same way if this becomes a recurring issue.
- **`dist/`** is a gitignored build output — stale copies live there; ignore, it rebuilds.
- **arXiv PDF (A2)** needs regeneration from source, and a fresh drift-check against the 2026-09-04 dissertation pass, before submission.
- **This registry itself goes stale fast.** The 2026-07-03 version was ~2 months out of date by the time anyone checked it against the live system. Re-verify categories D, E, and H against the actual repo contents before trusting this file at a glance in the future, rather than assuming it was kept current.
