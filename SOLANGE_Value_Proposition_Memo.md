**MEMORANDUM**

**To:** Prof. Itamar Shabtai
**From:** Doron Cohen
**Re:** SOLANGE™ — Value Proposition Summary (for review prior to forwarding to Molly)
**Program:** IST 697, Claremont Graduate University

---

Itamar,

Below is a summary of SOLANGE's value proposition — who it serves, why it is differentiated from existing platforms, and how its two core layers (governance and translational) fit together. I'd appreciate your review before this goes to Molly. It is not yet incorporated into the dissertation body; it is intended as a standalone positioning document.

*Revision note: Section 6.3 ("Commercialization Potential") has been revised in response to your comment that the answer to "Why would a pharma company spend money on this?" needed to be sharper — it now leads with business outcomes (reduced workflow fragmentation, reduced integration effort, improved reproducibility, audit-ready computational evidence, faster transition from computational results to actionable decisions) rather than with platform features.*

---

## 1. Conceptual Lineage: Why Feynman Matters Here

Richard Feynman (1918–1988), in his 1981/82 paper *"Simulating Physics with Computers,"* argued that classical computers cannot efficiently simulate quantum systems — the computational cost grows exponentially with particle count — but a computer built from quantum components could. This is the founding conceptual argument for quantum computing, and it anchors SOLANGE's literature lineage:

**Feynman (1982) → Lloyd (1996, universal quantum simulation) → Abrams & Lloyd (1999, phase estimation → eigenvalues) → Aspuru-Guzik (2005, applied to real molecules) → modern VQE/sqDRIFT methods used in SOLANGE today.**

## 2. End Users and Their Need

| End User | Why They Need SOLANGE |
|---|---|
| Computational/medicinal chemists (pharma & biotech) | Require ground-state energies and electronic structure for mutation sites that classical methods cannot handle past ~18 active electrons (e.g., TP53 C275F requires ~44e/88q). |
| Academic and translational cancer research labs | Study "non-druggable" tumor suppressors (TP53, STK11, KEAP1) that have no known small-molecule binding pocket, and need a computational foothold. |
| Regulatory/compliance teams within pharma | The P1–P9 provenance schema, aligned with FDA 21 CFR §11.10(e), means SOLANGE output could feed directly into an IND/NDA submission package rather than remaining exploratory research. |
| Quantum hardware partners (e.g., IBM) | Gain a flagship oncology use case demonstrating real-world value for their fault-tolerant hardware roadmap (Heron r3 → Starling → Blue Jay). |

**Core unmet need:** TP53, STK11, and KEAP1 mutations are classified "non-druggable" because no current drug binds them — there is no actionable pocket identifiable via classical structural methods. A large fraction of NSCLC patients carry these mutations and have limited targeted-therapy options today.

## 3. Differentiation: Why SOLANGE, and Is There Anything Comparable?

Adjacent precedents exist, but none combine all three elements SOLANGE does.

| Platform | What It Does | Why It Is Not SOLANGE |
|---|---|---|
| Qiskit Runtime (IBM) / Orquestra (Zapata) | First-generation, single-vendor quantum orchestration | Generic infrastructure — no oncology target, no regulatory layer |
| IonQ + AstraZeneca + AWS + Nvidia (June 2025) | Hybrid pipeline for Suzuki–Miyaura coupling | Generic small-molecule chemistry; no compliance framework |
| Quantinuum + Nvidia (ADAPT-GQE, imipramine) | Cross-vendor hybrid quantum chemistry demonstration | Generic chemistry; no oncology specificity; no audit trail |
| Cleveland Clinic + IBM (March 2026, *Science*) | QCSC workflow, Trp-cage miniprotein, embedded SQD on Heron r2 | Closest architectural precedent, but generic protein folding with zero Part 11 / regulatory provenance |
| Microsoft + Quantinuum (*Nature*, 800× error reduction) | Fault-tolerant QEC + chemistry demonstration | Pure hardware/algorithm demonstration — no disease application, no clinical framing |

**What is actually unique about SOLANGE** is the intersection of three elements that no other platform combines:

1. **Oncology-mutation-specific targeting** — TP53/STK11/KEAP1, not generic chemistry.
2. **Cross-vendor, bifurcation-aware orchestration** — IBM-native *and* Nvidia-substrate, not locked to a single vendor.
3. **Regulatory-grade provenance built into the workflow itself** (P1–P9, FDA 21 CFR Part 11) — the piece absent everywhere else, including the Cleveland Clinic/IBM result.

**Honest framing for the committee:** SOLANGE does not claim to out-compute IBM, Cleveland Clinic, or Microsoft on raw quantum chemistry — they have far greater resources. The contribution is one of **systems and information architecture**, applied to a disease-specific target that no one else addresses. This is why the dissertation is filed as an **Information Systems / DTech** contribution rather than a chemistry or physics contribution.

**Anticipated committee pushback:** *"Couldn't a major player simply add a compliance layer to their existing platform tomorrow?"*
**Response:** This has not happened yet; the gap is real today. The contribution is the formal specification of *how* such a layer should work — the P1–P9 schema itself — not a claim that it is technically impossible for another organization to replicate.

## 4. The "Single Platform" Argument

Today, replicating what SOLANGE does manually would require stitching together six separate systems:

1. A classical chemistry tool (PySCF, Gaussian) for the molecular model and integrals.
2. A quantum SDK (Qiskit, PennyLane, Cirq — and a different one again for Nvidia's CUDA-Q).
3. A separate noise-characterization/error-mitigation toolchain (Phase 3B / sqDRIFT).
4. A separate HPC job scheduler for classical simulation at scale.
5. A separate, typically manual, compliance/documentation system (spreadsheets, lab notebooks) to reconstruct an audit trail after the fact.
6. Manual translation of raw energy/Hamiltonian output into a form a drug-discovery decision-maker can act on.

This argument operates on two levels:

- **Practical/usability level:** one platform, one login, one data format, no manual re-entry, faster iteration.
- **Compliance/scientific-integrity level (the stronger argument for the committee):** every manual hand-off between disconnected tools is a point where the audit trail can break — a point where a regulator could reasonably ask, *"How do we know this number came from this run, on this date, on this hardware?"* with no satisfactory answer.

**Key framing:** the single-platform argument is not separate from the regulatory argument — it is the mechanism that makes the regulatory argument credible. Nothing falls through the cracks between systems, because there are no hand-offs between disconnected tools to create cracks in the first place.

## 5. Governance Layer vs. Translational Layer

- **Governance layer** — tracks, records, and audits every step of the quantum computation: who ran what, when, on what hardware, with which code version (the P1–P9 provenance schema), satisfying FDA 21 CFR §11.10(e) requirements for legally valid, tamper-evident, auditable records.
- **Translational layer** — bridges raw scientific output (ground-state energies, Hamiltonians) to a form usable in drug development, converting physics results into something a drug-discovery decision-maker can interpret and act on, and connecting that output to an actual regulated development pipeline rather than a merely publishable result.

Together, these two layers constitute SOLANGE's real contribution. The quantum chemistry techniques themselves — VQE, CASSCF, Jordan-Wigner encoding — are standard and already in use by major industry players. The contribution is the architecture that makes quantum-derived oncology data **trustworthy and usable** in a regulated context.

---

## 6. Beyond NSCLC: Generalization, Additional Interfaces, and Commercialization

*Addendum — for internal/committee discussion only. Not intended for the dissertation body.*

### 6.1 Is SOLANGE a General Platform, or NSCLC-Specific?

Architecturally, SOLANGE is already pan-cancer, not NSCLC-only. The backend's expansion gene map already carries config entries — Jordan-Wigner Hamiltonian source, PDB structure, active electron/orbital counts — for more than twenty additional tumor-suppressor/loss-of-function genes spanning multiple cancer types (e.g., VHL, BAP1, PBRM1, BRCA1/2, IDH1/2, ATM, RB1, NF1/2, TERT, CDKN2A), all driven through the same schema used for TP53/STK11/KEAP1. NSCLC is the *validated flagship use case* for the dissertation, not a technical ceiling.

**Recommended framing:** present NSCLC as the rigorously validated end-to-end demonstration (richest clinical and structural data, clearest "non-druggable" narrative), while stating explicitly that the underlying architecture generalizes to any loss-of-function mutation with a resolvable structure. This is more defensible to a committee than an unqualified pan-cancer claim, while still capturing the generalizability point reviewers look for in a systems/architecture contribution.

### 6.2 Candidate Additional Interfaces

The highest-value additions extend SOLANGE's actual differentiator — the governance/translational layers — rather than simply adding more disease coverage:

1. **Regulatory export interface** — a structured P1–P9 → IND/NDA submission-package template. This operationalizes the "regulatory-grade provenance" claim rather than leaving it as an assertion.
2. **Clinical/variant-interpretation interface** — extends the existing NGS-report upload pathway from 3D visualization into a ranked druggability/feasibility decision-support summary, reaching the translational-layer audience (clinicians, tumor boards) directly.
3. **Programmatic API/SDK** — allows a pharma data-science team to call SOLANGE from their own pipeline rather than through the web UI. This is typically the interface that actually drives partnership and licensing conversations, more so than the front-end.

Lower-priority additions (e.g., expanded visualization/chart types) would not materially advance either the scientific or commercialization argument and are not recommended as a near-term focus.

### 6.3 Commercialization Potential

*Revised per Prof. Shabtai's review comment: sharpen the answer to "Why would a pharma company spend money on this?" by leading with business outcomes — reduced workflow fragmentation, reduced integration effort, improved reproducibility, audit-ready computational evidence, and faster transition from computational results to actionable drug-discovery decisions — rather than leading with platform features (P1–P9, orchestration, cross-vendor support).*

#### 6.3.1 Why Would a Pharma Company Spend Money on This?

A pharma R&D organization is not paying for VQE itself — that algorithm is freely available. It would pay to remove **cost, delay, and risk** from the path between a computational result and a decision it can act on or defend to a regulator. SOLANGE converts each technical capability into a business outcome:

| Outcome | What it replaces today | Why it has business value |
|---|---|---|
| **Reduced workflow fragmentation** | Manually stitching together a classical chemistry tool, a quantum SDK, a noise-mitigation toolchain, an HPC scheduler, and a separate compliance system for every target | Fewer hand-offs means fewer person-hours per target and fewer places a project stalls waiting on a different team or tool |
| **Reduced integration effort** | Every new hire or partner team re-learning multiple vendor SDKs (Qiskit, PennyLane, CUDA-Q) and re-building glue code between them | Lower onboarding cost and faster time-to-first-result for a computational chemistry team adopting the platform |
| **Improved reproducibility** | Re-running an analysis and getting a different answer because the exact code version, hardware, or parameters used originally were not captured | Reproducibility failures are a recognized, measurable cost driver in pharma R&D (wasted analyst time, delayed go/no-go decisions); a platform that records every run by construction removes that failure mode rather than managing it after the fact |
| **Audit-ready computational evidence** | A separate, often retrospective effort — typically weeks of work per study — to reconstruct an audit trail before a regulatory filing | Removes a known bottleneck and source of risk immediately before a submission deadline, rather than discovering a documentation gap under regulatory review |
| **Faster transition from computational result to actionable decision** | Raw energies/Hamiltonians sitting in a research notebook, requiring manual translation before a chemistry or regulatory team can use them | Shortens the distance between "the computation finished" and "a decision-maker can act on it," compressing the discovery-to-decision cycle |

**The net pitch:** SOLANGE does not ask a pharma buyer to evaluate a new quantum algorithm. It asks them to evaluate whether removing integration overhead, reproducibility risk, and audit-trail risk from their existing computational-chemistry workflow is worth paying for — a much easier internal business case to make than "should we invest in quantum computing."

#### 6.3.2 Commercialization Paths

Realistic paths, in approximate order of near-term plausibility — each is a route to monetizing the outcomes above, not the underlying technical features:

- **B2B SaaS/API licensing** to pharma and biotech computational-chemistry teams — priced and sold on reduced integration effort and audit-readiness, not on access to a VQE implementation; this is the defensible moat, since a generic notebook cannot be sold as a compliance product.
- **Academic/CRO licensing** to translational cancer research labs without in-house quantum computing expertise — sold on removing the need to build and maintain a multi-vendor toolchain in-house.
- **Co-development or data partnership with a quantum hardware vendor** (e.g., IBM) — a flagship oncology workload of this kind is exactly what hardware vendors look to co-fund or showcase; the existing Heron r3 relationship is a natural starting point.
- **Direct-to-clinic regulatory submission product** — plausible in principle but on a five-to-ten-year regulatory-approval horizon; not a near-term commercialization claim.

**Recommended framing for any external-facing material:** lead with the outcomes a buyer already budgets for — reduced integration cost, reduced reproducibility risk, audit-readiness, faster decision cycles — and treat the P1–P9 schema and cross-vendor orchestration as the *mechanism* that delivers them, not as the headline. The NSCLC application remains the proof of concept demonstrating that the mechanism works, not the product itself.

---

Please let me know if anything here needs adjustment before it goes to Molly.

Best,
Doron
