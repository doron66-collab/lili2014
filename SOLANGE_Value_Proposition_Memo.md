**MEMORANDUM**

**To:** Prof. Itamar Shabtai
**From:** Doron Cohen
**Re:** SOLANGE™ — Value Proposition Summary (for review prior to forwarding to Molly)
**Program:** IST 697, Claremont Graduate University

---

Itamar,

Below is a summary of SOLANGE's value proposition — who it serves, why it is differentiated from existing platforms, and how its two core layers (governance and translational) fit together. I'd appreciate your review before this goes to Molly. It is not yet incorporated into the dissertation body; it is intended as a standalone positioning document.

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

Please let me know if anything here needs adjustment before it goes to Molly.

Best,
Doron
