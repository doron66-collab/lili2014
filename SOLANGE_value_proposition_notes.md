# SOLANGE™ — Value Proposition Notes
*Working notes for committee/advisor discussion (Molly). Not part of the dissertation — for talking-points and positioning only.*

---

## 1. Who is Feynman, and why does he matter to SOLANGE?

**Richard Feynman** (1918–1988) — Nobel laureate physicist (QED, 1965), inventor of Feynman diagrams.

In 1981/82 he argued, in *"Simulating Physics with Computers,"* that classical computers cannot efficiently simulate quantum systems (the cost explodes exponentially with particle count), but a computer *built from quantum components* could. This is the founding conceptual argument for quantum computing.

SOLANGE's literature lineage traces directly back to this insight:
**Feynman (1982) → Lloyd (1996, universal quantum simulation) → Abrams & Lloyd (1999, phase estimation → eigenvalues) → Aspuru-Guzik (2005, applied to real molecules) → modern VQE/SqDRIFT methods used in SOLANGE today.**

---

## 2. Who are the end users, and why would they want SOLANGE?

| End user | Why they'd want it |
|---|---|
| **Computational/medicinal chemists (pharma & biotech)** | Need ground-state energies / electronic structure for mutation sites classical methods can't handle past ~18 active electrons (e.g., TP53 C275F needs ~44e/88q). |
| **Academic/translational cancer research labs** | Studying "non-druggable" tumor suppressors (TP53, STK11, KEAP1) with no known small-molecule binding pocket — need a computational foothold. |
| **Regulatory/compliance teams inside pharma** | P1–P9 provenance + FDA 21 CFR Part 11 alignment means SOLANGE output could feed directly into an IND/NDA submission package, not just exploratory research. |
| **Quantum hardware partners (IBM, etc.)** | A flagship oncology use case demonstrating real-world value of their fault-tolerant roadmap (Heron r3 → Starling → Blue Jay). |

**Core unmet need:** TP53/STK11/KEAP1 mutations are "non-druggable" because no current drug binds them — there's no actionable pocket via classical structural methods. A large fraction of NSCLC patients have these mutations and limited targeted-therapy options today.

---

## 3. Why SOLANGE vs. the big players? Is there anything like it?

Not really — adjacent precedents exist, but none combine all three things SOLANGE does.

| Platform | What it does | Why it's not SOLANGE |
|---|---|---|
| Qiskit Runtime (IBM) / Orquestra (Zapata) | First-gen single-vendor quantum orchestration | Generic infrastructure — no oncology target, no regulatory layer |
| IonQ + AstraZeneca + AWS + Nvidia (June 2025) | Hybrid pipeline for Suzuki–Miyaura coupling | Generic small-molecule chemistry; no compliance framework |
| Quantinuum + Nvidia (ADAPT-GQE, imipramine) | Cross-vendor hybrid quantum chemistry demo | Generic chemistry; no oncology specificity; no audit trail |
| **Cleveland Clinic + IBM** (March 2026, *Science*) | QCSC workflow, Trp-cage miniprotein, embedded SQD on Heron r2 | **Closest precedent architecturally** — but generic protein folding, **zero Part 11 / regulatory provenance** |
| Microsoft + Quantinuum (*Nature*, 800× error reduction) | Fault-tolerant QEC + chemistry demo | Pure hardware/algorithm demo — no disease application, no clinical framing |

**What's actually unique about SOLANGE** is the intersection of three things nobody else combines:
1. **Oncology-mutation-specific targeting** (TP53/STK11/KEAP1, not generic chemistry)
2. **Cross-vendor, bifurcation-aware orchestration** (IBM-native *and* Nvidia-substrate, not locked to one vendor)
3. **Regulatory-grade provenance built into the workflow itself** (P1–P9, FDA 21 CFR Part 11) — the piece that is *absent everywhere else*, including the Cleveland Clinic/IBM result.

**Honest framing for the committee:** SOLANGE isn't claiming to out-compute IBM/Cleveland Clinic/Microsoft on raw quantum chemistry — they have far more resources. The contribution is **systems/information-architecture**, applied to a disease-specific target nobody else addresses. This is why the dissertation is filed as an **Information Systems / DTech** contribution, not chemistry or physics.

**Anticipated pushback to prepare for:** *"Couldn't a big player just bolt a compliance layer onto their existing platform tomorrow?"*
→ Answer: this hasn't happened yet, the gap is real today, and the contribution is formally specifying *how* that layer should work (the P1–P9 schema itself) — not a claim that it's technically impossible for someone else to do it.

---

## 4. The "one platform, no jumping between systems" argument

Today, doing what SOLANGE does manually would require stitching together:
1. A classical chemistry tool (PySCF, Gaussian) for the molecular model/integrals
2. A quantum SDK (Qiskit, PennyLane, Cirq — different one again for Nvidia's CUDA-Q)
3. A separate noise-characterization/error-mitigation toolchain (Phase 3B / SqDRIFT)
4. A separate HPC job scheduler for classical-simulation-at-scale
5. A separate, usually manual, compliance/documentation system (spreadsheets, lab notebooks) to reconstruct an audit trail after the fact
6. Manual translation of raw energy/Hamiltonian outputs into something a drug-discovery decision-maker can act on

**Two layers to this argument — use both:**
- **Practical/usability layer:** one platform, one login, one data format, no manual re-entry, faster iteration.
- **Compliance/scientific-integrity layer (the stronger one for the committee):** every manual hand-off between disconnected tools is a place the audit trail can break — a place a regulator could ask "how do we know this number came from this run, on this date, on this hardware?" with no good answer.

**Key framing:** *The single-platform argument is not separate from the regulatory argument — it's the mechanism that makes the regulatory argument credible.* Nothing falls through the cracks between systems because there are no cracks — no hand-offs between disconnected tools.

---

## 5. Governance layer vs. translational layer (terminology clarification)

- **Governance layer** — tracks, records, and audits every step of the quantum computation: who ran what, when, on what hardware, with what code version (P1–P9 provenance), satisfying FDA 21 CFR §11.10(e) — legally valid, tamper-evident, auditable records.
- **Translational layer** — bridges *raw scientific output* (ground-state energies, Hamiltonians) to *something usable in drug development* — converting physics results into a form a drug-discovery decision-maker can interpret and act on, and connecting that into an actual regulated development pipeline (not just a publishable result).

Together, these two layers are SOLANGE's real contribution — not the quantum chemistry techniques themselves (VQE, CASSCF, Jordan-Wigner are standard, already used by the big players), but the architecture that makes quantum-derived oncology data **trustworthy and usable** in a regulated context.

---

*Status: not yet incorporated into the dissertation body — kept separate per request, for talking-points use ahead of meeting with Molly.*
