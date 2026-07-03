# The Classical–Quantum Architectural Continuum
### Theoretical Foundations — literature review, TOGAF lens, and research focus
*Dissertation · IST 697 · Doron Cohen · CGU 2026 · Supervisor: Prof. Itamar Shabtai*

> **Note.** All seed citations are marked `[verify]` and must be confirmed against
> the literature database before inclusion.

---

## 1. The Focus — a Venn Synthesis

The contribution of this dissertation is not located in any single mature field,
but precisely in the region where **three** of them intersect:

1. **Classical HPC** — CCSD(T), CASSCF, and full-CI electronic-structure methods
   executed on research-computing infrastructure (USC CARC, NSF ACCESS). This pole
   is bounded by a hard tractability ceiling at roughly **18 active electrons**.
2. **Quantum computing** — variational (VQE) and sample-based quantum diagonalization (sqDRIFT) algorithms
   on near-term hardware (IBM Heron r3), spanning the NISQ-to-fault-tolerant
   transition.
3. **Computational oncology** — the non-druggable NSCLC tumour-suppressor targets
   (TP53, STK11, KEAP1) and the C275F ground-state problem, for which no
   small-molecule binding pocket exists.

**One-sentence focus (paste-ready):**

> *SOLANGE™ is the first provenance-complete, enterprise-architected pipeline that
> transports a non-druggable NSCLC tumour-suppressor mutation (C275F) from the
> classical ground-state ceiling — CCSD(T)/CASSCF at ~18 active electrons — onto
> quantum hardware (VQE today, sqDRIFT sample-based quantum diagonalization on IBM Heron r3 tomorrow) under a
> 21 CFR Part 11 audit trail.*

*(Figure 1 — the three-set Venn diagram — is rendered in the companion
`theoretical_foundations.html`.)*

---

## 2. An Enterprise-Architecture Lens — TOGAF

TOGAF (The Open Group Architecture Framework, Standard 10th Edition) is the most
widely adopted enterprise-architecture methodology. It contributes two things to
this dissertation. First, a **recognised vocabulary** that decomposes a complex
socio-technical system into discrete architecture domains. Second, the
**Architecture Development Method (ADM)** — a cyclical process that models
architecture not as a static artefact but as something that *evolves* through
governed iterations. This is the precise conceptual instrument needed to describe
SOLANGE's Phase 3A → 3B trajectory as a principled architectural migration rather
than an ad-hoc engineering progression.

### 2.1 SOLANGE mapped onto TOGAF's four domains

| TOGAF domain | Framework definition | SOLANGE™ instantiation | Provenance |
|---|---|---|---|
| **Business** | Strategy, drivers, stakeholders, the problem solved | Non-druggable NSCLC tumour-suppressor targets (TP53, STK11, KEAP1); the unmet need to characterise mutations with no druggable pocket | — |
| **Data** | Logical/physical data assets and their governance | Jordan–Wigner-encoded molecular Hamiltonians, CASSCF reference datasets, immutable `P1–P9` provenance records (Supabase, 21 CFR §11.10(e)) | P1 · P2 · P5 |
| **Application** | Applications, their behaviour and interactions | CASSCF/PySCF (classical reference), VQE on PennyLane (live Phase 3A), sqDRIFT sample-based quantum diagonalisation (Phase 3B) | P1 · P6 |
| **Technology** | Hardware, networks, platform infrastructure | The classical→quantum bridge: HPC (USC CARC / NSF ACCESS, Render cloud) ↔ IBM Heron r3 QPU | P3 · P4 |

### 2.2 The ADM cycle as the Phase 3A → 3B narrative

A single turn of the ADM cycle maps directly onto SOLANGE's maturation:

```
Phase 3A (minimal)        Phase 3A (intermediate)                    Phase 3B
4-qubit CAS(2e,2o) VQE  →  scientific ceiling: ~24e/48q (~4.5 PB)  →  94+ qubit full site
live on classical cloud    operational run:    ~20e/40q (~17.6 TB)     IBM Heron r3, sqDRIFT
[architecture validated]   VQE *simulated* on HPC (CARC/ACCESS)        [true quantum regime]
                           [2 orders of magnitude above minimal]
```

A clarification that matters scientifically: the intermediate step is **not**
classical chemistry — that is impossible past ~18e for CCSD(T). It is a *quantum*
algorithm (VQE) executed on a classical HPC **simulator** of the quantum circuit.
The scientific ceiling target is ~24e/48q — precisely where both classical
chemistry and classical state-vector simulation (~4.5 PB) break down. The
operational run uses ~20e/40q (~17.6 TB, within NSF ACCESS / USC CARC reach),
which fully demonstrates orchestration-pipeline scaling at meaningful
quantum-circuit width without approaching the memory limit.

Throughout this migration the **provenance layer (P1–P9)** plays the role of
TOGAF's *Architecture Governance* — the continuity-of-record discipline that makes
each iteration auditable.

### 2.3 Basis-set scope and Phase 3A target selection

The STO-3G minimal basis set was selected to maintain a 4-qubit active space
compatible with current quantum hardware, with the understanding that chemically
accurate results require larger basis sets in Phase 3B.

Phase 3A results are reported exclusively for the eight validated per-site
targets — TP53 C275F, TP53 Y220C, KEAP1 LOF, KEAP1 G333C, KEAP1 R320Q,
STK11 LKB1, STK11 F354L, and STK11 D194N — each associated with a distinct
model compound chosen to represent the key sidechain interaction at the mutation
site (toluene for Phe, methanethiol for Cys, acetic acid for Asp, acetamide for
Gln, isobutane for Leu). All remaining patient-report variants are designated
Phase 3B and flagged for IBM Heron r3 computation, as no per-mutation Phase 3A
Hamiltonian exists for them. Presenting residue-class surrogates as per-mutation
results would constitute scientific misrepresentation and is explicitly excluded
from the SOLANGE™ pipeline by design.

---

## 3. Literature-Review Skeleton

Four movements — physical theory → hybrid algorithms → systems/infrastructure
evolution → enterprise-architecture lens — closing on the research gap.

### 3.1 Physical foundations of quantum simulation of chemistry
Establishes *why* a quantum computer is the natural instrument for ground-state
energetics, and where classical methods break down.
- Feynman — *Simulating Physics with Computers*. `[verify: Feynman 1982, Int. J. Theor. Phys.]`
- Lloyd — *Universal Quantum Simulators*. `[verify: Lloyd 1996, Science]`
- Classical tractability ceiling — CCSD(T), CASSCF, full-CI scaling; the ~18-electron wall. `[verify: Helgaker, Jørgensen & Olsen, Molecular Electronic-Structure Theory]`
- Comprehensive survey. `[verify: Cao et al. 2019, Chem. Rev.]`

### 3.2 The hybrid classical–quantum algorithmic layer
The methods that straddle both paradigms — the algorithmic content of the Venn
intersection.
- VQE. `[verify: Peruzzo et al. 2014, Nat. Commun.]`
- Theory of variational hybrid quantum-classical algorithms. `[verify: McClean et al. 2016, New J. Phys.]`
- NISQ-era framing and limits. `[verify: Preskill 2018, Quantum]`
- qDRIFT — stochastic compiler for Hamiltonian simulation (basis of sqDRIFT). `[verify: Campbell 2019, Phys. Rev. Lett.]`
- VQE review — methods and best practices. `[verify: Tilly et al. 2022, Phys. Rep.]`

### 3.3 Architectural evolution — from HPC to quantum-accelerated pipelines
The systems/infrastructure lineage: the classical HPC pole and the migration toward
fault tolerance.
- Quantum algorithms for chemistry/materials — resource estimates. `[verify: Bauer et al. 2020, Chem. Rev.]`
- Academic HPC as the classical pole — USC CARC and the NSF ACCESS allocation model
  as exemplars of the research-computing substrate. `[verify: USC CARC; NSF ACCESS]`
- Quantum-HPC integration / hybrid orchestration patterns. `[verify: recent quantum-HPC integration literature]`
- Fault-tolerance roadmaps and the transitional window (Heron → Starling → Blue Jay). `[verify: IBM Quantum roadmap, primary sources]`

### 3.4 Enterprise architecture as a design discipline
The lens itself, plus the provenance/compliance architecture that governs the whole.
- TOGAF Standard, 10th Edition — ADM and the four domains. `[verify: The Open Group 2022]`
- EA modelling for complex systems. `[verify: Lankhorst, Enterprise Architecture at Work; ArchiMate]`
- Provenance & reproducibility — FAIR data principles. `[verify: Wilkinson et al. 2016, Sci. Data]`
- Regulatory-grade audit trails. `[verify: FDA 21 CFR Part 11]`

### 3.5 Synthesis & research gap
The literature treats each circle of Figure 1 maturely **in isolation**. What is
absent — the gap this dissertation fills — is a single architecture that:
- **(a)** is framed by an established enterprise-architecture method (TOGAF/ADM);
- **(b)** carries a problem end-to-end across the classical→quantum technology boundary;
- **(c)** maintains a regulatory-grade provenance chain (P1–P9); and
- **(d)** targets specifically the *non-druggable* tumour-suppressor class for which
  no pharmacological pocket exists.

SOLANGE™ is that intersection.

---

*Companion: `theoretical_foundations.html` (rendered Venn diagram + styled tables).*
