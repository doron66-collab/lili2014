# SOLANGE Copilot: A Privacy-Preserving Local-LLM Assistant for Auditable Quantum-Oncology Provenance

**Author:** Doron Cohen
**Instructor:** Prof. Wallace Chipidza
**Course:** IST 362 — Emerging Technologies
**Institution:** Claremont Graduate University, 2026
**Repository:** https://github.com/doron66-collab/ist362-solange-copilot (run instructions in Appendix A)

> *Evaluation numbers are from real runs on the author's machine (Apple Silicon, Ollama); tables are reproducible with the harness in Appendix A.*

---

## Abstract

Quantum-oncology platforms such as SOLANGE produce results that are simultaneously **complex** (variational quantum eigensolver energies, provenance records with cryptographic seals) and **stochastic** (quantum sampling, variational noise). Explaining these results to reviewers, and doing so over sensitive genomic and clinical data, is a bottleneck. This project asks a focused question in emerging technologies: *are small language models now good enough to do useful scientific document work entirely locally, so that no data leaves the machine?* We built the **SOLANGE Copilot**, a retrieval-augmented local-LLM application that (A) turns a SOLANGE P1–P9 provenance record into an auditable plain-English narrative and (B) explains why a tumor-suppressor mutation is non-druggable, grounded in real peer-reviewed literature with citations. The copilot runs four local models via **Ollama** and never contacts a cloud API. We evaluate the models on extraction accuracy against a ground-truth record, grounding accuracy, hallucination-refusal, per-mutation fact coverage, and a repeated-run sampling-temperature (randomness) study. A 3-billion-parameter model (`llama3.2`) matched an 8-billion reasoning model on every quality metric while running ~12× faster, and enriching the retrieval corpus lifted the smallest (0.6 B) model to perfect grounding. Our organizing thesis: **reduce non-determinism where it is reducible (temperature = 0 at the language-model layer) and notarize it where it is irreducible (cryptographic seals over stochastic quantum sampling)** — together yielding a reproducible, auditable pipeline suited to regulated electronic-record settings.

---

## 1. Introduction

### 1.1 Context

SOLANGE is a 3D quantum-simulation platform for **non-druggable NSCLC tumor-suppressor mutations** (TP53, STK11, KEAP1). Its Phase 3A computes molecular ground-state energies for mutation proxies with the Variational Quantum Eigensolver (VQE) on a Jordan-Wigner Hamiltonian built from a PySCF CASSCF active space; Phase 3B targets real quantum hardware (IBM Heron r3) for problems beyond the classical limit. Every run emits a nine-part provenance record (P1–P9) whose P8 field is a SHA-256 seal, recomputed and checked by a notarization component (LEON) under a "verify, don't trust" principle aligned with FDA 21 CFR §11.10(e).

### 1.2 Problem

Two frictions motivate this project. First, provenance records and quantum results are **hard to read**: a reviewer must parse nested JSON, energies in Hartree, active-space notation, and seals. Second, the data is **sensitive**; genomic and clinical context should not be sent to a third-party cloud API. A tool that explains results must therefore be both *trustworthy* (grounded, auditable, non-hallucinating) and *private*.

### 1.3 The emerging-technology question

The emerging capability we investigate is that **small language models are now good enough to run locally** — on a laptop, with no data egress — and still perform extraction, summarization, and grounded question-answering usefully. Local execution converts the classic privacy/compliance trade-off: a capability a cloud API cannot offer. This directly engages the course themes of the benefits, downsides, and risks of AI and quantum computing (the risk we confront head-on is **hallucination**, which we measure and mitigate), and of **choosing among LLMs** — including local models — by use case, cost, speed, and reliability.

Two forces make this timely. First, capable open-weight models now run on commodity hardware through runtimes such as Ollama, so "local" no longer means "toy": the models we test range from 0.6 to 8 billion parameters and all run on a single laptop. Second, the *marginal cost* of a local query is essentially zero — there is no per-token API billing and no monthly minimum — and there is no cloud data-processing agreement to negotiate, because nothing leaves the machine. That combination matters most exactly where budgets and data-governance capacity are thinnest: a single researcher, a small lab, or a clinic handling regulated genomic data. The open question this project interrogates is whether the smaller models that local execution makes practical are *accurate and disciplined enough* for scientific use — a question we answer with measurement rather than assertion.

### 1.4 Thesis

A quantum→AI pipeline is a chain of stochastic components. Our organizing claim is that auditability comes from treating each layer correctly:

> **Determinism where reducible, notarization where irreducible.** Where randomness is optional — the language-model sampling layer — we remove it (temperature = 0), yielding bit-for-bit reproducible explanations. Where randomness is physical — quantum sampling — we cannot remove it, so we *notarize* a specific result with a cryptographic seal. Factual correctness is handled by a third, orthogonal mechanism: retrieval grounding.

### 1.5 Contributions

1. A working, dependency-free local-LLM application with two grounded modes plus free-form Q&A, running on real SOLANGE provenance records at three scales (4, 24, and 88 qubits).
2. A reproducible evaluation harness measuring extraction accuracy, grounding accuracy, hallucination-refusal, per-mutation fact coverage, and a repeated-run temperature study.
3. Empirical findings: a 3 B model is the cost/quality sweet spot; corpus enrichment lifts a 0.6 B model to perfect grounding; sampling temperature increases variability without improving accuracy, so temperature = 0 is the right default for auditable use.

---

## 2. Choice of Model(s) [Section 1]

We use **Ollama** to serve models locally over its HTTP API (`localhost:11434`); no SDK and no network egress are required. We deliberately compare four models spanning three size tiers and two families (general vs. reasoning), because *model choice is itself a finding*:

| Model | Params | Family | Role in the study |
|---|---|---|---|
| `qwen3:0.6b` | 0.6 B | general (thinking) | tiny/fast baseline |
| `deepseek-r1:1.5b` | 1.5 B | reasoning | small reasoning |
| `llama3.2:latest` | 3 B | general | mid, general-purpose |
| `deepseek-r1:8b` | 8 B | reasoning | largest / strongest |

For retrieval we use a dependency-free **BM25** ranker (no embedding model required), with an optional `nomic-embed-text` semantic path. Rationale: a resource-constrained local deployment benefits from a retriever that needs no model and works instantly, while still allowing a semantic upgrade. All generation runs at temperature 0.2 by default and temperature 0 for auditable production use (Section 5).

---

## 3. How the Application Works [Section 2]

### 3.1 Architecture

```
ist362_copilot/
├── backend/
│   ├── serve.py          zero-dependency stdlib HTTP server (recommended)
│   ├── app.py            equivalent FastAPI server (optional)
│   ├── llm_adapter.py    Ollama backend + deterministic mock; auto-detects
│   ├── retriever.py      BM25 retrieval (pure stdlib) + optional embeddings
│   ├── prompts.py        system prompt + per-mode templates (grounding rules)
│   └── data/
│       ├── corpus.json       19-passage literature corpus (real PubMed refs)
│       ├── sample_run.json   evaluation answer key (P1–P9)
│       └── sample_runs.json  real runs at 4 / 24 / 88 qubits
├── frontend/index.html   single-page UI (no build step)
└── eval/                 evaluation harness + datasets
```

A design goal was **zero external dependencies**: the adapter and retriever use only the Python standard library, and `serve.py` is built on `http.server`, so the whole application — including the evaluation — runs on any machine with Python and no `pip install`, and falls back to a deterministic mock backend when no model is installed. This makes the demo reproducible and CI-friendly, and it side-steps package-compilation failures on bleeding-edge Python versions.

### 3.2 Backends: Ollama and mock

`llm_adapter.py` defines one interface with two implementations. `OllamaBackend` POSTs to the local Ollama server and strips reasoning models' `<think>…</think>` traces from the final answer (keeping them available separately). `MockBackend` returns deterministic, extractive text so the app and evaluation run with no model installed. `get_backend()` auto-detects a live Ollama server, else uses the mock. Every API response reports which backend and model produced it — honest provenance for the tool's own outputs.

### 3.3 Retrieval and grounding

For Modes B and the free-form chat, the user's query retrieves the top-*k* passages from a 19-document corpus (15 summarized from real PubMed-indexed articles, each carrying its DOI; the remainder are SOLANGE method facts). The system prompt enforces a strict discipline: **answer only from the provided context, cite bracketed sources, and say "not in the provided context" rather than guess.** This is what makes answers auditable and is exactly what the evaluation measures. The full retrieved passages are returned to the UI so a reader can verify every claim against its source.

### 3.4 The two modes plus chat

- **Mode A — Explain a run.** Input: a SOLANGE P1–P9 provenance record (JSON). The model produces a plain-English narrative — what mutation, what method and settings, the resulting ground-state energy (with units), and how the record is notarized — and checks internal consistency. The record is the sole context; the model must quote exact values. The UI ships three real example records at increasing scale (Section 3.5).
- **Mode B — Druggability.** Input: a mutation (e.g., `TP53 C275F`). The app retrieves literature and the model explains *why* the alteration is non-druggable and what indirect strategies exist, citing sources.
- **Chat.** Free-form grounded Q&A over the corpus with citations.

### 3.5 Multi-scale real records

To demonstrate the copilot beyond a toy 4-qubit proxy, Mode A includes three real SOLANGE records:

| Scale | Active space | Qubits | Method | Backend | Energy |
|---|---|---|---|---|---|
| Small | CASSCF(2,2) | 4 | VQE, 27-term JW Hamiltonian | laptop simulator | −266.49074 Ha |
| Mid | CASSCF(12,12) | **24** | exact statevector diagonalisation | NVIDIA L40S GPU (HPC/Laguna) | −437.76124 Ha |
| Large | CASSCF(44,44) | 88 | sqDRIFT (sample-based quantum diagonalisation) | IBM Heron r3 | *planned* |

The 24-qubit record is a genuine platform output included verbatim with its real seal; the 4-qubit record carries the real Phase 3A CASSCF values and Hamiltonian re-expressed in the current P1-P9 schema and re-sealed with the platform's own LEON code (so both verify live); the 88-qubit Phase 3B record is a clearly-marked plan (its final energy and seal are intentionally null, since that run has not yet been executed). A subtle but important point the copilot handles correctly: **CASSCF(12,12) maps to 24 qubits**, not 12 — twelve spatial orbitals become 24 spin-orbitals under the Jordan-Wigner transform. Because the copilot operates on the provenance record rather than the physics, it is *scale-agnostic*: the same tool explains a 4-qubit and an 88-qubit run.

### 3.6 Determinism

Language-model sampling temperature is the one knob that injects avoidable randomness. Setting temperature = 0 makes generation deterministic — identical input yields identical output — which is a precondition for reproducible, auditable records. Section 5 measures this directly.

### 3.7 A request, end to end

To make the mechanics concrete, consider a Mode B query for `TP53 C275F`. Five steps run entirely against `localhost`:

1. **Retrieve.** The mutation string, expanded with the terms "undruggable non-druggable therapy," is tokenized and scored against all 19 corpus passages by BM25; the top four are selected. BM25 rewards passages that contain the query's rarer terms and normalizes for passage length, so a short, on-topic passage outranks a long, loosely related one.
2. **Assemble context.** The retrieved passages are rendered as a numbered block — each with its title, full text, and source citation (author, venue, year, DOI) — so the model sees exactly what the reader will later see.
3. **Build the prompt.** A fixed system prompt (the grounding contract, §3.8) is combined with the numbered context and a task instruction ("explain why this alteration is non-druggable and summarize indirect strategies; cite [n]").
4. **Generate locally.** The assembled prompt is POSTed to the local Ollama server; reasoning-model "think" traces are stripped from the final answer and kept separately.
5. **Return with provenance.** The response carries the answer, the exact sources used (so a reader can verify each claim), and the backend, model, and latency that produced it.

The same five-step pipeline serves Mode A — with the provenance record as the sole context instead of retrieved passages — and the free-form chat. Because retrieval and generation are separate stages, the retriever can be swapped (BM25 ↔ embeddings) and the model can be swapped (any Ollama model) without touching the rest of the pipeline.

### 3.8 The grounding contract

The system prompt encodes four rules that convert a general chatbot into an auditable assistant:

1. **Use only the provided context.** The model may not draw on its pretrained knowledge. This both improves faithfulness and makes behavior *consistent across models* — a 0.6 B and an 8 B model are asked to do the same bounded task, which is what makes the head-to-head comparison fair.
2. **Cite bracketed sources inline.** Every claim is traceable to a passage the reader can see and check.
3. **Refuse rather than guess.** If the answer is absent, respond "that is not in the provided context" — the behavior our trap questions measure directly.
4. **Quote exact values with units.** Critical for provenance records, where an energy in Hartree or a qubit count must be reproduced precisely, not paraphrased.

For Mode A we add interpretation rules (the §3.6 consistency guide) so the model explains what values *mean* and checks them against one another rather than echoing them. This contract is the single most important design element in the system: it is what makes a sub-billion-parameter model usable for grounded scientific work, as the evaluation confirms (Section 5.3, where retrieval quality — not model size — moved the smallest model to perfect grounding).

---

## 4. Potential and Demonstrated Capabilities [Section 3]

**Demonstrated.** Running against local models, the copilot: (1) explains the real 24-qubit HPC record correctly — identifying the mutation, exact diagonalisation method, 24-qubit count, ground-state energy of −437.76 Ha, the NVIDIA L40S backend, and the SHA-256/LEON notarization — from the record alone; (2) explains why TP53 C275F is non-druggable, correctly attributing this to its nature as a transcription factor lacking a compact catalytic site, and lists indirect strategies (MDM2 inhibition, restoring mutant p53) with citations; (3) refuses out-of-context questions instead of fabricating answers. Notably, capability (2) was produced by `qwen3:0.6b`, a 522 MB model, grounded and cited.

**Live seal verification (LEON, locally).** Mode A includes a **Verify seal** action that reproduces the platform's exact P8 sealing algorithm — SHA-256 over the canonical serialization of the P1-P7 and P9 fields (excluding the calibration-epoch timestamp), identical to the platform's own `build_p8_seal` in `backend/routes/leon.py` and confirmed bit-for-bit against a real HPC record — and recomputes it on demand. For records the platform itself cannot re-attest (an older compact-format export, or a planned run whose seal is null) the copilot returns the same honest verdicts the production notary does, rather than a false pass or fail. A valid record verifies against its stored hash; editing any sealed value (e.g. an energy) in the record box makes the recomputed hash diverge, so the tool reports a mismatch and, when the sealed snapshot is present, names the altered field. This runs entirely offline and turns "verify, don't trust" from a slogan into a working, local re-attestation.

**Internal-consistency checking.** Beyond summarizing, Mode A *audits* the record. On the 24-qubit run, `llama3.2` correctly reasoned that "this method does not involve a quantum circuit, so it is expected to have a gate count of 0 and a circuit depth of 0" — interpreting the zero values rather than flagging them — and independently verified the qubit count against the active space ("CASSCF(12,12) uses 12 orbitals … a qubit count of 24, which matches the provided value") under the Jordan-Wigner rule, before returning a consistency verdict. This turns the copilot from a summarizer into a lightweight auditor of provenance records — the practical embodiment of "verify, don't trust." (The 3 B model's phrasing of the qubit check was slightly clumsy — it framed a confirmation as a "potential inconsistency" — a minor fluency limitation, though the logic and conclusion were correct.)

**Potential.** Because the tool is record-driven and private, it generalizes to any lab that must explain and audit computational results over sensitive data: batch-narrating provenance records for review, onboarding, or regulatory dossiers; live literature-grounded Q&A during analysis; and — with the seal payload present in real records — automated seal re-verification (a natural next step, Section 7).

---

## 5. Evaluation [Section 4]

### 5.1 Methodology

Every task has a built-in ground truth, so scoring is objective and reproducible (harness in Appendix A). We report four instruments:

1. **Extraction accuracy (Mode A):** does the model read exact fields (basis set, qubit count, energy, notary, algorithm) correctly from the provenance record? Ground truth = the record itself (7 fields).
2. **Grounding accuracy (Mode B / chat):** does it surface the correct fact from retrieved literature (5 questions)?
3. **Hallucination refusal:** for questions whose answer is *not* in the corpus, does it correctly refuse (3 trap questions)?
4. **Latency:** wall-clock seconds per query.

Scoring is keyword-set containment (an answer is correct if it contains any accepted form of the expected value), which is conservative and deterministic.

### 5.2 Main benchmark (enriched 19-document corpus)

| Model | Extraction | Grounding | Hallucination refusal | Mean latency (s) |
|---|---|---|---|---|
| **llama3.2:latest (3 B)** | **1.00** | **1.00** | **1.00** | **1.08** |
| qwen3:0.6b | 0.86 | 1.00 | 1.00 | 1.62 |
| deepseek-r1:1.5b | 0.71 | 0.80 | 1.00 | 4.77 |
| deepseek-r1:8b | 1.00 | 1.00 | 1.00 | 13.62 |

**Findings.** (i) `llama3.2` (3 B) is the **cost/quality sweet spot**: perfect on all three quality metrics and the *fastest* model, beating the 8 B reasoning model on quality (tie) at ~12.6× lower latency. Bigger is not better for this workload. (ii) `deepseek-r1:1.5b` is the weakest (0.71 extraction, 0.80 grounding) *and* slow — the small reasoning model is a poor fit here. (iii) **All four models refuse all trap questions (1.00)** — the grounding discipline holds across every size. (iv) A reasoning model can be *slower than a smaller general model* because it emits hidden "thinking" tokens; size alone does not predict latency.

### 5.3 Effect of corpus enrichment

Enriching the corpus from 10 to 19 real passages produced a clean, controlled result:

| Model | Grounding (10 docs) | Grounding (19 docs) | Extraction (unchanged) |
|---|---|---|---|
| llama3.2 | 1.00 | 1.00 | 1.00 |
| qwen3:0.6b | 0.80 | **1.00** | 0.86 |
| deepseek-r1:1.5b | 1.00 | 0.80 | 0.71 |
| deepseek-r1:8b | 1.00 | 1.00 | 1.00 |

Extraction accuracy is **identical** before and after (it depends on the record, not the corpus) — a control confirming the experiment is well-behaved. Grounding for the smallest model, `qwen3:0.6b`, rose to **1.00**: better retrieval context compensated for small model size, the central promise of retrieval-augmented generation. The `deepseek-r1:1.5b` grounding movement (1.00→0.80) is one question on a five-item set, within the single-question noise band established in Section 5.5, and should not be read as a trend.

### 5.4 Per-mutation drill-down (Mode B)

We also score how completely each model explains a *specific* mutation (fraction of expected fact-groups covered, and whether it cites sources). Representative results (initial 10-document corpus):

| Mutation | llama3.2 | qwen3:0.6b | deepseek-r1:1.5b | deepseek-r1:8b |
|---|---|---|---|---|
| TP53 C275F | 0.75 (cited) | 0.50 (cited) | 0.50 (cited) | 0.75 (cited) |
| STK11 loss | 1.00 (cited) | 1.00 (cited) | 0.50 (cited) | 1.00 (cited) |
| KEAP1 loss | 1.00 (cited) | 1.00 (cited) | 1.00 (**not** cited) | 1.00 (cited) |

Fact coverage is a *completeness* score over 3–4 expected topics, not a correctness score; because the denominator is small, each miss moves it a lot. The one genuinely informative cell is `deepseek-r1:1.5b` on KEAP1, the single case in the entire study where a model failed to cite its sources — concrete evidence that the small reasoning model is the weak link. `llama3.2` is strong across both the benchmark and the drill-down.

### 5.5 Randomness study: effect of sampling temperature [addresses Topic 1]

Sampling temperature is the randomness knob of a language model. Because output is itself random, each temperature was run **5 times**; we report mean ± standard deviation (grounding refusal was 1.00 ± 0.00 for both models at all temperatures and is omitted). Initial corpus:

| Temp | llama3.2 grounding | qwen3:0.6b grounding |
|---|---|---|
| 0.0 | 1.00 ± 0.00 | 0.80 ± 0.00 |
| 0.7 | 0.88 ± 0.11 | 0.76 ± 0.17 |
| 1.3 | 1.00 ± 0.00 | 0.80 ± 0.14 |

**Findings.** (i) At temperature 0 the standard deviation is exactly **0.00** for both models — generation is deterministic and therefore reproducible; this is the empirical basis of our thesis. (ii) Temperature does not change *mean* accuracy significantly (all differences are within one standard deviation), but it **increases variability** (sd rises from 0 to ≈0.15). Randomness here buys unpredictability, not accuracy. (iii) The stronger model (`llama3.2`) is essentially immune; the weaker model (`qwen3:0.6b`) is more variable — **randomness harms consistency more as the model gets smaller.** (iv) Refusal is rock-solid at every temperature. For an auditable scientific tool the implication is direct: **temperature = 0 is the correct default**, giving full reproducibility with no loss of grounding and no increase in hallucination.

### 5.6 Threats to validity

The question sets are small (5 grounding, 3 trap, 3–4 facts per mutation), so single-item flips can masquerade as trends; we mitigate by repeated runs in Section 5.5 and by treating ≤ one-question differences as noise. Keyword-containment scoring can, in principle, credit a right keyword in a wrong sentence; manual spot-checks did not surface such cases. All numbers are from one hardware/software configuration.

### 5.7 Cost and accessibility

Accuracy is only one axis of model choice; cost and access are others, and they are where the local design is strongest. Every query in this study cost nothing beyond electricity: there is no per-token API charge, no monthly minimum, and no cloud data-processing agreement to negotiate, because no data leaves the machine. The complete evaluation — hundreds of model calls across four models and three studies — ran on a single laptop with no cloud account. On disk and in memory the models are modest: the largest we used (`deepseek-r1:8b`) is roughly 5 GB and runs without a discrete GPU, while the sweet-spot model (`llama3.2`, 3 B) is about 2 GB. Taken together with the accuracy results, this puts a capable, private, auditable assistant within reach of a single researcher or a small lab with no cloud budget — precisely the resource-constrained scenario the course's small-business topic describes, transposed to a scientific setting. The cost of the approach is raw model capability relative to frontier cloud models, which we quantify (Table 5.2) rather than assume, and which grounding substantially offsets (Table 5.3).

---

## 6. Discussion

### 6.1 Determinism where reducible, notarization where irreducible

The evaluation makes the organizing thesis concrete. A quantum→AI pipeline chains several stochastic layers, and each is handled by the mechanism appropriate to it. At the **language-model layer**, randomness is optional: setting temperature to 0 makes generation bit-for-bit reproducible (Section 5.5, standard deviation exactly 0.00) with no measured cost to grounding or refusal — so we remove it. At the **quantum layer**, randomness is physical: quantum sampling cannot be "turned off," so instead of eliminating it we *pin* a specific result with a SHA-256 seal that a notary (LEON) recomputes on every read. A third, orthogonal mechanism — **retrieval grounding** — handles factual correctness. Separating these concerns is what lets one tool be private, reproducible, and auditable at once.

The distinction we are careful *not* to blur is that determinism buys **reproducibility, not truth**. Temperature 0 guarantees the same output for the same input; it does not make a wrong answer right. Correctness rests on grounding and, ultimately, on human peer review of the underlying science. Framed against FDA 21 CFR §11.10(e), reproducibility supports the "consistent intended performance" and record-integrity expectations, but is not by itself a claim of regulatory compliance — a boundary we state plainly rather than overselling.

### 6.2 Relationship to the four suggested topics

Although self-proposed, the project engages all four topics offered in the assignment. It explores the benefits, downsides, and risks of AI and quantum computing as **emerging technologies (Topic 2)**: the benefit is privacy and near-zero marginal cost, the risk is hallucination — which we measure and mitigate through the grounding contract — and quantum computing is the substrate the copilot serves. It provides a concrete framework for **choosing among LLMs, including local models (Topic 3)**, ranking four models by extraction accuracy, grounding, hallucination-resistance, speed, and cost; our recommendation of a 3 B general model over a larger reasoning model is a worked example of that selection. It demonstrates the utility of **randomness in deep learning (Topic 1)** through the repeated-run temperature study, showing that sampling randomness trades reproducibility for no accuracy gain on this task. And it targets **efficiency and effectiveness in a healthcare-research setting (Topic 4)** — faster, auditable review of oncology simulation results. The unifying contribution is methodological: a discipline for making a stochastic scientific pipeline trustworthy enough to explain and audit on a private, local machine.

## 7. Limitations [Section 5]

1. **Small local models are weaker than frontier cloud models.** `qwen3:0.6b` misses ~1 in 7 extraction fields; `deepseek-r1:1.5b` is weakest overall. The tool mitigates but does not eliminate this via grounding and refusal discipline.
2. **Determinism ≠ correctness.** Temperature 0 gives reproducibility, which supports (but does not by itself satisfy) 21 CFR §11 controls; factual correctness rests on the retrieval grounding, and regulatory validation is out of scope.
3. **Curated, English-only corpus.** The 19-passage BM25 corpus is intentionally small and hand-built; recall is bounded by what it contains, and BM25 is lexical.
4. **Evaluation scale.** Small question sets and single-configuration measurement (Section 5.6).
5. **Phase 3B is illustrative.** The 88-qubit record is a plan, not an executed run.

---

## 8. Future Work [Section 6]

1. **Chain-of-custody verification.** Seal verification is implemented (Section 4); a natural extension is to verify an entire batch of records at once and to check *lineage* across successive P1–P9 versions of the same run, flagging any break in the custody chain.
2. **Semantic retrieval.** Enable the `nomic-embed-text` path and compare lexical vs. semantic grounding as an added evaluation axis.
3. **Live literature.** Pull fresh PubMed abstracts per mutation at query time to keep Mode B current.
4. **Larger, repeated evaluation.** Expand question sets and average over runs for tighter confidence intervals; extend the temperature study to all four models on the enriched corpus.
5. **Batch narration & export.** Auto-generate reviewer-ready narratives for a folder of provenance records, with citations, for regulatory dossiers.
6. **Multimodal explanation.** Read the platform's 3D structure/energy figures alongside the record.

---

## 9. Conclusion

The SOLANGE Copilot shows that capable, grounded, auditable scientific assistance can run **entirely on a laptop**, with no data leaving the machine — the privacy dividend of local models. Empirically, a 3-billion-parameter model was the sweet spot, retrieval enrichment lifted a sub-billion model to perfect grounding, and sampling randomness cost consistency without buying accuracy. These results support a simple engineering doctrine for stochastic quantum→AI pipelines: **remove randomness where you can, notarize it where you can't, and ground every claim in a citable source.**

---

## References

Biomedical passages in the corpus are summarized from articles indexed in **PubMed**; each is cited below with its DOI.

1. Hassin O, Oren M. *Drugging p53 in cancer: one protein, many targets.* Nat Rev Drug Discov, 2022. DOI: 10.1038/s41573-022-00571-8
2. Wang H, et al. *Targeting p53 pathways: mechanisms, structures, and advances in therapy.* Signal Transduct Target Ther, 2023. DOI: 10.1038/s41392-023-01347-1
3. Chen X, et al. *Mutant p53 in cancer: from molecular mechanism to therapeutic modulation.* Cell Death Dis, 2022. DOI: 10.1038/s41419-022-05408-1
4. Mehri A, et al. *Dihydropyrimidine derivatives as MDM2 inhibitors.* Chem Biol Drug Des, 2023. DOI: 10.1111/cbdd.14399
5. Gourisankar S, et al. *Rewiring cancer drivers to activate apoptosis.* Nature, 2023. DOI: 10.1038/s41586-023-06348-2
6. Békés M, et al. *PROTAC targeted protein degraders: the past is prologue.* Nat Rev Drug Discov, 2022. DOI: 10.1038/s41573-021-00371-6
7. Dale B, et al. *Advancing targeted protein degradation for cancer therapy.* Nat Rev Cancer, 2021. DOI: 10.1038/s41568-021-00365-x
8. Skoulidis F, et al. *STK11/LKB1 mutations and PD-1 inhibitor resistance in KRAS-mutant lung adenocarcinoma.* Cancer Discov, 2018. DOI: 10.1158/2159-8290.CD-18-0099
9. Qian Y, et al. *MCT4-dependent lactate secretion suppresses antitumor immunity in LKB1-deficient lung adenocarcinoma.* Cancer Cell, 2023. DOI: 10.1016/j.ccell.2023.05.015
10. Gao Y, et al. *LKB1 in lung cancerigenesis: a serine/threonine kinase as tumor suppressor.* Protein Cell, 2011. DOI: 10.1007/s13238-011-1021-6
11. Kim J, et al. *CPS1 maintains pyrimidine pools and DNA synthesis in KRAS/LKB1-mutant lung cancer cells.* Nature, 2017. DOI: 10.1038/nature22359
12. Singh A, et al. *Small molecule inhibitor of NRF2 selectively intervenes therapeutic resistance in KEAP1-deficient NSCLC tumors.* ACS Chem Biol, 2016. DOI: 10.1021/acschembio.6b00651
13. Weiss-Sadan T, et al. *NRF2 activation induces NADH-reductive stress, providing a metabolic vulnerability in lung cancer.* Cell Metab, 2023. DOI: 10.1016/j.cmet.2023.01.012
14. Janes MR, et al. *Targeting KRAS mutant cancers with a covalent G12C-specific inhibitor.* Cell, 2018. DOI: 10.1016/j.cell.2018.01.006
15. Yamamoto S, et al. *WEE1 confers resistance to KRAS G12C inhibitors in non-small cell lung cancer.* Cancer Lett, 2024. DOI: 10.1016/j.canlet.2024.217414
16. Lebow ES, et al. *Analysis of tumor mutational burden, PFS, and local-regional control in locally advanced NSCLC treated with chemoradiation and durvalumab.* JAMA Netw Open, 2023. DOI: 10.1001/jamanetworkopen.2022.49591

**Technical references.** Sun Q, et al. *PySCF 2.0.* J Chem Phys, 2020. DOI: 10.1063/5.0006074 · Jordan P, Wigner E. *Über das Paulische Äquivalenzverbot.* Z Phys, 1928 · McClean JR, et al. *OpenFermion.* Quantum Sci Technol, 2020. DOI: 10.1088/2058-9565/ab8ebc · Ollama (local model runtime), ollama.com · Robertson S, Zaragoza H. *The Probabilistic Relevance Framework: BM25 and Beyond.* 2009.

---

## Appendix A — Reproducing the results

```bash
# From the repository root:
# 1. Run the app (no dependencies)
./start.sh                                 # open http://localhost:8000/app

# 2. Main benchmark (Table 5.2)
cd eval
SOLANGE_BACKEND=ollama python3 run_eval.py \
  --models llama3.2:latest qwen3:0.6b deepseek-r1:1.5b deepseek-r1:8b

# 3. Per-mutation drill-down (Table 5.4)
SOLANGE_BACKEND=ollama python3 run_eval.py --mutation "TP53 C275F" \
  --models llama3.2:latest qwen3:0.6b deepseek-r1:1.5b deepseek-r1:8b

# 4. Randomness study (Table 5.5)
SOLANGE_BACKEND=ollama python3 run_eval.py --temperature-sweep \
  --models qwen3:0.6b --temps 0.0 0.7 1.3 --repeats 5
```

Prerequisites: Python 3.10+, Ollama with the four models pulled (`ollama pull …`). With no models installed, everything still runs on a deterministic mock backend (numbers will not reflect model quality).
