# SOLANGE Copilot — Presentation Q&A Prep (IST 362)

Anticipated questions with short, defensible answers. Golden rule: **sell the
AI/software/auditability contribution; present the quantum layer honestly as the
substrate, not as new physics.** Never overclaim.

---

## A. Framing & scope

**Q: In one sentence, what is this?**
A local-LLM assistant that turns quantum-oncology provenance records into
plain-English, auditable narratives and explains mutation druggability — running
entirely on a laptop so no sensitive data leaves the machine.

**Q: What is the actual contribution — you didn't invent new quantum methods?**
Correct. The contribution is the **emerging-tech software layer**: showing that
small *local* models can do trustworthy scientific document work privately, and
adding cryptographic auditability (LEON) on top. The quantum results are the
substrate my tool reads and explains, not a claim of new physics.

**Q: Why is this "emerging technologies"?**
Two converging trends: capable open-weight LLMs now run on commodity hardware,
and privacy/compliance make *local* execution valuable. The project measures
whether that combination is good enough for real scientific use.

---

## B. Model choice & the LLM

**Q: Why these four models?**
Three size tiers (0.6B, 1.5B, 3B, 8B) across two families (general vs.
reasoning), so model choice becomes an experiment rather than a guess.

**Q: Why local instead of GPT-4/Claude?**
Privacy (no data egress over regulated genomic data), zero marginal cost, and
full reproducibility. A frontier cloud model would be more capable but can't
offer those three — which is the whole point.

**Q: Why did the 3B model win over the 8B?**
On this bounded, grounded task the 3B (`llama3.2`) matched the 8B on every
quality metric while running ~12× faster. For retrieval-grounded work, the model
mostly has to read and cite — raw capacity matters less than latency.

---

## C. RAG, grounding & hallucination

**Q: How do you stop the model from hallucinating?**
A strict system prompt ("use only the provided context, cite sources, refuse if
absent"), plus we *measure* it: trap questions with no answer in the corpus.
Every model refused 100% of them.

**Q: Isn't the corpus tiny (19 passages)?**
Yes — intentionally curated, and I disclose it as a limitation. Recall is bounded
by the corpus; the design goal was grounding discipline, not coverage. Future
work adds live PubMed retrieval.

**Q: Why BM25 and not embeddings?**
BM25 needs no model, runs instantly, and fits a resource-constrained local
deployment. The embedding path exists as an optional upgrade and a future eval
axis.

---

## D. Evaluation & methodology

**Q: How do you know the numbers are meaningful with such small test sets?**
I treat ≤ one-question differences as noise, and for the randomness study I ran
each condition 5× and report mean ± standard deviation. I state the small-n
threat to validity explicitly rather than hiding it.

**Q: What's your single most important result?**
Bigger isn't better: a 3B model beat an 8B on speed at equal quality, and
enriching the corpus lifted the 0.6B model to perfect grounding — retrieval
quality beats model size for this task.

**Q: How is "extraction accuracy" scored objectively?**
Against a ground-truth provenance record: the correct field values are known, so
scoring is deterministic keyword-set containment, not human judgement.

**Q: The randomness study looked flat — why include it?**
That *is* the finding: temperature 0 is bit-for-bit reproducible (std-dev 0.00)
with no accuracy cost, while higher temperature only adds variance. It's the
evidence behind "reduce randomness where it's reducible."

---

## E. LEON, seals & auditability

**Q: What does "Verify seal" actually do?**
It recomputes the record's P8 SHA-256 seal locally and compares it to the stored
hash — a faithful reproduction of the platform's own `build_p8_seal`. Edit any
sealed field and the hash diverges, so it names the tampered field. Fully offline.

**Q: Isn't the seal just re-hashing what you already have — what's the value?**
The value is "verify, don't trust": integrity is derived from a recomputable
proof, not from trusting the sender. It detects post-notarization tampering,
which is exactly what an auditor (or 21 CFR §11) needs.

**Q: Determinism means it's correct, right?**
No — and I'm careful about this. Determinism buys *reproducibility*, not truth.
Correctness comes from retrieval grounding; the underlying science still needs
peer review. I separate those three concerns deliberately.

---

## F. The quantum substrate (expect scrutiny — stay honest)

**Q: Is a 4-qubit CAS(2,2) toluene proxy a scientifically meaningful model of
TP53 C275F?**
It's a **minimal proof-of-concept proxy** for the sidechain's frontier orbitals,
not a claim about actual drug efficacy. Its role here is to give my *software* a
real, sealed provenance record to explain and audit — the copilot is
scale-agnostic and treats 4, 24, or 88 qubits identically.

**Q: Why do you even need a quantum computer — you ran 24 qubits on a GPU?**
Exactly — 4 and 24 qubits are classical benchmarks (exact statevector). The
quantum-hardware need is for the *large* active spaces (e.g. 88 qubits / 44
correlated electrons) that exceed classical exact methods; that Phase 3B run is
marked as planned, not executed.

**Q: Where does the classical method actually break down?**
Classical exact treatment of strong correlation becomes intractable as the active
space grows (CCSD(T) around ~18 correlated electrons in our setting); the full
C275F problem (44e / 88 qubits) is beyond that, which motivates the quantum
approach — but I don't claim to have run it yet.

**Q: What is Jordan-Wigner / CASSCF(12,12) → 24 qubits?**
CASSCF(12,12) means 12 electrons in 12 spatial orbitals; Jordan-Wigner maps each
spin-orbital to one qubit, so 12 orbitals → 24 spin-orbitals → 24 qubits. The
copilot verifies this consistency automatically.

**Q: Are these mutations really "non-druggable"?**
"Historically hard to drug" is the accurate phrasing, and it's grounded in the
peer-reviewed literature the copilot cites (e.g. p53 as a transcription factor
lacking a compact pocket), not my opinion.

---

## G. Privacy & compliance

**Q: How is privacy actually guaranteed?**
Nothing leaves the machine: models run via local Ollama on `localhost`, there is
no network call, no API key, no cloud account. That's verifiable — you can pull
the network cable and it still works.

**Q: Does this make you 21 CFR §11 compliant?**
No — I'm precise about this. Reproducibility and seal re-verification *support*
the record-integrity and "consistent performance" expectations, but full
compliance is a validation process out of scope for a course project.

---

## H. Limitations & honesty

**Q: What's the biggest weakness?**
Small local models are weaker than frontier models — they can miss a field or
mis-phrase. Grounding and refusal discipline mitigate it, and I quantify rather
than assume the gap.

**Q: What would you do with more time?**
Chain-of-custody verification across record versions, semantic retrieval,
live PubMed, and a larger repeated evaluation for tighter confidence intervals.

**Q: What did you learn?**
That trustworthiness in a stochastic AI+quantum pipeline is an engineering
discipline: remove randomness where you can (temperature 0), notarize it where
you can't (seals), and ground every claim in a citable source.

---

## I. Presentation delivery — meta

**Q (to self): How do I open?**
"This is a private AI assistant for quantum-oncology results that runs entirely
on my laptop — no data leaves the machine — and it never answers without citing
its source."

**Q (to self): If the live demo fails?**
Stay calm: the app falls back to a deterministic mock (say so), and I have
screenshots on the slides as backup. A slow model *is* the point of the
evaluation — talk over it.

**Q (to self): If asked something I don't know (esp. quantum)?**
Answer honestly: "That's beyond what I've validated — here's what I can defend,
and here's what would need peer review." Honesty reads as competence to an expert.
