# SOLANGE Copilot — Final Presentation Script (IST 362)

**Goal:** ~9-9.5 minute talk + live demo, slide-by-slide with exact narration
and exact demo steps.

---

## Before you present (checklist)
- [ ] **Ollama running** and models pulled: `llama3.2:latest`, `qwen3:0.6b`,
      `deepseek-r1:1.5b`, `deepseek-r1:8b`.
- [ ] **Restart the server** if it was already running, so it picks up the
      latest retrieval fix: `cd ist362_copilot && ./start.sh`
- [ ] Browser open at **http://localhost:8000/app**; badge shows **backend: ollama**.
- [ ] Zoom the browser to ~125% so the class can read it.
- [ ] Presentation deck open, GitHub repo link ready on the last slide.

---

## SLIDE 1 — Title (10s)
> "This is a private AI assistant for quantum-oncology results that runs
> entirely on my laptop — no data ever leaves the machine — and it never
> makes up an answer without citing its source."

**→ advance to Slide 2**

---

## SLIDE 2 — The problem (30s)
> "SOLANGE produces results that are complex — nine-part provenance records,
> quantum energies, cryptographic seals — and stochastic, because quantum
> sampling is non-reproducible by default. And it's sensitive genomic data,
> so a cloud LLM is off the table — we don't want to send that out."

**→ advance to Slide 3**

---

## SLIDE 3 — The idea + thesis (30s)
> "The emerging capability here is that small local models are now good
> enough to do this work — extract, summarize, answer — entirely on a
> laptop, with nothing leaving the machine. Every answer is grounded in
> retrieved sources and cites them; it refuses rather than guesses."
>
> "Our thesis: determinism where reducible, notarization where irreducible.
> We remove randomness at the LLM layer with temperature zero. Where
> randomness is physical — quantum sampling — we don't remove it, we
> notarize it with cryptographic seals."

**→ advance to Slide 4**

---

## SLIDE 4 — Model choice (15s)
> "We tested four local models across three size tiers — from qwen3 at half
> a gigabyte up to an 8-billion-parameter reasoning model. llama3.2, the
> 3-billion general model, turned out to be our sweet spot — and model
> choice is itself one of our findings, which I'll show live."

**→ advance to Slide 5**

---

## SLIDE 5 — How it works (15s)
> "The pipeline is five steps: query, retrieve over 19 real PubMed-sourced
> references with BM25, ground under a strict contract — use only the
> retrieved context — generate locally via Ollama, and answer with cited
> sources."

**→ advance to Slide 6**

---

## SLIDE 6 — Capabilities (10s)
> "Two modes: explaining a provenance run, and explaining why a mutation is
> non-druggable — plus live LEON seal verification, and it handles real
> records at 4, 24, and 88 qubits."

**→ switch to the browser (localhost:8000/app)**

---

## DEMO 1 — Mode A: Explain a run (90s)

**Exactly what to demo:**
1. Tab **A · Explain a run**.
2. Dropdown → select **"24 qubits — CASSCF(12,12) … (REAL run)"**.
3. Model → **llama3.2:latest**.
4. Click **Explain run**.
5. Read the output aloud — confirm it names the mutation, method, exact
   energy (−437.76 Ha), hardware, and seal.

**Talking points:**
> "This is a real record from our HPC cluster."
>
> "Notice it reasoned that gate count 0 and depth 0 are expected because
> this used exact diagonalisation, not a circuit — it audits, it doesn't
> just summarize."

---

## DEMO 2 — LEON seal verification (60s) ⭐ the memorable moment

**Exactly what to demo:**
1. Click **🔐 Verify seal (LEON)** → shows **✅ Seal VERIFIED**.
2. In the record text box, change `-437.7612416122472` to any other number,
   e.g. `-430.0`.
3. Click **🔐 Verify seal** again → shows **❌ Seal MISMATCH — Changed field:
   p7_energy_ha**.

**Talking points:**
> "It recomputed the SHA-256 seal locally and it matches — the record is
> intact."
>
> "I changed one digit; it caught exactly which field. This is 'verify,
> don't trust,' running entirely offline."

---

## DEMO 3 — Mode B: Druggability (60s)

**Exactly what to demo:**
1. Tab **B · Druggability**.
2. Input field → type **TP53 C275F**.
3. Model → **qwen3:0.6b** (on purpose — the smallest model).
4. Click **Explain druggability**.
5. Scroll down to show the sources list.

**Talking point:**
> "This is a half-gigabyte model, yet the answer is grounded and cited —
> see the sources below, real PubMed references."

---

## DEMO 4 — Compare models (60-75s, includes the limitation moment) ⭐

**Exactly what to demo:**
1. Tab **⚖️ Compare models**.
2. Question field → type **"Why is TP53 non-druggable?"** (exactly this
   phrasing, not "undruggable").
3. Check only **llama3.2:latest** and **qwen3:0.6b** (two models only, not
   all four — keeps the contrast clean).
4. Click **Run comparison**.
5. Once results appear, point at both cards together: llama3.2's grounded
   answer (with a bracketed source citation) **versus** qwen3's refusal.
6. **Critical**: explicitly point at the line "**grounded in 4 source(s)**"
   under qwen3's answer before explaining the refusal.

**Talking points (the limitation moment):**
> "Same question, both models retrieved the exact same four sources — you
> can see it says 'grounded in 4 sources' right there under qwen3's answer.
> llama3.2 used them and cited its source. qwen3 — the half-gigabyte model —
> refused anyway, even though the right information was sitting right in
> front of it."
>
> "This isn't a retrieval bug — it's a capability limit we documented in the
> paper: smaller models can be over-cautious, misjudging when to apply the
> 'refuse if not grounded' rule even when grounding is present. It's exactly
> why model choice mattered, and exactly why we evaluate reliability
> per-model instead of assuming one grounding prompt works equally well
> everywhere."

---

**→ switch back to slides, advance to Slide 8**

## SLIDE 8 — Evaluation (40s)
> "Bigger isn't better: llama3.2 is perfect on every quality metric and the
> fastest — it ties the 8-billion model at about twelve times lower
> latency. Retrieval beats size: enriching the corpus lifted the
> 0.6-billion model to 100% grounding. And every model refused 100% of the
> formal out-of-context trap questions — no hallucination on those."

**→ advance to Slide 9**

---

## SLIDE 9 — Randomness · Topic 1 (30s)
> "At temperature zero, standard deviation is exactly 0.00 — fully
> deterministic, bit-for-bit reproducible — with 100% refusal at every
> temperature tested. Higher temperature only adds variance, no accuracy
> gain. That's the language-model half of the thesis: reduce randomness
> where it's reducible."

**→ advance to Slide 10**

---

## SLIDE 10 — Limitations · Future (25s)
> "Honestly: small local models are weaker than frontier models — you just
> saw that live with qwen3's over-caution. Determinism buys reproducibility
> but not correctness, our evaluation set is small, and the 88-qubit Phase
> 3B record is a plan, not yet run. Future work: chain-of-custody
> verification, semantic retrieval, live PubMed lookup, and a larger
> repeated evaluation."

**→ advance to Slide 11**

---

## SLIDE 11 — Close (10s)
> "A capable, private, auditable assistant — on a laptop. Code and run
> instructions are here:" *(point at the GitHub link)*

---

## Timing cheat-sheet

| Segment | Time |
|---|---|
| Problem + idea + thesis (Slides 2-3) | 1 min |
| Model choice + how it works + capabilities (Slides 4-6) | 40s |
| Mode A + LEON (Demo 1-2) | 2.5 min |
| Mode B + Compare (Demo 3-4) | 2-2.25 min |
| Evaluation + randomness (Slides 8-9) | 1.2 min |
| Limits + future + close (Slides 10-11) | 35s |
| **Total** | **~9-9.5 min** |

## If a demo fails (backup)
- Ollama not responding → badge shows `mock`; say *"falling back to the
  built-in deterministic backend"* and keep going.
- A model is slow → talk over it; that latency IS the point of the
  evaluation.
- Keep the paper's screenshots on a slide as backup in case the laptop
  misbehaves.
