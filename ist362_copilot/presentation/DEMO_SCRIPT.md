# SOLANGE Copilot — Live Demo Script (IST 362)

**Goal:** ~8–10 minute talk + live demo. Lead with the *core* (two grounded modes
+ evaluation + LEON seal), keep the model-comparison tab as a short bonus.

---

## 0. Before you present (checklist)
- [ ] **Ollama running** (the app icon is open) and models pulled:
      `llama3.2:latest`, `qwen3:0.6b`, `deepseek-r1:1.5b`, `deepseek-r1:8b`.
- [ ] App started: from `ist362_copilot/` run `./start.sh`.
- [ ] Browser open at **http://localhost:8000/app**; badge shows **backend: ollama**.
- [ ] Zoom the browser to ~125% so the class can read it.
- [ ] Have the GitHub repo link ready on the last slide.

**One-line framing to open with:**
> "This is a private AI assistant for quantum-oncology results that runs
> entirely on my laptop — no data ever leaves the machine — and it never makes
> up an answer without citing its source."

---

## 1. The problem (30s, slide)
- SOLANGE produces results that are **complex** (quantum energies, 9-part
  provenance records with cryptographic seals) and **stochastic** (quantum
  sampling). Explaining and auditing them — over **sensitive genomic data** — is
  a bottleneck.
- A cloud LLM would mean sending that data out. We don't want that.

## 2. The idea + thesis (30s, slide)
- Emerging capability: **small local models are now good enough** to do this work
  on a laptop, privately.
- Thesis: **"Determinism where reducible, notarization where irreducible."**
  Remove randomness at the LLM layer (temperature 0); notarize it where it's
  physical (quantum sampling), via cryptographic seals.

## 3. DEMO — Mode A: Explain a run (90s)
1. Tab **A · Explain a run**. From the dropdown pick **"24 qubits — CASSCF(12,12)
   … (REAL run)"**.
2. Say: *"This is a real record from our HPC cluster."* Pick model **llama3.2**.
3. Click **Explain run**.
4. Read the output aloud — point out it names the mutation, the method, the exact
   energy (−437.76 Ha), the hardware, and the seal.
5. **Punchline:** *"Notice it reasoned that gate count 0 and depth 0 are expected
   because this used exact diagonalisation, not a circuit — it audits, it doesn't
   just summarize."*

## 4. DEMO — LEON seal verification (60s) ⭐ the memorable moment
1. Click **🔐 Verify seal (LEON)** → shows **✅ Seal VERIFIED**.
   Say: *"It recomputed the SHA-256 seal locally and it matches — the record is
   intact."*
2. In the record box, change the energy `-437.7612416122472` to e.g. `-430.0`.
3. Click **🔐 Verify seal** again → **❌ Seal MISMATCH — Changed field:
   p7_energy_ha**.
   Say: *"I changed one digit; it caught exactly which field. This is
   'verify, don't trust,' running entirely offline."*

## 5. DEMO — Mode B: Druggability (60s)
1. Tab **B · Druggability**. Enter **TP53 C275F**. Pick **qwen3:0.6b** (the 522 MB
   model) on purpose.
2. Click **Explain druggability**.
3. Say: *"This is a half-gigabyte model, yet the answer is grounded and cited —
   see the sources below, real PubMed references."* Scroll the sources.

## 6. DEMO — Compare models (45s, bonus)
1. Tab **⚖️ Compare models**. Question: *"Why is TP53 non-druggable?"* Check
   **llama3.2** and **qwen3:0.6b**. Click **Run comparison**.
2. Say: *"Same question, two models, side by side — with a warm-up so the timing
   is fair. This is how we picked our model."*

## 7. Evaluation (60s, slide — the 25-point section)
- Show the main results table. Key findings, in one breath:
  - **llama3.2 (3B) is the sweet spot** — perfect on every quality metric and the
    fastest; it ties the 8B model on quality at ~12× lower latency.
  - **Corpus enrichment lifted the 0.6B model to perfect grounding** — better
    retrieval beats a bigger model.
  - **Every model refused every out-of-context question** — grounding holds.
- Randomness study (Topic 1): *"Temperature 0 is deterministic — reproducible —
  with no accuracy cost. Randomness only adds variance. That's the language-model
  half of our thesis."*

## 8. Limitations + Future + close (30s, slide)
- Honest limits: small models are weaker; determinism ≠ correctness; small eval set.
- Future: chain-of-custody verification, semantic retrieval, live PubMed.
- **Close:** *"A capable, private, auditable assistant — on a laptop. Code and
  run instructions are here:"* → show the GitHub link.

---

## Timing cheat-sheet
| Segment | Time |
|---|---|
| Problem + idea + thesis | 1.5 min |
| Mode A + LEON | 2.5 min |
| Mode B + Compare | 2 min |
| Evaluation + randomness | 1.5 min |
| Limits + future + close | 0.5 min |
| **Total** | **~8 min** |

## If a demo fails (backup)
- Ollama not responding → badge shows `mock`; say *"falling back to the built-in
  deterministic backend"* and keep going; the flow still works.
- A model is slow → talk over it; that latency IS the point of the evaluation.
- Keep the paper's screenshots on a slide as a backup in case the laptop misbehaves.
