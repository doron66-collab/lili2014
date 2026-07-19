# SOLANGE Copilot

**A privacy-preserving, local-LLM copilot for auditable quantum-oncology results.**
IST 362 — Emerging Technologies · Doron Cohen · Claremont Graduate University, 2026

The copilot runs language models **entirely on your own machine** (via
[Ollama](https://ollama.com)) so that sensitive genomic, clinical, and
provenance data never leaves it — a privacy guarantee a cloud API cannot make.
It connects to the [SOLANGE](https://solange-platform.bio) quantum-simulation
platform for non-druggable NSCLC tumor-suppressor mutations (TP53, STK11, KEAP1).

Two capabilities, one shared engine (local LLM + retrieval-augmented grounding):

| Mode | What you give it | What it does |
|------|------------------|--------------|
| **A · Explain a run** | A SOLANGE P1–P9 provenance record (VQE / CASSCF results) | Turns it into a plain-English, auditable narrative and sanity-checks internal consistency |
| **B · Druggability** | A mutation (e.g. `TP53 C275F`) | Explains *why* it is non-druggable, grounded in retrieved literature, with citations |

Every answer reports which **backend** and **model** produced it and lists the
exact **sources** it used — "verify, don't trust," applied to the language model.

**LEON seal verification.** In Mode A, the **🔐 Verify seal** button recomputes the
record's P8 SHA-256 seal locally and compares it to the stored hash — a faithful
reproduction of the platform's LEON notary. A valid record shows ✅ *Seal
VERIFIED*; edit any sealed value (e.g. an energy) in the record box and re-verify
to see ❌ *Seal MISMATCH* — live tamper detection, entirely offline.

---

## Why this project (emerging-technologies angle)

The emerging capability on display is that **small, local models are now good
enough** to do useful scientific document work — extraction, summarization,
grounded Q&A — on a laptop, with **no data egress**. That reframes the classic
privacy/compliance trade-off (FDA 21 CFR §11) for regulated research data. The
paper evaluates how well several local models actually deliver on this promise.

---

## Architecture

```
ist362_copilot/
├── backend/
│   ├── app.py            FastAPI API (also serves the frontend at /app)
│   ├── llm_adapter.py    Ollama backend + deterministic mock; auto-detects
│   ├── retriever.py      BM25 retrieval (pure stdlib) + optional embeddings
│   ├── prompts.py        system prompt + per-mode templates (grounding rules)
│   └── data/
│       ├── corpus.json       literature corpus (real PubMed-sourced refs + method facts)
│       └── sample_run.json   example P1–P9 record = evaluation answer key
├── frontend/
│   └── index.html        single-page UI (SOLANGE styling), no build step
├── eval/
│   ├── eval_dataset.json extraction / grounding / trap questions with ground truth
│   └── run_eval.py       scores every model → Markdown comparison table
├── requirements.txt
└── README.md
```

**No third-party ML dependencies.** The adapter and retriever use only the
Python standard library, so the whole app — including the evaluation — runs even
with **zero models installed**, falling back to a deterministic mock backend.
Install Ollama to switch on real local models.

---

## Quick start

### Easiest path — no installation at all (recommended)

The server has a **zero-dependency** version built only on the Python standard
library. Nothing to `pip install`, nothing to compile — ideal on very new Python
versions (e.g. 3.14) where prebuilt packages may not exist yet.

```bash
cd ist362_copilot/backend
python3 serve.py
```
Then open **<http://localhost:8000/app>**. Stop with `Ctrl+C`.
Change the port with `PORT=8010 python3 serve.py`.

That's it — skip to *Install Ollama* below to switch on real local models. The
rest of this section (FastAPI, pip) is an optional alternative.

### Optional — the FastAPI version
```bash
cd ist362_copilot
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cd backend && uvicorn app:app --reload --port 8000
```

### Install Ollama and pull models
Download Ollama from <https://ollama.com/download>, then:
```bash
ollama pull llama3.2:latest      # 3B, general
ollama pull qwen3:0.6b           # tiny/fast baseline
ollama pull deepseek-r1:1.5b     # small reasoning
ollama pull deepseek-r1:8b       # largest / strongest
ollama pull nomic-embed-text     # optional: semantic retrieval
```
Ollama serves on `http://localhost:11434` automatically. The app detects it.

### 3. Run the backend (which also serves the UI)
```bash
cd backend
uvicorn app:app --reload --port 8000
```
Open **<http://localhost:8000/app>** for the UI, or **<http://localhost:8000/docs>**
for the interactive API.

- With Ollama running → real local models.
- Without Ollama → the app still runs on the mock backend (great for a quick
  look or CI). Force a backend with `SOLANGE_BACKEND=ollama` or `=mock`.

---

## Run the evaluation

```bash
cd eval
# head-to-head across your installed local models:
SOLANGE_BACKEND=ollama python run_eval.py \
  --models llama3.2:latest qwen3:0.6b deepseek-r1:1.5b deepseek-r1:8b
```
This prints and saves a comparison table (`results.json`, `results_table.md`):

| Model | Backend | Extraction acc. | Grounding acc. | Hallucination refusal | Mean latency (s) |
|---|---|---|---|---|---|
| … | ollama | … | … | … | … |

Metrics:
- **Extraction accuracy** — did the model read exact provenance fields correctly
  (ground truth = `sample_run.json`).
- **Grounding accuracy** — did it surface the right fact from retrieved literature.
- **Hallucination refusal rate** — for questions whose answer is *not* in the
  corpus, did it correctly refuse instead of inventing one.
- **Mean latency** — wall-clock per query, per model (the local cost/speed axis).

Offline sanity check (no models): `python run_eval.py`.

### Per-mutation evaluation (Mode B focus)

To compare the models **on one specific mutation** (Druggability mode only),
scoring how many key facts each model covers for that mutation and whether it
cites its sources:
```bash
SOLANGE_BACKEND=ollama python run_eval.py --mutation "TP53 C275F" \
  --models llama3.2:latest qwen3:0.6b deepseek-r1:1.5b deepseek-r1:8b
```
Defined mutations: `TP53 C275F`, `STK11 loss`, `KEAP1 loss`. Results save to
`results_<mutation>.json` / `.txt`.

### Randomness study (effect of sampling temperature)

Sampling **temperature** is the randomness knob of a language model: higher
temperature = more random token choices. This study measures how that randomness
affects faithfulness — grounding accuracy on answerable questions and refusal
rate on out-of-context questions:
```bash
SOLANGE_BACKEND=ollama python run_eval.py --temperature-sweep \
  --models llama3.2:latest --temps 0.0 0.5 1.0
```
Saves `results_temperature.json` / `.txt`. Expected trend: lower temperature →
more faithful, better-grounded, more consistent answers; higher temperature →
more variability and more hallucination. (The mock backend is deterministic, so
run this against Ollama to see the effect.)

---

## Configuration

| Env var | Default | Meaning |
|---|---|---|
| `SOLANGE_BACKEND` | `auto` | `auto` \| `ollama` \| `mock` |
| `OLLAMA_HOST` | `http://localhost:11434` | Ollama server URL |
| `SOLANGE_MODEL` | `llama3.2:latest` | default model when none is requested |
| `SOLANGE_LLM_TIMEOUT` | `120` | per-request timeout (seconds) |

---

## Limitations (see paper §5)

Local 0.6–8B models are weaker than frontier cloud models: they can still miss a
field or mis-cite. The BM25 corpus is intentionally small and curated. The
sample provenance record is a Phase 3A classical proxy, not live hardware output.
The mock backend is for plumbing only and does not represent model quality.

## Attribution

Biomedical passages in `corpus.json` are summarized from articles indexed in
**PubMed**; each carries its DOI. This tool is a research/education prototype and
is **not** a medical device or a source of clinical advice.
