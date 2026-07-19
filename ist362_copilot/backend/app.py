"""
app.py — FastAPI backend for the SOLANGE Copilot (IST 362 midterm project).

Endpoints
  GET  /                      service banner + active backend
  GET  /api/health           liveness + which LLM backend/model is active
  GET  /api/models           models available on the local Ollama server
  GET  /api/sample-run       the bundled example provenance record
  POST /api/explain-run      Mode A: explain a P1-P9 provenance run
  POST /api/druggability     Mode B: explain why a mutation is non-druggable
  POST /api/chat             grounded Q&A over the literature corpus

Every response reports the backend ("ollama" or "mock") and model that produced
it, plus the retrieved sources, so results stay auditable end-to-end.
"""
from __future__ import annotations

import json
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from llm_adapter import get_backend, OllamaBackend
from retriever import get_retriever
import prompts
import leon_verify

_DATA_DIR = Path(__file__).resolve().parent / "data"
_FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"

app = FastAPI(
    title="SOLANGE Copilot API",
    description="Privacy-preserving local-LLM copilot for auditable quantum-oncology "
                "results — IST 362 Emerging Technologies · Doron Cohen · CGU 2026",
    version="0.1.0",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
)


# ── request models ───────────────────────────────────────────────────────────
class ExplainRunRequest(BaseModel):
    run: dict | None = None          # a P1-P9 record; if None, use bundled sample
    question: str | None = None
    model: str | None = None


class DruggabilityRequest(BaseModel):
    mutation: str
    question: str | None = None
    model: str | None = None
    top_k: int = 4


class ChatRequest(BaseModel):
    question: str
    model: str | None = None
    top_k: int = 4


class VerifySealRequest(BaseModel):
    record: dict | None = None
    record_text: str | None = None


# ── read endpoints ───────────────────────────────────────────────────────────
@app.get("/")
def root():
    b = get_backend()
    return {
        "service": "SOLANGE Copilot API",
        "version": "0.1.0",
        "backend": b.name,
        "note": "Local-LLM RAG copilot. Mode A: explain provenance runs. "
                "Mode B: explain mutation non-druggability.",
    }


@app.get("/api/health")
def health():
    b = get_backend()
    info = {"status": "ok", "backend": b.name}
    if isinstance(b, OllamaBackend):
        try:
            info["models"] = b.list_models()
        except Exception as e:
            info["models_error"] = str(e)
    else:
        info["models"] = b.list_models()
        info["hint"] = ("Running on the deterministic mock backend. Start Ollama "
                        "and set SOLANGE_BACKEND=ollama (or leave =auto) to use "
                        "real local models.")
    return info


@app.get("/api/models")
def models():
    b = get_backend()
    return {"backend": b.name, "models": b.list_models()}


@app.get("/api/sample-run")
def sample_run():
    return json.loads((_DATA_DIR / "sample_run.json").read_text(encoding="utf-8"))


@app.get("/api/sample-runs")
def sample_runs():
    """A set of example provenance records at increasing scale (4 -> 24 -> 88 qubits)."""
    data = json.loads((_DATA_DIR / "sample_runs.json").read_text(encoding="utf-8"))
    # raw string so the browser doesn't mangle floats (0.0 -> 0) before verify
    for run in data.get("runs", []):
        run["record_text"] = json.dumps(run["record"], indent=2, ensure_ascii=False)
    return data


# ── Mode A: explain a provenance run ─────────────────────────────────────────
@app.post("/api/explain-run")
def explain_run(req: ExplainRunRequest):
    run = req.run or json.loads((_DATA_DIR / "sample_run.json").read_text("utf-8"))
    run_json = json.dumps(run, indent=2)
    prompt = prompts.explain_run_prompt(run_json, req.question)
    backend = get_backend()
    result = backend.generate(prompt, model=req.model, system=prompts.SYSTEM_PROMPT)
    return {
        "mode": "explain-run",
        "mutation_id": run.get("mutation_id"),
        "answer": result.text,
        **result.to_dict(),
    }


# ── Mode B: explain non-druggability, grounded in literature ─────────────────
@app.post("/api/druggability")
def druggability(req: DruggabilityRequest):
    retriever = get_retriever()
    hits = retriever.search(req.mutation + " undruggable non-druggable therapy",
                            top_k=req.top_k)
    context = retriever.context_block(hits)
    prompt = prompts.druggability_prompt(req.mutation, context, req.question)
    backend = get_backend()
    result = backend.generate(prompt, model=req.model, system=prompts.SYSTEM_PROMPT)
    return {
        "mode": "druggability",
        "mutation": req.mutation,
        "answer": result.text,
        "sources": [h.to_dict() for h in hits],
        **result.to_dict(),
    }


# ── LEON: live seal verification ─────────────────────────────────────────────
@app.post("/api/verify-seal")
def verify_seal(req: VerifySealRequest):
    """Re-attest a provenance record's P8 seal locally (verify, don't trust)."""
    if req.record_text is not None:
        try:
            record = json.loads(req.record_text)
        except Exception as e:
            return {"mode": "verify-seal", "verifiable": False,
                    "reason": f"Invalid JSON in record: {e}"}
    else:
        record = req.record or {}
    return {"mode": "verify-seal", **leon_verify.recompute_seal(record)}


# ── grounded free-form Q&A ───────────────────────────────────────────────────
@app.post("/api/chat")
def chat(req: ChatRequest):
    retriever = get_retriever()
    hits = retriever.search(req.question, top_k=req.top_k)
    context = retriever.context_block(hits)
    prompt = prompts.chat_prompt(context, req.question)
    backend = get_backend()
    result = backend.generate(prompt, model=req.model, system=prompts.SYSTEM_PROMPT)
    return {
        "mode": "chat",
        "answer": result.text,
        "sources": [h.to_dict() for h in hits],
        **result.to_dict(),
    }


# ── serve the frontend (so `uvicorn app:app` gives a one-URL demo) ───────────
if _FRONTEND_DIR.exists():
    app.mount("/app", StaticFiles(directory=str(_FRONTEND_DIR), html=True),
              name="frontend")
