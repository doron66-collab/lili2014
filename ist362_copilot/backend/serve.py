"""
serve.py — zero-dependency server for the SOLANGE Copilot.

Same endpoints as app.py (the FastAPI version) but built entirely on Python's
standard library (http.server), so it runs with **no pip install** — nothing to
compile, no third-party packages. This is the recommended way to run the demo,
especially on very new Python versions where prebuilt wheels may not yet exist.

    python3 serve.py            # then open http://localhost:8000/app

Endpoints
  GET  /                     service banner + active backend
  GET  /api/health           liveness + active backend/models
  GET  /api/models           models on the local Ollama server
  GET  /api/sample-run       bundled example provenance record
  POST /api/explain-run      Mode A: explain a P1-P9 provenance run
  POST /api/druggability     Mode B: explain non-druggability (grounded)
  POST /api/chat             grounded Q&A over the literature corpus
  GET  /app/ ...             the single-page frontend
"""
from __future__ import annotations

import json
import os
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

from llm_adapter import get_backend, OllamaBackend
from retriever import get_retriever
import prompts
import leon_verify

_DATA_DIR = Path(__file__).resolve().parent / "data"
_FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"
_MIME = {".html": "text/html", ".css": "text/css", ".js": "application/javascript",
         ".json": "application/json", ".svg": "image/svg+xml"}


# ── endpoint handlers (return plain dicts) ───────────────────────────────────
def h_root():
    return {"service": "SOLANGE Copilot API", "version": "0.1.0",
            "backend": get_backend().name,
            "note": "Local-LLM RAG copilot (stdlib server). Mode A: explain "
                    "provenance runs. Mode B: explain non-druggability."}


def h_health():
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
                        "to use real local models.")
    return info


def h_models():
    b = get_backend()
    return {"backend": b.name, "models": b.list_models()}


def h_sample_run():
    return json.loads((_DATA_DIR / "sample_run.json").read_text("utf-8"))


def h_sample_runs():
    data = json.loads((_DATA_DIR / "sample_runs.json").read_text("utf-8"))
    # Also send each record as a raw Python-serialized string so the browser can
    # display/edit it WITHOUT JSON.parse mangling floats (e.g. 0.0 -> 0), which
    # would otherwise break seal verification.
    for run in data.get("runs", []):
        run["record_text"] = json.dumps(run["record"], indent=2, ensure_ascii=False)
    return data


def h_explain_run(body: dict):
    run = body.get("run") or json.loads((_DATA_DIR / "sample_run.json").read_text("utf-8"))
    prompt = prompts.explain_run_prompt(json.dumps(run, indent=2), body.get("question"))
    res = get_backend().generate(prompt, model=body.get("model"),
                                 system=prompts.SYSTEM_PROMPT)
    return {"mode": "explain-run", "mutation_id": run.get("mutation_id"),
            "answer": res.text, **res.to_dict()}


def h_druggability(body: dict):
    mutation = body.get("mutation", "")
    r = get_retriever()
    hits = r.search(mutation + " undruggable non-druggable therapy",
                    top_k=int(body.get("top_k", 4)))
    prompt = prompts.druggability_prompt(mutation, r.context_block(hits),
                                         body.get("question"))
    res = get_backend().generate(prompt, model=body.get("model"),
                                 system=prompts.SYSTEM_PROMPT)
    return {"mode": "druggability", "mutation": mutation, "answer": res.text,
            "sources": [h.to_dict() for h in hits], **res.to_dict()}


def h_warmup(body: dict):
    """Load a model into memory with a trivial generation, so the *next* call
    measures generation time only (fair latency comparison, no load penalty)."""
    model = body.get("model")
    try:
        res = get_backend().generate("Reply with: OK", model=model, temperature=0)
        return {"ok": True, "model": res.model, "backend": res.backend}
    except Exception as e:
        return {"ok": False, "error": str(e)}


def h_verify_seal(body: dict):
    # Prefer the raw record_text (parsed here in Python, which preserves float
    # formatting like 0.0); fall back to a pre-parsed record object.
    text = body.get("record_text")
    if text is not None:
        try:
            record = json.loads(text)
        except Exception as e:
            return {"mode": "verify-seal", "verifiable": False,
                    "reason": f"Invalid JSON in record: {e}"}
    else:
        record = body.get("record") or {}
    return {"mode": "verify-seal", **leon_verify.recompute_seal(record)}


def h_chat(body: dict):
    q = body.get("question", "")
    r = get_retriever()
    hits = r.search(q, top_k=int(body.get("top_k", 4)))
    prompt = prompts.chat_prompt(r.context_block(hits), q)
    res = get_backend().generate(prompt, model=body.get("model"),
                                 system=prompts.SYSTEM_PROMPT)
    return {"mode": "chat", "answer": res.text,
            "sources": [h.to_dict() for h in hits], **res.to_dict()}


class Handler(BaseHTTPRequestHandler):
    def _send_json(self, obj, status=200):
        payload = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _send_file(self, path: Path):
        data = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", _MIME.get(path.suffix, "application/octet-stream"))
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self._send_json({}, status=204)

    def do_GET(self):
        route = urlparse(self.path).path
        try:
            if route == "/":
                return self._send_json(h_root())
            if route == "/api/health":
                return self._send_json(h_health())
            if route == "/api/models":
                return self._send_json(h_models())
            if route == "/api/sample-run":
                return self._send_json(h_sample_run())
            if route == "/api/sample-runs":
                return self._send_json(h_sample_runs())
            if route == "/app" or route == "/app/":
                return self._send_file(_FRONTEND_DIR / "index.html")
            if route.startswith("/app/"):
                target = (_FRONTEND_DIR / route[len("/app/"):]).resolve()
                if _FRONTEND_DIR.resolve() in target.parents and target.is_file():
                    return self._send_file(target)
                return self._send_file(_FRONTEND_DIR / "index.html")
            self._send_json({"error": "not found", "path": route}, status=404)
        except Exception as e:
            self._send_json({"error": str(e)}, status=500)

    def do_POST(self):
        route = urlparse(self.path).path
        try:
            length = int(self.headers.get("Content-Length", 0))
            body = json.loads(self.rfile.read(length) or b"{}") if length else {}
            if route == "/api/explain-run":
                return self._send_json(h_explain_run(body))
            if route == "/api/druggability":
                return self._send_json(h_druggability(body))
            if route == "/api/chat":
                return self._send_json(h_chat(body))
            if route == "/api/verify-seal":
                return self._send_json(h_verify_seal(body))
            if route == "/api/warmup":
                return self._send_json(h_warmup(body))
            self._send_json({"error": "not found", "path": route}, status=404)
        except Exception as e:
            self._send_json({"error": str(e)}, status=500)

    def log_message(self, fmt, *args):  # quieter console
        print("  %s - %s" % (self.address_string(), fmt % args))


def main():
    port = int(os.environ.get("PORT", "8000"))
    backend = get_backend()
    print("=" * 62)
    print("  SOLANGE Copilot — stdlib server (no dependencies)")
    print(f"  Active LLM backend: {backend.name}")
    print(f"  Open:  http://localhost:{port}/app")
    print(f"  API:   http://localhost:{port}/api/health")
    print("  Stop:  press Ctrl+C")
    print("=" * 62)
    ThreadingHTTPServer(("0.0.0.0", port), Handler).serve_forever()


if __name__ == "__main__":
    main()
