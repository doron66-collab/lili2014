"""
llm_adapter.py — pluggable local-LLM backend for the SOLANGE Copilot.

Design goal (IST 362, emerging technologies): the *emerging* capability we are
demonstrating is running capable language models entirely **locally**, so that
sensitive genomic / clinical / provenance data never leaves the machine. That
privacy guarantee is something a cloud API cannot offer, and it is the core
argument of the accompanying paper.

Two interchangeable backends implement one interface:

  * OllamaBackend — talks to a local Ollama server (http://localhost:11434).
    This is what runs on the author's machine for the graded demo and the
    quantitative evaluation.
  * MockBackend  — deterministic, dependency-free responses so the whole
    application (API, frontend, and evaluation harness) runs and can be tested
    on any machine, including CI, with no model installed.

`get_backend()` auto-detects: if a live Ollama server is reachable it uses it,
otherwise it transparently falls back to the mock. The active backend is always
reported through the API so the UI/paper can state exactly what produced a
result (honest provenance — the same principle SOLANGE applies to its science).
"""
from __future__ import annotations

import json
import os
import re
import time
import urllib.request
import urllib.error
from dataclasses import dataclass, field

OLLAMA_HOST = os.environ.get("OLLAMA_HOST", "http://localhost:11434")
DEFAULT_MODEL = os.environ.get("SOLANGE_MODEL", "llama3.2:latest")
REQUEST_TIMEOUT = float(os.environ.get("SOLANGE_LLM_TIMEOUT", "120"))

# deepseek-r1 (and other reasoning models) emit their chain-of-thought wrapped
# in <think>...</think>. For a grounded, auditable answer we keep only the final
# response and surface the reasoning separately.
_THINK_RE = re.compile(r"<think>(.*?)</think>", re.DOTALL)


@dataclass
class LLMResult:
    """One generation, plus the metadata the evaluation harness needs."""
    text: str
    model: str
    backend: str
    latency_s: float
    reasoning: str = ""            # extracted <think> content, if any
    raw: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "text": self.text,
            "model": self.model,
            "backend": self.backend,
            "latency_s": round(self.latency_s, 3),
            "reasoning": self.reasoning,
        }


def _split_reasoning(text: str) -> tuple[str, str]:
    """Return (clean_answer, reasoning). Strips <think> blocks from the answer."""
    reasoning = "\n".join(m.strip() for m in _THINK_RE.findall(text)).strip()
    clean = _THINK_RE.sub("", text).strip()
    return clean, reasoning


class LLMBackend:
    name = "base"

    def generate(self, prompt: str, *, model: str | None = None,
                 system: str | None = None, temperature: float = 0.2) -> LLMResult:
        raise NotImplementedError

    def list_models(self) -> list[str]:
        raise NotImplementedError

    def available(self) -> bool:
        return True


class OllamaBackend(LLMBackend):
    """Calls a local Ollama server via its native HTTP API — no SDK, no cloud."""
    name = "ollama"

    def __init__(self, host: str = OLLAMA_HOST):
        self.host = host.rstrip("/")

    def available(self) -> bool:
        try:
            self._get("/api/tags", timeout=2)
            return True
        except Exception:
            return False

    def list_models(self) -> list[str]:
        data = self._get("/api/tags", timeout=5)
        return sorted(m["name"] for m in data.get("models", []))

    def generate(self, prompt: str, *, model: str | None = None,
                 system: str | None = None, temperature: float = 0.2) -> LLMResult:
        model = model or DEFAULT_MODEL
        payload = {
            "model": model,
            "prompt": prompt,
            "stream": False,
            "options": {"temperature": temperature},
        }
        if system:
            payload["system"] = system
        t0 = time.perf_counter()
        data = self._post("/api/generate", payload, timeout=REQUEST_TIMEOUT)
        latency = time.perf_counter() - t0
        clean, reasoning = _split_reasoning(data.get("response", ""))
        return LLMResult(text=clean, model=model, backend=self.name,
                         latency_s=latency, reasoning=reasoning, raw=data)

    def embed(self, text: str, model: str = "nomic-embed-text") -> list[float]:
        """Optional semantic-search path. Requires `ollama pull nomic-embed-text`."""
        data = self._post("/api/embeddings", {"model": model, "prompt": text},
                          timeout=REQUEST_TIMEOUT)
        return data.get("embedding", [])

    # ── low-level HTTP helpers (stdlib only) ─────────────────────────────
    def _get(self, path: str, timeout: float) -> dict:
        req = urllib.request.Request(self.host + path)
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8"))

    def _post(self, path: str, body: dict, timeout: float) -> dict:
        req = urllib.request.Request(
            self.host + path,
            data=json.dumps(body).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8"))


class MockBackend(LLMBackend):
    """
    Deterministic stand-in so the app + evaluation run with zero models installed.
    It is intentionally *extractive*: it echoes grounded facts from the prompt /
    retrieved context. That keeps demos honest (it never invents citations) and
    lets the evaluation harness exercise the full pipeline offline.
    """
    name = "mock"
    MODELS = ["mock-small", "mock-large"]

    def list_models(self) -> list[str]:
        return list(self.MODELS)

    def generate(self, prompt: str, *, model: str | None = None,
                 system: str | None = None, temperature: float = 0.2) -> LLMResult:
        model = model or "mock-small"
        t0 = time.perf_counter()
        text = self._synthesize(prompt)
        time.sleep(0.01)  # keep latency numbers non-zero for the harness
        return LLMResult(text=text, model=model, backend=self.name,
                         latency_s=time.perf_counter() - t0)

    @staticmethod
    def _synthesize(prompt: str) -> str:
        """Pull the most answer-bearing lines out of the assembled prompt."""
        facts: list[str] = []
        for line in prompt.splitlines():
            s = line.strip("-• \t")
            low = s.lower()
            if any(k in low for k in (
                "energy", "qubit", "mutation", "casscf", "vqe", "hartree",
                "undruggable", "non-druggable", "nrf2", "keap1", "stk11",
                "tp53", "notariz", "p8", "seal", "converged")):
                facts.append(s)
        head = facts[:6] if facts else [
            "No grounded facts were found in the provided context."]
        body = " ".join(head)
        return ("[mock backend — deterministic, no model installed] "
                "Grounded summary from provided context: " + body)


_BACKEND_CACHE: LLMBackend | None = None


def get_backend(force: str | None = None) -> LLMBackend:
    """
    Resolve the active backend. Order:
      1. explicit `force` argument ("ollama" | "mock")
      2. SOLANGE_BACKEND env var
      3. auto-detect a live Ollama server, else mock.
    Cached so we don't probe the socket on every request.
    """
    global _BACKEND_CACHE
    choice = (force or os.environ.get("SOLANGE_BACKEND", "auto")).lower()

    if choice == "ollama":
        return OllamaBackend()
    if choice == "mock":
        return MockBackend()

    if _BACKEND_CACHE is not None and force is None:
        return _BACKEND_CACHE

    ollama = OllamaBackend()
    backend: LLMBackend = ollama if ollama.available() else MockBackend()
    _BACKEND_CACHE = backend
    return backend
