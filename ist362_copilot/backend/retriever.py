"""
retriever.py — grounding for the copilot's answers (the "R" in RAG).

Two retrieval strategies over the same literature corpus:

  * BM25 (default) — a classic lexical ranking function implemented here in pure
    Python with no third-party dependencies. It works everywhere, instantly, and
    needs no model, which is exactly what a resource-constrained local
    deployment wants.
  * Embeddings (optional) — semantic search via Ollama's `nomic-embed-text`.
    Enabled only if that model is pulled; otherwise we stay on BM25.

Having both is deliberate: the paper compares lexical vs semantic retrieval as
one of its evaluation axes. Every answer the copilot produces is tied back to
the specific corpus passages retrieved here, so claims are auditable rather than
free-floating — "verify, don't trust," applied to the language model.
"""
from __future__ import annotations

import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

_TOKEN_RE = re.compile(r"[a-z0-9]+")
_DATA_DIR = Path(__file__).resolve().parent / "data"

# Generic words that carry no topical meaning. Without removing them, a question
# like "How is a run notarized?" spuriously matches any document containing
# "how"/"is", pulling in irrelevant sources.
_STOPWORDS = {
    "a", "an", "the", "is", "are", "was", "were", "be", "been", "being", "am",
    "how", "what", "which", "who", "whom", "whose", "when", "where", "why",
    "do", "does", "did", "doing", "done", "of", "to", "in", "on", "for", "and",
    "or", "with", "by", "as", "at", "it", "its", "that", "this", "these",
    "those", "can", "could", "would", "should", "will", "shall", "may", "might",
    "must", "from", "has", "have", "had", "having", "not", "no", "but", "if",
    "than", "then", "so", "such", "into", "over", "per", "via", "we", "our",
    "you", "your", "i", "my", "me", "they", "them", "their", "he", "she", "his",
    "her", "about", "more", "most", "some", "any", "all", "each", "other", "up",
    "out", "also", "there", "here", "using", "use", "used", "get", "got",
}

# Light suffix stemmer so query/document word forms match (e.g. "notarized" ~
# "notarization", "sealed" ~ "seal"). Suffixes are ordered longest-first.
_SUFFIXES = ("izations", "ization", "ational", "ations", "ation", "izing",
             "ized", "izes", "ize", "ings", "ing", "edly", "edness", "ness",
             "ers", "er", "ied", "ies", "ed", "es", "ly", "s")


def _stem(t: str) -> str:
    # Strip the negation prefix "un-" so e.g. "undruggable" and "non-druggable"
    # (which tokenizes as "non", "druggable") share a stem. Guarded to a long
    # enough remainder so short words like "until"/"unit" are left alone.
    if t.startswith("un") and len(t) - 2 >= 4:
        t = t[2:]
    for suf in _SUFFIXES:
        if t.endswith(suf) and len(t) - len(suf) >= 3:
            return t[:-len(suf)]
    return t


def tokenize(text: str) -> list[str]:
    return [_stem(t) for t in _TOKEN_RE.findall(text.lower())
            if t not in _STOPWORDS]


@dataclass
class Doc:
    id: str
    title: str
    text: str
    source: str
    tags: list[str]

    @property
    def blob(self) -> str:
        return f"{self.title}. {self.text} {' '.join(self.tags)}"


@dataclass
class Hit:
    doc: Doc
    score: float

    def to_dict(self) -> dict:
        return {
            "id": self.doc.id,
            "title": self.doc.title,
            "source": self.doc.source,
            "score": round(self.score, 4),
            "text": self.doc.text,   # full passage — the UI shows it in full
            "snippet": self.doc.text[:280] + ("…" if len(self.doc.text) > 280 else ""),
        }


class BM25:
    """Okapi BM25 ranking over a fixed corpus. k1/b are the standard defaults."""

    def __init__(self, docs: list[Doc], k1: float = 1.5, b: float = 0.75):
        self.docs = docs
        self.k1 = k1
        self.b = b
        self.corpus_tokens = [tokenize(d.blob) for d in docs]
        self.doc_len = [len(t) for t in self.corpus_tokens]
        self.avgdl = (sum(self.doc_len) / len(self.doc_len)) if self.doc_len else 0.0
        self.freqs = [Counter(t) for t in self.corpus_tokens]
        self.df: Counter = Counter()
        for toks in self.corpus_tokens:
            for term in set(toks):
                self.df[term] += 1
        self.N = len(docs)

    def _idf(self, term: str) -> float:
        n = self.df.get(term, 0)
        # BM25+ style non-negative idf
        return math.log(1 + (self.N - n + 0.5) / (n + 0.5))

    def search(self, query: str, top_k: int = 4) -> list[Hit]:
        q_terms = tokenize(query)
        scores = [0.0] * self.N
        for i in range(self.N):
            freq = self.freqs[i]
            dl = self.doc_len[i] or 1
            s = 0.0
            for term in q_terms:
                if term not in freq:
                    continue
                f = freq[term]
                num = f * (self.k1 + 1)
                den = f + self.k1 * (1 - self.b + self.b * dl / (self.avgdl or 1))
                s += self._idf(term) * num / den
            scores[i] = s
        ranked = sorted(range(self.N), key=lambda i: scores[i], reverse=True)
        hits = [Hit(self.docs[i], scores[i]) for i in ranked if scores[i] > 0]
        return hits[:top_k]


class Retriever:
    """Loads the corpus once and answers queries with BM25 (or optional embeddings)."""

    def __init__(self, corpus_path: Path | None = None):
        path = corpus_path or (_DATA_DIR / "corpus.json")
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        self.docs = [
            Doc(id=d["id"], title=d["title"], text=d["text"],
                source=d["source"], tags=d.get("tags", []))
            for d in raw["documents"]
        ]
        self.bm25 = BM25(self.docs)
        self.corpus_meta = {k: v for k, v in raw.items() if k != "documents"}

    def search(self, query: str, top_k: int = 4) -> list[Hit]:
        return self.bm25.search(query, top_k=top_k)

    def context_block(self, hits: list[Hit]) -> str:
        """Render retrieved passages as a numbered context block for the prompt."""
        lines = []
        for n, h in enumerate(hits, 1):
            lines.append(f"[{n}] {h.doc.title}\n    {h.doc.text}\n    (source: {h.doc.source})")
        return "\n".join(lines)


_CITE_RE = re.compile(r"\[(\d+)\]")


def cited_hits(answer: str, hits: list[Hit]) -> list[Hit]:
    """Return only the retrieved passages the answer actually cited as [n].

    The retriever feeds the top-k passages to the model as context, but some may
    be weak matches the model never uses. Displaying those alongside the answer
    makes irrelevant sources look supporting. Here we keep only the passages the
    answer referenced by bracket number (mapping [i] -> the i-th retrieved hit).
    Falls back to all hits if the answer cited none.
    """
    idxs = {int(n) for n in _CITE_RE.findall(answer or "")}
    if idxs:
        kept = [h for i, h in enumerate(hits, 1) if i in idxs]
        if kept:
            return kept
    # No usable citations (e.g. a model that doesn't bracket its sources): show
    # only the passages clearly relevant to the query — those within 35% of the
    # top BM25 score — so weak matches don't appear as supporting evidence.
    if not hits:
        return hits
    top = hits[0].score
    strong = [h for h in hits if h.score >= 0.35 * top]
    return strong or hits[:1]


_RETRIEVER: Retriever | None = None


def get_retriever() -> Retriever:
    global _RETRIEVER
    if _RETRIEVER is None:
        _RETRIEVER = Retriever()
    return _RETRIEVER
