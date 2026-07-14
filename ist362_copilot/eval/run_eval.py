"""
run_eval.py — quantitative evaluation of the SOLANGE Copilot.

Produces the numbers the paper's evaluation section [25 pts] is built on. Three
machine-scorable tasks, each with a built-in ground truth:

  1. Extraction accuracy (Mode A): does the model read exact fields
     (basis set, qubit count, energy, notary, …) correctly out of a provenance
     record? Ground truth = the record itself, so scoring is objective.

  2. Grounding accuracy (Mode B / chat): does the model surface the correct
     fact from the retrieved literature?

  3. Hallucination resistance (traps): for questions whose answer is NOT in the
     corpus, does the model correctly refuse instead of inventing an answer?

Run against every installed model to produce a head-to-head comparison table:

    # all local models, real Ollama:
    SOLANGE_BACKEND=ollama python run_eval.py --models llama3.2:latest qwen3:0.6b deepseek-r1:1.5b deepseek-r1:8b

    # offline sanity check with the deterministic mock:
    python run_eval.py

Results are printed as a Markdown table and saved to eval/results.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

# make the backend package importable when run from anywhere
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "backend"))

from llm_adapter import get_backend           # noqa: E402
from retriever import get_retriever           # noqa: E402
import prompts                                 # noqa: E402

_HERE = Path(__file__).resolve().parent
_DATA = _ROOT / "backend" / "data"


def contains_any(text: str, options: list[str]) -> bool:
    low = text.lower()
    return any(opt.lower() in low for opt in options)


def eval_extraction(backend, model, dataset, run_record):
    """Mode A: ask each field question against the provenance record."""
    run_json = json.dumps(run_record, indent=2)
    rows, correct, lat = [], 0, []
    for f in dataset["extraction"]["fields"]:
        prompt = prompts.explain_run_prompt(run_json, f["question"])
        res = backend.generate(prompt, model=model, system=prompts.SYSTEM_PROMPT)
        ok = contains_any(res.text, f["expect_any"])
        correct += ok
        lat.append(res.latency_s)
        rows.append({"field": f["name"], "ok": ok, "answer": res.text[:120]})
    n = len(dataset["extraction"]["fields"])
    return {"accuracy": correct / n, "n": n, "mean_latency_s": _mean(lat), "rows": rows}


def eval_grounding(backend, model, dataset):
    """Mode B / chat: questions answerable from the corpus."""
    retriever = get_retriever()
    rows, correct, lat = [], 0, []
    for item in dataset["grounding"]["items"]:
        hits = retriever.search(item["question"], top_k=4)
        prompt = prompts.chat_prompt(retriever.context_block(hits), item["question"])
        res = backend.generate(prompt, model=model, system=prompts.SYSTEM_PROMPT)
        ok = contains_any(res.text, item["expect_any"])
        correct += ok
        lat.append(res.latency_s)
        rows.append({"q": item["question"][:60], "ok": ok})
    n = len(dataset["grounding"]["items"])
    return {"accuracy": correct / n, "n": n, "mean_latency_s": _mean(lat), "rows": rows}


def eval_traps(backend, model, dataset):
    """Out-of-context questions: correct behavior is to refuse."""
    retriever = get_retriever()
    markers = dataset["traps"]["refusal_markers"]
    rows, refused, lat = [], 0, []
    for item in dataset["traps"]["items"]:
        hits = retriever.search(item["question"], top_k=4)
        prompt = prompts.chat_prompt(retriever.context_block(hits), item["question"])
        res = backend.generate(prompt, model=model, system=prompts.SYSTEM_PROMPT)
        ok = contains_any(res.text, markers)
        refused += ok
        lat.append(res.latency_s)
        rows.append({"q": item["question"][:60], "refused": ok, "answer": res.text[:100]})
    n = len(dataset["traps"]["items"])
    return {"refusal_rate": refused / n, "n": n, "mean_latency_s": _mean(lat), "rows": rows}


def _mean(xs):
    return round(sum(xs) / len(xs), 3) if xs else 0.0


def run_for_model(backend, model, dataset, run_record):
    print(f"\n=== Evaluating model: {model or '(default)'} on backend '{backend.name}' ===")
    ext = eval_extraction(backend, model, dataset, run_record)
    grd = eval_grounding(backend, model, dataset)
    trp = eval_traps(backend, model, dataset)
    all_lat = ext["mean_latency_s"] + grd["mean_latency_s"] + trp["mean_latency_s"]
    return {
        "model": model or "(default)",
        "backend": backend.name,
        "extraction_accuracy": round(ext["accuracy"], 3),
        "grounding_accuracy": round(grd["accuracy"], 3),
        "hallucination_refusal_rate": round(trp["refusal_rate"], 3),
        "mean_latency_s": round(all_lat / 3, 3),
        "detail": {"extraction": ext, "grounding": grd, "traps": trp},
    }


def markdown_table(results: list[dict]) -> str:
    hdr = ("| Model | Backend | Extraction acc. | Grounding acc. | "
           "Hallucination refusal | Mean latency (s) |\n"
           "|---|---|---|---|---|---|\n")
    rows = ""
    for r in results:
        rows += (f"| {r['model']} | {r['backend']} | {r['extraction_accuracy']:.2f} | "
                 f"{r['grounding_accuracy']:.2f} | {r['hallucination_refusal_rate']:.2f} | "
                 f"{r['mean_latency_s']:.2f} |\n")
    return hdr + rows


def main():
    ap = argparse.ArgumentParser(description="Evaluate the SOLANGE Copilot.")
    ap.add_argument("--models", nargs="*", default=[None],
                    help="model names to compare (default: backend default)")
    ap.add_argument("--out", default=str(_HERE / "results.json"))
    args = ap.parse_args()

    dataset = json.loads((_HERE / "eval_dataset.json").read_text("utf-8"))
    run_record = json.loads((_DATA / "sample_run.json").read_text("utf-8"))
    backend = get_backend()

    results = [run_for_model(backend, m, dataset, run_record) for m in args.models]

    table = markdown_table(results)
    print("\n" + table)
    out = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "backend": backend.name, "results": results, "table_markdown": table}
    Path(args.out).write_text(json.dumps(out, indent=2), encoding="utf-8")
    Path(_HERE / "results_table.md").write_text(table, encoding="utf-8")
    print(f"\nSaved: {args.out}")
    print(f"Saved: {_HERE / 'results_table.md'}")


if __name__ == "__main__":
    main()
