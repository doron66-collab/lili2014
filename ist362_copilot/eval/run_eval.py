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


def eval_grounding(backend, model, dataset, temperature=0.2):
    """Mode B / chat: questions answerable from the corpus."""
    retriever = get_retriever()
    rows, correct, lat = [], 0, []
    for item in dataset["grounding"]["items"]:
        hits = retriever.search(item["question"], top_k=4)
        prompt = prompts.chat_prompt(retriever.context_block(hits), item["question"])
        res = backend.generate(prompt, model=model, system=prompts.SYSTEM_PROMPT,
                               temperature=temperature)
        ok = contains_any(res.text, item["expect_any"])
        correct += ok
        lat.append(res.latency_s)
        rows.append({"q": item["question"][:60], "ok": ok})
    n = len(dataset["grounding"]["items"])
    return {"accuracy": correct / n, "n": n, "mean_latency_s": _mean(lat), "rows": rows}


def eval_traps(backend, model, dataset, temperature=0.2):
    """Out-of-context questions: correct behavior is to refuse."""
    retriever = get_retriever()
    markers = dataset["traps"]["refusal_markers"]
    rows, refused, lat = [], 0, []
    for item in dataset["traps"]["items"]:
        hits = retriever.search(item["question"], top_k=4)
        prompt = prompts.chat_prompt(retriever.context_block(hits), item["question"])
        res = backend.generate(prompt, model=model, system=prompts.SYSTEM_PROMPT,
                               temperature=temperature)
        ok = contains_any(res.text, markers)
        refused += ok
        lat.append(res.latency_s)
        rows.append({"q": item["question"][:60], "refused": ok, "answer": res.text[:100]})
    n = len(dataset["traps"]["items"])
    return {"refusal_rate": refused / n, "n": n, "mean_latency_s": _mean(lat), "rows": rows}


def _mean(xs):
    return round(sum(xs) / len(xs), 3) if xs else 0.0


_CITE_RE = __import__("re").compile(r"\[\d+\]")


def eval_mutation(backend, model, mutation, probes):
    """Mode B for ONE mutation: how well does each model explain *this* mutation?"""
    retriever = get_retriever()
    hits = retriever.search(mutation + " undruggable non-druggable therapy", top_k=4)
    prompt = prompts.druggability_prompt(mutation, retriever.context_block(hits))
    res = backend.generate(prompt, model=model, system=prompts.SYSTEM_PROMPT)
    covered = [p["fact"] for p in probes if contains_any(res.text, p["expect_any"])]
    facts_score = len(covered) / len(probes)
    cited = bool(_CITE_RE.search(res.text))
    return {
        "model": model or "(default)",
        "backend": backend.name,
        "facts_covered": round(facts_score, 3),
        "facts_detail": f"{len(covered)}/{len(probes)}",
        "cited_sources": cited,
        "latency_s": round(res.latency_s, 3),
        "answer": res.text,
    }


_MUT_COLS = [
    ("Model", lambda r: r["model"]),
    ("Facts covered", lambda r: f"{r['facts_covered']:.2f} ({r['facts_detail']})"),
    ("Cited sources", lambda r: "yes" if r["cited_sources"] else "no"),
    ("Latency (s)", lambda r: f"{r['latency_s']:.2f}"),
]


def _aligned(cols, results):
    headers = [c[0] for c in cols]
    rows = [[get(r) for _, get in cols] for r in results]
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) if rows
              else len(headers[i]) for i in range(len(headers))]

    def line(l, m, rr):
        return l + m.join("─" * (w + 2) for w in widths) + rr

    def fmt(cells):
        return "│ " + " │ ".join(c.center(widths[i]) for i, c in enumerate(cells)) + " │"

    out = [line("┌", "┬", "┐"), fmt(headers), line("├", "┼", "┤")]
    out += [fmt(row) for row in rows]
    out.append(line("└", "┴", "┘"))
    return "\n".join(out)


_TEMP_COLS = [
    ("Temperature", lambda r: f"{r['temperature']:.1f}"),
    ("Grounding acc.", lambda r: f"{r['grounding_accuracy']:.2f}"),
    ("Hallucination refusal", lambda r: f"{r['hallucination_refusal_rate']:.2f}"),
    ("Mean latency (s)", lambda r: f"{r['mean_latency_s']:.2f}"),
]


def run_temperature_eval(backend, model, dataset, temps):
    """
    Randomness study (Topic 1): sampling temperature IS the randomness knob in a
    language model. Higher temperature = more random token choices. We measure how
    that randomness affects faithfulness — grounding accuracy on answerable
    questions, and refusal rate on out-of-context 'trap' questions.
    """
    model = model or None
    print(f"\n=== Randomness (temperature) study — model '{model or '(default)'}' "
          f"on backend '{backend.name}' ===")
    print("    (higher temperature = more randomness in the model's sampling)")
    results = []
    for t in temps:
        grd = eval_grounding(backend, model, dataset, temperature=t)
        trp = eval_traps(backend, model, dataset, temperature=t)
        results.append({
            "temperature": t,
            "grounding_accuracy": round(grd["accuracy"], 3),
            "hallucination_refusal_rate": round(trp["refusal_rate"], 3),
            "mean_latency_s": round((grd["mean_latency_s"] + trp["mean_latency_s"]) / 2, 3),
        })
    table = _aligned(_TEMP_COLS, results)
    print("\n" + table + "\n")
    out = {"model": model or "(default)", "backend": backend.name,
           "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "results": results}
    Path(_HERE / "results_temperature.json").write_text(json.dumps(out, indent=2), "utf-8")
    Path(_HERE / "results_temperature.txt").write_text(table, "utf-8")
    print(f"Saved: {_HERE / 'results_temperature.json'}")
    print(f"Saved: {_HERE / 'results_temperature.txt'}")
    if backend.name == "mock":
        print("\nNote: the mock backend is deterministic, so temperature has no "
              "effect here. Run with SOLANGE_BACKEND=ollama for real variation.")


def run_mutation_eval(backend, models, mutation, dataset):
    probes = dataset["mutation_probes"].get(mutation)
    if probes is None:
        avail = ", ".join(dataset["mutation_probes"].keys() - {"description"})
        raise SystemExit(f"No probes defined for '{mutation}'. Available: {avail}")
    print(f"\n=== Druggability (Mode B) evaluation for: {mutation} "
          f"on backend '{backend.name}' ===")
    results = [eval_mutation(backend, m, mutation, probes) for m in models]
    table = _aligned(_MUT_COLS, results)
    print("\n" + table + "\n")
    out = {"mutation": mutation, "backend": backend.name,
           "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "results": results}
    slug = mutation.replace(" ", "_").replace("/", "-")
    Path(_HERE / f"results_{slug}.json").write_text(json.dumps(out, indent=2), "utf-8")
    Path(_HERE / f"results_{slug}.txt").write_text(table, "utf-8")
    print(f"Saved: {_HERE / ('results_' + slug + '.json')}")
    print(f"Saved: {_HERE / ('results_' + slug + '.txt')}")


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


_COLS = [
    ("Model", lambda r: r["model"]),
    ("Backend", lambda r: r["backend"]),
    ("Extraction acc.", lambda r: f"{r['extraction_accuracy']:.2f}"),
    ("Grounding acc.", lambda r: f"{r['grounding_accuracy']:.2f}"),
    ("Hallucination refusal", lambda r: f"{r['hallucination_refusal_rate']:.2f}"),
    ("Mean latency (s)", lambda r: f"{r['mean_latency_s']:.2f}"),
]


def markdown_table(results: list[dict]) -> str:
    hdr = "| " + " | ".join(c[0] for c in _COLS) + " |\n"
    sep = "|" + "|".join("---" for _ in _COLS) + "|\n"
    rows = "".join("| " + " | ".join(get(r) for _, get in _COLS) + " |\n"
                   for r in results)
    return hdr + sep + rows


def aligned_table(results: list[dict]) -> str:
    """Fixed-width, box-drawn table — readable directly in the terminal."""
    headers = [c[0] for c in _COLS]
    rows = [[get(r) for _, get in _COLS] for r in results]
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) if rows
              else len(headers[i]) for i in range(len(headers))]

    def line(left, mid, right):
        return left + mid.join("─" * (w + 2) for w in widths) + right

    def fmt(cells):
        return "│ " + " │ ".join(c.center(widths[i]) for i, c in enumerate(cells)) + " │"

    out = [line("┌", "┬", "┐"), fmt(headers), line("├", "┼", "┤")]
    out += [fmt(row) for row in rows]
    out.append(line("└", "┴", "┘"))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description="Evaluate the SOLANGE Copilot.")
    ap.add_argument("--models", nargs="*", default=[None],
                    help="model names to compare (default: backend default)")
    ap.add_argument("--mutation", default=None,
                    help="run the per-mutation Druggability eval for this mutation "
                         "(e.g. 'TP53 C275F') instead of the full benchmark")
    ap.add_argument("--temperature-sweep", action="store_true",
                    help="run the randomness study: sweep sampling temperature "
                         "and measure its effect on grounding & hallucination")
    ap.add_argument("--temps", nargs="*", type=float, default=[0.0, 0.5, 1.0],
                    help="temperatures for --temperature-sweep (default: 0.0 0.5 1.0)")
    ap.add_argument("--out", default=str(_HERE / "results.json"))
    args = ap.parse_args()

    dataset = json.loads((_HERE / "eval_dataset.json").read_text("utf-8"))
    run_record = json.loads((_DATA / "sample_run.json").read_text("utf-8"))
    backend = get_backend()

    # Randomness study: one model, several temperatures.
    if args.temperature_sweep:
        run_temperature_eval(backend, args.models[0], dataset, args.temps)
        return

    # Mutation-focused mode: one mutation, all models, Mode B only.
    if args.mutation:
        run_mutation_eval(backend, args.models, args.mutation, dataset)
        return

    results = [run_for_model(backend, m, dataset, run_record) for m in args.models]

    table_md = markdown_table(results)
    table_txt = aligned_table(results)
    print("\n" + table_txt + "\n")            # readable in the terminal
    out = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "backend": backend.name, "results": results,
           "table_markdown": table_md, "table_aligned": table_txt}
    Path(args.out).write_text(json.dumps(out, indent=2), encoding="utf-8")
    Path(_HERE / "results_table.md").write_text(table_md, encoding="utf-8")
    Path(_HERE / "results_table.txt").write_text(table_txt, encoding="utf-8")
    print(f"Saved: {args.out}")
    print(f"Saved: {_HERE / 'results_table.md'}  (paste into the paper)")
    print(f"Saved: {_HERE / 'results_table.txt'}  (aligned, for the terminal)")


if __name__ == "__main__":
    main()
