"""
prompts.py — system prompt and per-mode prompt templates.

Both modes share one discipline: the model must answer **only** from the
provided context (the provenance record and/or retrieved literature passages),
cite the bracketed sources it used, and say "not in the provided context" rather
than guess. That discipline is what makes the output auditable and is what the
evaluation harness measures (grounding / hallucination rate).
"""

SYSTEM_PROMPT = (
    "You are the SOLANGE Copilot, a careful scientific assistant for a quantum-"
    "oncology research platform. You run entirely locally so that sensitive data "
    "never leaves the machine. Rules you must follow:\n"
    "1. Use ONLY the information in the CONTEXT provided. Do not add facts from "
    "outside it.\n"
    "2. When you use a numbered source like [1], cite it inline.\n"
    "3. If the answer is not in the context, say 'That is not in the provided "
    "context.' Never invent numbers, citations, or drug names.\n"
    "4. Be concise and precise. Prefer exact values (with units) from the context."
)


# Interpretation rules so the model explains what the numbers *mean* and checks
# internal consistency, instead of merely echoing field values.
_CONSISTENCY_GUIDE = (
    "When you describe the method, interpret the numbers — do not just repeat "
    "them — and check the record for internal consistency:\n"
    "- If the method is an exact or classical diagonalisation (ansatz is 'none' "
    "or similar), there is NO quantum circuit, so a gate count of 0 and a "
    "circuit depth of 0 are EXPECTED and consistent — say so explicitly rather "
    "than presenting 0 as surprising.\n"
    "- If the method uses a quantum circuit (e.g. VQE with an ansatz, or "
    "sample-based / sqDRIFT sampling), then gate count and circuit depth should "
    "be greater than 0; flag it as an inconsistency if they are 0.\n"
    "- If both a qubit count and an active space or encoding are given, note "
    "whether they agree (under Jordan-Wigner, qubits = 2 x number of spatial "
    "orbitals; e.g. CASSCF(12,12) uses 12 orbitals -> 24 qubits).\n"
    "- State whether the result carries a notarization seal (a hash) and who "
    "sealed it, and whether required fields for the stated phase are present.\n"
    "Conclude with a one-line verdict: consistent, or list the specific fields "
    "that conflict."
)


def explain_run_prompt(run_json: str, question: str | None = None) -> str:
    """Mode A — turn a P1-P9 provenance record into a plain-English narrative."""
    task = (question or
            "Explain this simulation run in plain English for a reviewer: what "
            "mutation was studied, what method and settings were used, what the "
            "resulting ground-state energy was, and how the result is notarized. "
            "Then state whether the record appears internally consistent.")
    return (
        "CONTEXT — SOLANGE provenance record (JSON):\n"
        f"{run_json}\n\n"
        f"TASK: {task}\n\n"
        f"{_CONSISTENCY_GUIDE}\n\n"
        "Answer using only the fields above. Quote exact numeric values with "
        "their units (Ha for energies). If a field is absent, say so."
    )


def druggability_prompt(mutation: str, context_block: str,
                        question: str | None = None) -> str:
    """Mode B — explain why a mutation is non-druggable, grounded in literature."""
    task = (question or
            f"Explain why the tumor-suppressor alteration '{mutation}' is "
            "considered non-druggable (or hard to drug), and summarize any "
            "indirect therapeutic strategies reported.")
    return (
        "CONTEXT — retrieved literature passages:\n"
        f"{context_block}\n\n"
        f"TASK: {task}\n\n"
        "Ground every claim in the numbered passages above and cite them like "
        "[1], [2]. If the passages do not address part of the question, say that "
        "explicitly instead of guessing."
    )


def chat_prompt(context_block: str, question: str) -> str:
    """Free-form grounded Q&A over the literature corpus."""
    return (
        "CONTEXT — retrieved passages:\n"
        f"{context_block}\n\n"
        f"QUESTION: {question}\n\n"
        "Answer only from the passages above and cite sources like [1]. If the "
        "answer is not present, say 'That is not in the provided context.'"
    )
