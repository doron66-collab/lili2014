#!/usr/bin/env python3
"""
solange_qpu.py — Phase 3B scaffold: IBM Heron r3 / sqDRIFT execution interface.

HONEST STATUS: this is a SCAFFOLD, not a working quantum execution path. It
defines and tests the parts that don't require real hardware access — the CLI
shape (mirroring solange_hpc.py / solange_dmrg.py), the credential check, and
the P9 ML-decoder provenance record structure specified in the dissertation
(§06.iii) — and seals that record the same way LEON already seals a DMRG
classification (backend/routes/leon.py :: notarize_generic — no backend changes
needed when Phase 3B goes live; the generic seal path already supports any
record shape).

What this script deliberately does NOT do: fabricate a quantum result. Building
the actual sqDRIFT Trotter circuit from a Jordan-Wigner Hamiltonian (the qDRIFT
sampling procedure over the JW terms solange_hpc.py already builds) and
executing it via qiskit-ibm-runtime is real algorithmic + hardware-access work
that cannot be honestly written, let alone claimed correct, without IBM Quantum
credentials to test against. run_sqdrift_circuit() raises NotImplementedError
with a precise description of what remains — a documented gap, not a silent one.

USAGE (credential check only — this is all that can be verified right now):
  python solange_qpu.py --key TP53_C275F --side mutant --check-credentials
"""

import argparse
import hashlib
import json
import os
import sys
import time
import uuid
from pathlib import Path

DEFAULT_BACKEND = "ibm_heron_r3"


def check_ibm_credentials():
    """Report whether qiskit-ibm-runtime is installed and an account/token is
    configured — WITHOUT attempting to connect (no network call here). Returns
    (available: bool, detail: str)."""
    try:
        import qiskit_ibm_runtime  # noqa: F401
    except ImportError:
        return False, "qiskit-ibm-runtime is not installed (pip install qiskit-ibm-runtime)"

    token = os.environ.get("QISKIT_IBM_TOKEN")
    saved_account = False
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        # Checking for a saved account file does not require network access.
        saved_account = bool(QiskitRuntimeService.saved_accounts())
    except Exception:
        pass

    if token or saved_account:
        return True, ("QISKIT_IBM_TOKEN set" if token else "saved IBM Quantum account found")
    return False, ("qiskit-ibm-runtime is installed, but no credentials found "
                    "(set QISKIT_IBM_TOKEN or run QiskitRuntimeService.save_account(...))")


def build_p9_stub(decoder="pyMatching", ml_decoder=False):
    """The nine-element provenance record's P9 field, per §06.iii of the
    dissertation. Deterministic decoding (pyMatching, BP-OSD) degrades to the
    minimal form; an ML-decoded QEC path (e.g. an Ising-class pre-decoder) would
    populate the additional fields. Nothing here is fabricated — fields this
    script cannot know yet (weights checksum, training-corpus hash) are None,
    not guessed.
    """
    if not ml_decoder:
        return {"decoder": decoder, "decoder_version": None, "ml_decoded": False}
    return {
        "decoder": decoder, "decoder_version": None, "ml_decoded": True,
        "pre_decoder_id": None, "pre_decoder_weights_checksum": None,
        "training_corpus_manifest_hash": None,
        "inference_config": {"precision": None, "batch_size": None, "threshold": None},
        "execution_environment": {"gpu_class": None, "cuda_runtime_version": None},
        "global_decoder": {"algorithm": decoder, "version": None},
    }


# ── LEON seal (self-contained, mirrors solange_dmrg.py / solange_hpc.py — this
# script must run standalone with no dependency on the backend package). When
# Phase 3B circuit execution is implemented, submitting a real result needs NO
# new LEON logic: leon.notarize_generic() already accepts any record shape.
def _seal_payload(record, exclude):
    return json.dumps({k: v for k, v in record.items() if k not in exclude},
                      sort_keys=True, default=str)


def _seal_hash(record, exclude):
    return hashlib.sha256(_seal_payload(record, exclude).encode()).hexdigest()


def run_sqdrift_circuit(jw_terms, backend, shots, decoder):
    """Execute the sqDRIFT sampling procedure on IBM Heron and return a
    provenance-complete result. NOT YET IMPLEMENTED.

    What remains, precisely:
      1. qDRIFT Trotter step sampling over `jw_terms` (weighted-random term
         selection per the qDRIFT algorithm — see arXiv:2508.02578) to build the
         sample circuits, reusing the JW term list solange_hpc.py already
         produces (build_jw_terms) as input.
      2. Submission via qiskit_ibm_runtime.SamplerV2 on a Session bound to
         `backend`, with `shots` per circuit.
      3. Classical diagonalization of the sampled configurations (the "SQD" half
         of sample-based quantum diagonalization) to recover the ground-state
         energy estimate.
      4. Populating P1 (circuit fingerprint), P3 (calibration epoch from the
         backend's live properties), P4 (per-gate error budget), P5 (raw shot
         histogram), P6 (error mitigation applied), P7 (estimator + 95% CI).
    Each step requires live IBM Quantum access to implement AND validate — it
    cannot be written honestly offline. This function exists so the call site
    (main()) and provenance shape are already correct; only this body remains.
    """
    raise NotImplementedError(
        "Phase 3B sqDRIFT circuit execution is not yet implemented — see this "
        "function's docstring for the precise remaining steps. This is a "
        "documented scope gap (§08 Limitations), not a silent one."
    )


def main():
    ap = argparse.ArgumentParser(
        description="SOLANGE Phase 3B scaffold — IBM Heron r3 / sqDRIFT (NOT YET EXECUTABLE).")
    ap.add_argument("--key", help="target key, e.g. TP53_C275F")
    ap.add_argument("--side", default="mutant", choices=["native", "mutant"])
    ap.add_argument("--backend", default=DEFAULT_BACKEND)
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--decoder", default="pyMatching",
                    help="deterministic global decoder (pyMatching, BP-OSD, ...)")
    ap.add_argument("--ml-decoder", action="store_true",
                    help="mark this run as using a learned pre-decoder (Ising-class) — "
                         "populates the extended P9 fields (still None until Phase 3B "
                         "circuit execution is implemented)")
    ap.add_argument("--out", default="./out")
    ap.add_argument("--check-credentials", action="store_true",
                    help="only check for qiskit-ibm-runtime + IBM Quantum credentials, "
                         "then exit — the one thing this script can honestly verify today")
    args = ap.parse_args()

    print("=" * 68)
    print("SOLANGE Phase 3B scaffold — IBM Heron r3 / sqDRIFT")
    available, detail = check_ibm_credentials()
    print(f"IBM Quantum credentials: {'AVAILABLE' if available else 'NOT AVAILABLE'} — {detail}")

    if args.check_credentials:
        print("=" * 68)
        sys.exit(0 if available else 1)

    if not args.key:
        ap.error("--key is required (unless --check-credentials)")

    p9 = build_p9_stub(decoder=args.decoder, ml_decoder=args.ml_decoder)
    print(f"P9 provenance shape (ml_decoded={p9['ml_decoded']}): {json.dumps(p9)}")

    if not available:
        print("-" * 68)
        print("STOPPING — refusing to fabricate a result without IBM Quantum access.")
        print("This is the honest behavior: no local JSON is written, no result is")
        print("claimed. Configure credentials (see message above) and re-run.")
        print("=" * 68)
        sys.exit(1)

    # Credentials ARE available — but circuit synthesis itself is still a gap.
    try:
        run_sqdrift_circuit(jw_terms=None, backend=args.backend, shots=args.shots,
                            decoder=args.decoder)
    except NotImplementedError as e:
        print("-" * 68)
        print(f"NOT YET IMPLEMENTED: {e}")
        print("=" * 68)
        sys.exit(1)


if __name__ == "__main__":
    main()
