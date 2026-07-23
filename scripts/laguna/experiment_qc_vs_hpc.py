#!/usr/bin/env python3
"""
experiment_qc_vs_hpc.py — head-to-head: the SAME active-space Hamiltonian evaluated
CLASSICALLY (exact, on this Laguna node) and on REAL QUANTUM HARDWARE, then compared.

WHY THIS IS A CLEAN EXPERIMENT
Classical and quantum evaluate the *identical* quantity — ⟨H⟩ on the *identical*
fixed Hartree–Fock state — from the *identical* Jordan–Wigner Hamiltonian. So their
difference is **pure hardware noise**: no optimization variance, no ansatz mismatch,
no different active space. That is exactly what makes it a fair race.

Run it ON Laguna and it does both sides in one command: the exact classical
diagonalization runs here (CPU/largemem — trivial at 4 qubits, real work as the
active space grows), and the quantum measurement is dispatched to IBM Heron.

THREE NUMBERS (same H):
  E_ground   exact classical ground state          → the truth
  E_HF       exact classical ⟨H⟩ on the HF state    → classical reference
  E_QPU      measured ⟨H⟩ on the HF state, on IBM   → real hardware
TWO HONEST DELTAS:
  hardware noise = |E_QPU − E_HF|      (the cost of real hardware, isolated)
  HF−vs−ground   = |E_HF − E_ground|   (physics: what optimization must still close)

Order: --dry-run first (free, simulator → noise≈0 confirms the pipeline), then
--hardware (spends real quantum time → the real noise number).
"""
import argparse
import json
import os
import sys
from pathlib import Path

# Reuse the validated building blocks — one source of truth, no duplicated physics.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from solange_qpu import (jw_target, h2_target, _exact, measure, build_record,
                         submit, check_ibm_credentials)


def main():
    ap = argparse.ArgumentParser(description="QC vs classical head-to-head on the same Hamiltonian.")
    ap.add_argument("--key", help="SOLANGE target key (omit → H2 benchmark)")
    ap.add_argument("--side", default="native", choices=["native", "mutant"])
    ap.add_argument("--jw-file", default=str(Path(__file__).resolve().parents[2] / "jw_hamiltonians.json"))
    ap.add_argument("--hardware", action="store_true", help="run the quantum side on REAL IBM hardware")
    ap.add_argument("--dry-run", action="store_true", help="run the quantum side on a local simulator (free)")
    ap.add_argument("--backend", default="ibm_kingston", help="IBM backend (any online: ibm_fez/ibm_marrakesh/ibm_kingston)")
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--instance", default=os.environ.get("QISKIT_IBM_INSTANCE"))
    ap.add_argument("--out", default="./out")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the QPU row to SOLANGE (LEON notarizes it) so it lands in Rung 4")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    token = os.environ.get("QISKIT_IBM_TOKEN")

    print("=" * 76)
    print("SOLANGE experiment — same Hamiltonian: CLASSICAL (exact) vs REAL QUANTUM")
    available, detail = check_ibm_credentials()
    print(f"IBM credentials: {'AVAILABLE' if available else 'NOT AVAILABLE'} — {detail}")
    if args.hardware and not available:
        print("REFUSING hardware without credentials — no fabricated result. "
              "Set QISKIT_IBM_TOKEN (+ instance) and retry.")
        print("=" * 76); sys.exit(1)

    target = jw_target(args.key, args.side, args.jw_file) if args.key else h2_target()
    print(f"Target: {target['label']}  ·  {target['nq']} qubits")
    print(f"Quantum side: {'REAL HARDWARE · '+args.backend if args.hardware else 'DRY-RUN (local simulator)'}")
    print("-" * 76)

    # ── CLASSICAL side (exact diagonalization — runs here, on this Laguna node) ──
    E_HF, E_ground = _exact(target["obs"], target["circuit"])
    print(f"[classical · exact]  E_ground = {E_ground:.6f} Ha   (the truth)")
    print(f"[classical · exact]  E_HF     = {E_HF:.6f} Ha   (⟨H⟩ on the fixed HF state)")

    # ── QUANTUM side (same H, same HF state) ──
    E_QPU, backend_label, telemetry, meta = measure(target, args.hardware, args.backend,
                                                     args.shots, token, args.instance)
    print(f"[quantum  · {('hardware' if args.hardware else 'sim').ljust(8)}] E_QPU    = {E_QPU:.6f} Ha   (measured ⟨H⟩ on {backend_label})")

    noise_mha = abs(E_QPU - E_HF) * 1000.0
    gap_mha   = abs(E_HF - E_ground) * 1000.0
    print("-" * 76)
    print("RESULT")
    print(f"  hardware noise = |E_QPU − E_HF|    = {noise_mha:8.3f} mHa   "
          f"({'expected ≈0 on a simulator' if not args.hardware else 'the real cost of hardware'})")
    print(f"  HF − ground    = |E_HF − E_ground| = {gap_mha:8.3f} mHa   "
          f"(physics: what a full VQE/SQD would still need to close)")
    print("-" * 76)

    # ── benchmark artifact ──
    bench = {
        "experiment": "qc_vs_hpc_same_hamiltonian",
        "target": target["label"], "key": args.key, "side": args.side, "qubits": target["nq"],
        "classical_backend": "exact diagonalization (Laguna CPU)",
        "quantum_backend": backend_label,
        "E_ground_ha": E_ground, "E_HF_ha": E_HF, "E_QPU_ha": E_QPU,
        "hardware_noise_mha": noise_mha, "hf_vs_ground_mha": gap_mha,
        "mode": meta.get("mode"), "shots": args.shots,
        "telemetry": telemetry,
    }
    stem = f"experiment_{(args.key or 'H2')}_{args.side}_" + ("hw" if args.hardware else "sim")
    bpath = Path(args.out) / f"{stem}.json"
    bpath.write_text(json.dumps(bench, indent=2, default=str))
    print(f"WROTE benchmark → {bpath}")

    # ── optionally push the QPU row to SOLANGE (Rung 4), sealed by LEON ──
    if args.submit:
        record = build_record(target, E_QPU, E_HF, backend_label, telemetry, meta)
        submit(args.submit, record)

    print("=" * 76)
    if not args.hardware:
        print("Dry-run clean (noise≈0)? Re-run with --hardware for the real quantum-vs-classical number.")


if __name__ == "__main__":
    main()
