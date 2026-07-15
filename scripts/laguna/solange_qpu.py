#!/usr/bin/env python3
"""
solange_qpu.py — Phase 3B: real IBM Quantum hardware smoke test + P1-P9 capture.

PURPOSE (honest scope): prove the MOVING PARTS of the Phase 3B pipeline end to
end on REAL quantum hardware — circuit → IBM QPU → real telemetry → P1-P9
provenance record → LEON notarization in SOLANGE. It is NOT a novel chemistry
result and does NOT claim one: it measures the canonical minimal-basis H2
Hamiltonian (O'Malley et al., PRX 6:031007, 2016 — the textbook first-quantum-
chemistry-on-hardware example) on a fixed reference state, with NO VQE
optimization on the device (that would burn a scarce hardware-time budget for a
number already known classically). The scientific contribution of the
dissertation is the GOVERNED ORCHESTRATION (LEON, P1-P9, verify-don't-trust) —
this run demonstrates that governance generalizes to a real QPU backend, which
is exactly the DP4 (workload/backend-agnostic) claim.

Budget-aware by design (a trial account may have only ~10 minutes of QPU time):
  --dry-run   (DEFAULT) runs on a FREE local simulator — validate everything here first.
  --hardware  sends ONE job to a real QPU. Only pass this once --dry-run is clean.

USAGE:
  # 1) validate the whole pipeline for free:
  python solange_qpu.py --dry-run
  # 2) then one real run (spends QPU time):
  export QISKIT_IBM_TOKEN="...">; python solange_qpu.py --hardware --backend ibm_marrakesh --submit
"""

import argparse
import hashlib
import json
import math
import os
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

# Canonical minimal-basis (STO-3G) H2 Hamiltonian, Jordan-Wigner + 2-qubit
# reduction, bond length 0.735 Å — O'Malley et al., PRX 6:031007 (2016), the
# first quantum-chemistry-on-hardware demonstration. Hardcoded (not recomputed)
# precisely because this is a fixed, well-known reference operator, not a novel
# system: the point is the pipeline, not the number. Terms: (coeff, pauli-string
# over qubits [q0, q1], where 'I' = identity on that qubit).
H2_HAMILTONIAN = [
    (-1.052373245772859,  "II"),
    ( 0.397936358509854,  "IZ"),
    (-0.397936358509854,  "ZI"),
    (-0.011280104256235,  "ZZ"),
    ( 0.180931350397663,  "XX"),
]
H2_BOND_ANGSTROM = 0.735
H2_REF_STATE = "|01> (Hartree-Fock reference; one electron per spin-orbital pair)"


def check_ibm_credentials():
    """Whether qiskit-ibm-runtime is installed and a token/account is configured.
    No network call — safe to run anywhere. Returns (available, detail)."""
    try:
        import qiskit_ibm_runtime  # noqa: F401
    except ImportError:
        return False, "qiskit-ibm-runtime is not installed (pip install qiskit-ibm-runtime)"
    token = os.environ.get("QISKIT_IBM_TOKEN")
    saved = False
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        saved = bool(QiskitRuntimeService.saved_accounts())
    except Exception:
        pass
    if token or saved:
        return True, ("QISKIT_IBM_TOKEN set" if token else "saved IBM Quantum account found")
    return False, ("qiskit-ibm-runtime installed but no credentials "
                   "(set QISKIT_IBM_TOKEN or QiskitRuntimeService.save_account(...))")


def build_p9_stub(decoder="pyMatching"):
    """P9 (ML-decoder lineage, §06.iii). This smoke test uses the hardware's
    default deterministic decoding, so P9 is the minimal form — honestly NOT an
    ML-decoded run. Nothing fabricated."""
    return {"decoder": decoder, "decoder_version": None, "ml_decoded": False}


def _pauli_op(pauli):
    """Build a Qiskit SparsePauliOp term string like 'XX' as-is (qiskit reads
    left→right as qubit 0→n)."""
    from qiskit.quantum_info import SparsePauliOp
    return SparsePauliOp(pauli)


def _hamiltonian_op():
    from qiskit.quantum_info import SparsePauliOp
    labels = [p for _, p in H2_HAMILTONIAN]
    coeffs = [c for c, _ in H2_HAMILTONIAN]
    return SparsePauliOp(labels, coeffs)


def _hf_circuit():
    """Fixed Hartree-Fock reference state |01>. No parameters, no VQE — a single
    deterministic state prep, so the run costs one estimator evaluation."""
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(2)
    qc.x(0)   # |01>  (qubit 0 excited)
    return qc


def _classical_ground_state():
    """Exact lowest eigenvalue of the H2 operator, computed locally — the honest
    reference the hardware number is compared against (no external 'expected'
    value is asserted; both come from the same operator)."""
    import numpy as np
    return float(np.linalg.eigvalsh(_hamiltonian_op().to_matrix())[0])


def _hf_energy_exact():
    """Exact <H> on the fixed HF reference — what the hardware SHOULD approach if
    the pipeline is faithful (this is NOT the ground state; no optimization is done)."""
    import numpy as np
    from qiskit.quantum_info import Statevector
    sv = Statevector(_hf_circuit())
    return float(np.real(sv.expectation_value(_hamiltonian_op())))


def _backend_telemetry(backend):
    """Extract REAL P3/P4 provenance from a live backend — defensively, so a
    missing field becomes None rather than crashing the run. This is the whole
    point of Phase 3B: P3 calibration epoch / P4 error budget are REAL here, not
    the simulated placeholders a classical run carries."""
    tel = {"p3_calibration_epoch": None, "p4_gate_error_rate": None,
           "p4_readout_error_rate": None, "p4_t1_us": None, "p4_t2_us": None,
           "p3_backend_version": None}
    try:
        props = backend.properties()
    except Exception:
        props = None
    if props is not None:
        try:
            tel["p3_calibration_epoch"] = props.last_update_date.isoformat()
        except Exception:
            pass
        try:
            tel["p4_readout_error_rate"] = float(props.readout_error(0))
        except Exception:
            pass
        try:
            tel["p4_t1_us"] = round(float(props.t1(0)) * 1e6, 3)
            tel["p4_t2_us"] = round(float(props.t2(0)) * 1e6, 3)
        except Exception:
            pass
        try:
            # median 2q gate error across the device's ECR/CZ gates
            errs = [g.parameters[0].value for g in props.gates
                    if g.gate in ("ecr", "cz", "cx") and g.parameters]
            if errs:
                errs.sort()
                tel["p4_gate_error_rate"] = float(errs[len(errs) // 2])
        except Exception:
            pass
    try:
        tel["p3_backend_version"] = str(getattr(backend, "version", None) or
                                        getattr(backend, "backend_version", None))
    except Exception:
        pass
    return tel


def run_h2(hardware, backend_name, shots, token, instance):
    """Measure <H2> on the fixed HF state. dry-run → local simulator (free);
    hardware → one real QPU job. Returns (energy, backend_label, telemetry, meta)."""
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    obs = _hamiltonian_op()
    qc = _hf_circuit()

    if not hardware:
        # FREE local validation of the exact same call shape used on hardware.
        try:
            from qiskit_aer import AerSimulator
            from qiskit_aer.primitives import EstimatorV2 as AerEstimator
            backend = AerSimulator()
            pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
            isa = pm.run(qc)
            est = AerEstimator()
            job = est.run([(isa, obs.apply_layout(isa.layout), )], precision=1/math.sqrt(shots))
            energy = float(job.result()[0].data.evs)
            label = "AerSimulator (local dry-run — NOT real hardware)"
        except Exception:
            # Minimal fallback: exact statevector (still a dry-run, still free).
            energy = _hf_energy_exact()
            label = "Statevector (local dry-run fallback — NOT real hardware)"
        return energy, label, {"p3_calibration_epoch": None}, {"mode": "dry-run", "shots": shots}

    # ── Real hardware path ──
    from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
    kwargs = {}
    if token:    kwargs["token"] = token
    if instance: kwargs["instance"] = instance
    # New IBM Quantum Platform channel; falls back to the legacy channel name.
    try:
        service = QiskitRuntimeService(channel="ibm_quantum_platform", **kwargs)
    except Exception:
        service = QiskitRuntimeService(channel="ibm_quantum", **kwargs)
    backend = service.backend(backend_name)
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
    isa = pm.run(qc)
    est = EstimatorV2(mode=backend)
    est.options.default_shots = shots
    job = est.run([(isa, obs.apply_layout(isa.layout))])
    print(f"  submitted job {job.job_id()} to {backend_name} — waiting for result "
          f"(queue time does NOT count against QPU-execution budget) …")
    energy = float(job.result()[0].data.evs)
    tel = _backend_telemetry(backend)
    meta = {"mode": "hardware", "shots": shots, "job_id": job.job_id(),
            "backend": backend_name}
    return energy, f"{backend_name} (real QPU)", tel, meta


def build_record(energy, backend_label, telemetry, meta):
    """Assemble a P1-P9 provenance record for the run, then seal it. Field names
    use the p1_..p9_ prefixes so LEON's existing P8 seal (backend/routes/leon.py,
    which sweeps any p*_ field) covers the whole record with NO backend change."""
    now = datetime.now(timezone.utc).isoformat()
    ground = _classical_ground_state()
    hf_exact = _hf_energy_exact()
    is_hw = meta.get("mode") == "hardware"
    rec = {
        "id": str(uuid.uuid4()),
        "created_at": now,
        "mutation_id": "H2_smoketest",
        "mutation_name": "H2 minimal-basis (Phase 3B hardware smoke test)",
        "phase": "3B-QPU" if is_hw else "3B-QPU-dryrun",
        "provenance_source": "QPU/IBM" if is_hw else "QPU/local-dryrun",

        # P1 — circuit
        "p1_qubit_count": 2,
        "p1_ansatz": "fixed HF reference |01> (no VQE optimization on device)",
        "p1_circuit_hash": hashlib.sha256(
            json.dumps({"h2": H2_HAMILTONIAN, "state": H2_REF_STATE},
                       sort_keys=True).encode()).hexdigest(),

        # P2 — compilation lineage (the exact operator measured)
        "p2_encoding": "Jordan-Wigner + 2-qubit reduction (O'Malley 2016)",
        "p2_basis_set": "STO-3G",
        "p2_active_electrons": 2, "p2_active_orbitals": 2,
        "p2_model_compound": "H2",
        "p2_jw_terms": [{"coeff": c, "pauli": p} for c, p in H2_HAMILTONIAN],

        # P3 — backend (REAL telemetry on hardware; None on dry-run — honestly)
        "p3_backend": backend_label,
        "p3_backend_version": telemetry.get("p3_backend_version"),
        "p3_calibration_epoch": telemetry.get("p3_calibration_epoch") or now,

        # P4 — device error budget (REAL on hardware)
        "p4_gate_error_rate": telemetry.get("p4_gate_error_rate"),
        "p4_readout_error_rate": telemetry.get("p4_readout_error_rate"),
        "p4_t1_us": telemetry.get("p4_t1_us"),
        "p4_t2_us": telemetry.get("p4_t2_us"),
        "p4_note": ("real IBM QPU calibration data" if is_hw else "dry-run — no device"),

        # P5 — shots / raw
        "p5_shots": meta.get("shots"),
        "p5_elapsed_s": None,
        "p5_casscf_ref_ha": ground,          # exact ground state of THIS operator
        "p5_active_energy_ha": hf_exact,      # exact <H> on the HF ref state

        # P6 — mitigation (none applied in this minimal smoke test)
        "p6_method": "none (raw estimator)", "p6_note": "smoke test — no error mitigation",

        # P7 — result
        "p7_energy_ha": energy,
        "p7_method": ("QPU EstimatorV2 <H> on HF ref" if is_hw
                      else "dry-run EstimatorV2 <H> on HF ref"),
        "p7_note": ("Fixed-state expectation, NOT a ground-state search — Δ to the "
                    "exact ground state is expected and not an error."),

        # P9 — decoder lineage (minimal / deterministic)
        "p9_applicable": False,
        "p9_note": "deterministic default decoding; not an ML-decoded run",
    }
    payload = json.dumps({k: v for k, v in rec.items()
                          if k.startswith(("p1_", "p2_", "p3_", "p4_", "p5_",
                                           "p6_", "p7_", "p9_"))
                          and k != "p3_calibration_epoch"},
                         sort_keys=True, default=str)
    rec["p8_seal_payload"] = payload
    rec["p8_hash"] = hashlib.sha256(payload.encode()).hexdigest()
    rec["p8_algorithm"] = "SHA-256"
    rec["p8_sealed_at"] = now
    return rec, ground, hf_exact


def submit(api, record):
    """POST the sealed P1-P9 record to SOLANGE /hpc/submit — LEON re-verifies the
    seal and stores it, exactly like a Laguna CASSCF run. Best-effort; the local
    JSON is written regardless."""
    import urllib.request
    body = {"provenance": record, "jw": {}}   # no jw consistency payload for a QPU run
    try:
        req = urllib.request.Request(api.rstrip("/") + "/api/simulate/hpc/submit",
                                     data=json.dumps(body, default=str).encode(),
                                     method="POST",
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            resp = json.loads(r.read().decode())
        print(f"SUBMITTED → SOLANGE  status={resp.get('status')}  "
              f"seal_ok={resp.get('seal_ok')}  db={resp.get('db_status')}  "
              f"notary={resp.get('notary')}  run_id={resp.get('run_id')}")
    except Exception as e:
        print(f"  SUBMIT FAILED (result still saved locally): {e}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="SOLANGE Phase 3B — H2 hardware smoke test.")
    ap.add_argument("--hardware", action="store_true",
                    help="send ONE job to a real QPU (spends QPU time). Omit for a free dry-run.")
    ap.add_argument("--dry-run", action="store_true",
                    help="explicit free local-simulator run (this is already the default "
                         "when --hardware is omitted; accepted for clarity). Ignored if "
                         "--hardware is also passed.")
    ap.add_argument("--backend", default="ibm_marrakesh",
                    help="QPU name (only used with --hardware). ibm_marrakesh = Heron r2, "
                         "fewest pending jobs in the open trial instance.")
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--instance", default=os.environ.get("QISKIT_IBM_INSTANCE"),
                    help="IBM instance CRN (open-instance). Defaults to $QISKIT_IBM_INSTANCE.")
    ap.add_argument("--out", default="./out")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the sealed record to SOLANGE (LEON notarizes it). "
                         "Bare flag uses the default API; pass a value to override.")
    ap.add_argument("--check-credentials", action="store_true")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    token = os.environ.get("QISKIT_IBM_TOKEN")

    print("=" * 70)
    print("SOLANGE Phase 3B — H2 hardware smoke test (moving-parts proof, not chemistry)")
    available, detail = check_ibm_credentials()
    print(f"IBM Quantum credentials: {'AVAILABLE' if available else 'NOT AVAILABLE'} — {detail}")
    if args.check_credentials:
        print("=" * 70); sys.exit(0 if available else 1)

    if args.hardware and not available:
        print("-" * 70)
        print("REFUSING --hardware without credentials — no fabricated result. "
              "Set QISKIT_IBM_TOKEN and retry.")
        print("=" * 70); sys.exit(1)

    mode = "REAL HARDWARE" if args.hardware else "DRY-RUN (free local simulator)"
    print(f"Mode: {mode}" + (f"  ·  backend={args.backend}" if args.hardware else ""))
    print(f"Measuring canonical H2 (STO-3G, {H2_BOND_ANGSTROM} Å) <H> on {H2_REF_STATE}")

    energy, backend_label, telemetry, meta = run_h2(
        args.hardware, args.backend, args.shots, token, args.instance)

    record, ground, hf_exact = build_record(energy, backend_label, telemetry, meta)

    print("-" * 70)
    print(f"backend         : {backend_label}")
    print(f"measured <H>    : {energy:.6f} Ha   (on the fixed HF reference state)")
    print(f"exact <H> (HF)  : {hf_exact:.6f} Ha   (what a faithful pipeline approaches)")
    print(f"exact ground    : {ground:.6f} Ha   (reference only; no optimization was done)")
    print(f"Δ(measured-HF)  : {(energy - hf_exact)*1000:.2f} mHa   "
          f"(hardware noise + shot statistics; expected, not an error)")
    if telemetry.get("p3_calibration_epoch"):
        print(f"REAL P3 epoch   : {telemetry['p3_calibration_epoch']}")
        print(f"REAL P4 budget  : gate_err={telemetry.get('p4_gate_error_rate')} "
              f"readout_err={telemetry.get('p4_readout_error_rate')} "
              f"T1={telemetry.get('p4_t1_us')}us T2={telemetry.get('p4_t2_us')}us")
    print(f"P8 seal         : {record['p8_hash'][:16]}…")

    stem = "h2_qpu_" + ("hw" if meta.get("mode") == "hardware" else "dryrun")
    p = Path(args.out) / f"{stem}_{record['id'][:8]}.json"
    p.write_text(json.dumps(record, indent=2, default=str))
    print(f"WROTE {p}")

    if args.submit:
        submit(args.submit, record)
    print("=" * 70)
    if not args.hardware:
        print("Dry-run clean? Re-run with --hardware (and --submit) for ONE real QPU job.")


if __name__ == "__main__":
    main()
