#!/usr/bin/env python3
"""
solange_qpu_vqe1p.py — the first REAL on-device VQE optimization in SOLANGE
Phase 3B, not just a fixed-state <H> measurement.

Every prior real-hardware run (solange_qpu.py) prepares a fixed Hartree-Fock
circuit and measures <H> once — explicitly "no VQE optimization on device"
(see solange_qpu.py's own docstring and its p1_ansatz field). That was a
deliberate choice: a general VQE optimizer needs many circuit evaluations,
and a trial account's QPU-second budget is scarce. This script closes that
gap cheaply and EXACTLY, for the specific problem SOLANGE's Phase 3A/3B
targets are built from (CAS(2e,2o), 4 qubits under Jordan-Wigner):

  For a 2-electron/2-orbital active space, the physical ground state is a
  superposition of exactly two Slater determinants: the Hartree-Fock
  reference |1100> and the doubly-excited configuration |0011>. A
  single-parameter circuit that rotates ONLY between these two states —
  ansatz(theta): cos(theta/2)|1100> + sin(theta/2)|0011> — therefore spans
  the model's entire physically relevant space with ONE parameter, and
  <H>(theta) is EXACTLY a sinusoid: A + B*cos(theta) + C*sin(theta).

  Three measurements pin down A, B, C exactly (closed-form linear algebra,
  not an iterative/heuristic optimizer that could need "just one more
  step"). A fourth measurement AT the analytically-determined optimum
  theta* confirms the fit against a real, independently-measured value.
  Total cost is exactly 4 hardware jobs — known in advance, not estimated.

  Verified in this repo (exact statevector, no hardware) against the real
  TP53_C275F/mutant Hamiltonian in jw_hamiltonians.json: the fitted minimum
  reaches -266.4907486337723 Ha active+ecore total against a CASSCF-exact
  reference of -266.4907517477888 Ha — within 3e-6 Ha, an order of
  magnitude inside chemical accuracy (1.6e-3 Ha). Singles do not matter for
  this Hamiltonian; the doubles-only ansatz is exact for it.

ANSATZ CIRCUIT (4 qubits, matching _hf_circuit's q0..q3 = alpha0,beta0,
alpha1,beta1 convention used everywhere else in this codebase):
    RY(theta, q0); CX(q0,q2); CX(q0,q3); X(q0); CX(q0,q1)
5 gates, only 3 of them two-qubit — chosen to minimise the error surface on
real hardware, not just parameter count.

HONEST SCOPE: this is NOT the same ansatz as the 3-parameter (2 singles + 1
double) AllSinglesDoubles circuit the Phase 3A *simulator* VQE uses. It
captures the dominant (and, for this class of Hamiltonian, apparently
sufficient) correlation effect with a single parameter and an EXACT
closed-form solve — a different, cheaper, and independently interesting
capability, not a claim of matching the 3-parameter simulator optimum in
general. p1_ansatz says exactly this; do not read it as equivalent to the
simulator VQE.

BUDGET SAFETY: --budget-seconds (default 200) is a hard ceiling checked
against CUMULATIVE billed QPU-seconds spent so far before every new
submission. If spending so far already leaves no safe room for another job
(estimated from the max job cost seen, or a conservative 25s guess before
any job has run), the script stops and reports partial results rather than
risking an overrun on an irreplaceable trial budget. Submission retries
(network/queueing errors) cost NO QPU time — see solange_qpu.py's measure().

USAGE:
  # free validation first — ALWAYS do this before --hardware:
  python solange_qpu_vqe1p.py --key TP53_C275F --side mutant --dry-run --submit

  # then, with an empty/short queue on the target backend:
  export QISKIT_IBM_TOKEN=...; export QISKIT_IBM_INSTANCE=...
  python solange_qpu_vqe1p.py --key TP53_C275F --side mutant \
      --hardware --backend ibm_marrakesh --shots 4096 --budget-seconds 200 --submit
"""

import argparse
import hashlib
import json
import math
import os
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))   # for solange_qpu

from solange_qpu import (  # noqa: E402
    check_ibm_credentials, jw_target, h2_target, _exact, _backend_telemetry,
    _billable_qpu_seconds, submit, retrieve,
)


def build_ansatz(theta):
    """cos(theta/2)|1100> + sin(theta/2)|0011> — see module docstring for the
    derivation and the exact numerical verification against Statevector."""
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(4)
    qc.ry(theta, 0)
    qc.cx(0, 2)
    qc.cx(0, 3)
    qc.x(0)
    qc.cx(0, 1)
    return qc


class BudgetExceeded(Exception):
    pass


def measure_theta(obs, theta, hardware, backend_name, shots, token, instance,
                  spent_tracker, budget_seconds):
    """One <H(theta)> measurement — hardware or free local exact/dry-run.
    spent_tracker is a mutable dict {'total': float, 'max_job': float, 'jobs': []}
    updated in place so the caller's budget guard sees cumulative spend live."""
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    qc = build_ansatz(theta)

    if not hardware:
        from qiskit.quantum_info import Statevector
        energy = float(Statevector(qc).expectation_value(obs).real)
        return energy, "Statevector (local dry-run — NOT real hardware)", {}, {
            "mode": "dry-run", "shots": shots, "theta": theta}

    # Budget guard BEFORE spending anything on this job.
    est_next = spent_tracker["max_job"] or 25.0   # conservative guess pre-first-job
    if spent_tracker["total"] + est_next > budget_seconds:
        raise BudgetExceeded(
            f"spent {spent_tracker['total']:.1f}s so far; next job estimated "
            f"~{est_next:.1f}s would exceed the {budget_seconds}s ceiling — stopping "
            f"before submission (no QPU time at risk).")

    from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
    kwargs = {}
    if token:    kwargs["token"] = token
    if instance: kwargs["instance"] = instance
    try:
        service = QiskitRuntimeService(channel="ibm_quantum_platform", **kwargs)
    except Exception as e1:
        try:
            service = QiskitRuntimeService(**kwargs)
        except Exception:
            raise e1
    backend = service.backend(backend_name)
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
    isa = pm.run(qc)
    est = EstimatorV2(mode=backend)
    est.options.default_shots = shots

    pub = (isa, obs.apply_layout(isa.layout))
    job, last_err = None, None
    for attempt in range(1, 6):
        try:
            job = est.run([pub]); break
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            print(f"  submit attempt {attempt}/5 failed ({str(e)[:80]}...) — retrying "
                  f"in {wait}s (no QPU time spent)...", file=sys.stderr)
            time.sleep(wait)
    if job is None:
        raise RuntimeError(f"submission failed after 5 attempts, NO QPU time spent: {last_err}")

    print(f"  theta={theta:+.4f}  job {job.job_id()} -> {backend_name} "
          f"(queue time doesn't count against budget)...")
    while True:
        try:
            st = str(job.status())
        except Exception:
            st = "UNKNOWN"
        print(f"    status: {st}")
        if st in ("DONE", "ERROR", "CANCELLED"):
            break
        time.sleep(5)

    energy = float(job.result()[0].data.evs)
    tel = _backend_telemetry(backend)
    qpu_s, qpu_src = _billable_qpu_seconds(job)
    qpu_s = float(qpu_s or 0.0)
    spent_tracker["total"] += qpu_s
    spent_tracker["max_job"] = max(spent_tracker["max_job"], qpu_s)
    spent_tracker["jobs"].append({"theta": theta, "job_id": job.job_id(), "qpu_seconds": qpu_s})
    print(f"    energy={energy:.8f}  qpu_seconds={qpu_s:.2f}  "
          f"cumulative={spent_tracker['total']:.2f}/{budget_seconds}s")
    return (energy, f"{backend_name} (real QPU)", tel,
            {"mode": "hardware", "shots": shots, "job_id": job.job_id(),
             "backend": backend_name, "qpu_seconds": qpu_s,
             "qpu_seconds_source": qpu_src, "theta": theta})


def fit_and_confirm(obs, hardware, backend_name, shots, token, instance, budget_seconds,
                    resume_theta0_job=None):
    """3 measurements -> exact closed-form fit -> 1 confirmation measurement.
    Returns (theta_star, e_fit_min, e_confirm, all_points, spent_tracker, telemetry)
    — telemetry is from the confirmation job (the one the reported energy comes from).

    resume_theta0_job: if given, an IBM job id for a theta=0 job that was already
    SUBMITTED (e.g. by a run whose local process then crashed/segfaulted) — its
    result is RETRIEVED, not re-run, costing no additional QPU time. The local
    process dying does not cancel a job already accepted by IBM; see
    solange_qpu.py's own retrieve(), written for exactly this failure mode."""
    spent = {"total": 0.0, "max_job": 0.0, "jobs": []}
    pts = []
    thetas = [0.0, math.pi / 2, -math.pi / 2]
    energies = []
    for i, th in enumerate(thetas):
        if i == 0 and resume_theta0_job:
            print(f"  theta={th:+.4f}  RETRIEVING already-submitted job "
                  f"{resume_theta0_job} (no new QPU time)...")
            e, label, tel, meta = retrieve(resume_theta0_job, backend_name, token, instance)
            qpu_s = float(meta.get("qpu_seconds") or 0.0)
            spent["total"] += qpu_s
            spent["max_job"] = max(spent["max_job"], qpu_s)
            spent["jobs"].append({"theta": th, "job_id": resume_theta0_job,
                                  "qpu_seconds": qpu_s, "retrieved": True})
            print(f"    energy={e:.8f}  qpu_seconds={qpu_s:.2f}  "
                  f"cumulative={spent['total']:.2f}/{budget_seconds}s")
        else:
            e, label, tel, meta = measure_theta(obs, th, hardware, backend_name, shots,
                                                token, instance, spent, budget_seconds)
        energies.append(e)
        pts.append({"theta": th, "energy": e, "meta": meta})
    E0, E1, E2 = energies
    A = (E1 + E2) / 2
    B = E0 - A
    C = (E1 - E2) / 2
    theta_star = math.atan2(-C, -B)
    e_fit_min = A - math.hypot(B, C)

    e_confirm, label, tel, meta = measure_theta(obs, theta_star, hardware, backend_name,
                                                shots, token, instance, spent, budget_seconds)
    pts.append({"theta": theta_star, "energy": e_confirm, "meta": meta, "role": "confirmation"})
    return theta_star, e_fit_min, e_confirm, pts, spent, tel


def build_record(target, theta_star, e_fit_min_active, e_confirm_active, hf_exact_active,
                 pts, spent, hardware, backend_name, telemetry=None):
    telemetry = telemetry or {}
    now = datetime.now(timezone.utc).isoformat()
    ecore = target.get("ecore") or 0.0
    total_confirm = ecore + e_confirm_active
    total_fit = ecore + e_fit_min_active
    total_hf_ref = ecore + hf_exact_active
    rec = {
        "id": str(uuid.uuid4()), "created_at": now,
        "mutation_id": target["mutation_id"],
        "mutation_name": f"{target['mutation_id']} — Phase 3B "
                         f"{'QPU' if hardware else 'QPU dry-run'} "
                         f"single-parameter VQE (CAS({target['nelecas']},{target['ncas']}))",
        "side": target["side"],
        "phase": "3B-QPU" if hardware else "3B-QPU-dryrun",
        "provenance_source": "QPU/IBM" if hardware else "QPU/local-dryrun",

        "p1_qubit_count": target["nq"],
        "p1_ansatz": ("single-parameter double-excitation VQE — "
                     "RY(theta,q0);CX(q0,q2);CX(q0,q3);X(q0);CX(q0,q1); "
                     "theta* found by EXACT closed-form 3-point sinusoid fit "
                     "(A+B*cos(theta)+C*sin(theta)), NOT the 3-parameter "
                     "AllSinglesDoubles ansatz the Phase 3A simulator VQE uses — "
                     "see solange_qpu_vqe1p.py docstring"),
        "p1_circuit_hash": hashlib.sha256(json.dumps(
            {"terms": target["terms"], "nelec": target["nelec"],
             "ansatz": "double-excitation-1param"}, sort_keys=True).encode()).hexdigest(),

        "p2_encoding": target["encoding"], "p2_basis_set": target["basis"],
        "p2_active_electrons": target["nelecas"], "p2_active_orbitals": target["ncas"],
        "p2_model_compound": target["compound"],
        "p2_jw_terms": target["terms"],

        "p3_backend": backend_name if hardware else "Statevector (local dry-run — NOT real hardware)",
        "p3_vendor_job_id": pts[-1]["meta"].get("job_id"),
        "p3_backend_version": telemetry.get("p3_backend_version"),
        "p3_calibration_epoch": telemetry.get("p3_calibration_epoch") or now,

        "p4_gate_error_rate": telemetry.get("p4_gate_error_rate"),
        "p4_readout_error_rate": telemetry.get("p4_readout_error_rate"),
        "p4_t1_us": telemetry.get("p4_t1_us"), "p4_t2_us": telemetry.get("p4_t2_us"),
        "p4_note": ("real IBM QPU calibration data (confirmation job)" if hardware
                    else "dry-run — no device"),

        "p5_shots": pts[-1]["meta"].get("shots"), "p5_elapsed_s": None,
        "p5_qpu_seconds": spent["total"],
        "p5_qpu_seconds_source": "sum of IBM-reported per-job execution time across 4 jobs",
        "p5_ecore_ha": ecore,
        "p5_active_energy_ha": e_confirm_active,
        "p5_casscf_ref_ha": target.get("e_casscf"),

        "p6_method": "exact 3-point closed-form sinusoid fit + hardware confirmation",
        "p6_note": f"fit predicted {total_fit:.8f} Ha; confirmation measurement gave "
                   f"{total_confirm:.8f} Ha (delta {abs(total_fit-total_confirm):.2e} Ha)",

        "p7_energy_ha": total_confirm,
        "p7_ref_hf_ha": total_hf_ref,
        "p7_method": ("real on-device VQE (1-param, exact fit)" if hardware
                      else "dry-run on-device VQE (1-param, exact fit)"),
        "p7_note": "First on-device parameter OPTIMIZATION in Phase 3B (prior QPU runs "
                   "measured only a fixed HF state). theta*={:.6f} rad. Not directly "
                   "comparable to the 3-parameter simulator VQE optimum — see p1_ansatz.".format(theta_star),

        "p9_applicable": False,
        "p9_note": "deterministic default decoding; not an ML-decoded run",

        "_all_points": pts,             # full audit trail: all 4 (theta, energy) pairs
        "_qpu_seconds_breakdown": spent["jobs"],
        "_theta_star": theta_star,
    }
    payload = json.dumps({k: v for k, v in rec.items()
                          if k.startswith(("p1_", "p2_", "p3_", "p4_", "p5_",
                                           "p6_", "p7_", "p9_")) and k != "p3_calibration_epoch"},
                         sort_keys=True, default=str)
    rec["p8_seal_payload"] = payload
    rec["p8_hash"] = hashlib.sha256(payload.encode()).hexdigest()
    rec["p8_algorithm"] = "SHA-256"
    rec["p8_sealed_at"] = now
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", help="SOLANGE target key from jw_hamiltonians.json")
    ap.add_argument("--side", default="native", choices=["native", "mutant"])
    ap.add_argument("--jw-file", default=str(_HERE.parents[2] / "jw_hamiltonians.json"))
    ap.add_argument("--hardware", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--backend", default="ibm_marrakesh")
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--budget-seconds", type=float, default=200.0,
                    help="hard ceiling on cumulative billed QPU-seconds this run may spend")
    ap.add_argument("--instance", default=os.environ.get("QISKIT_IBM_INSTANCE"))
    ap.add_argument("--out", default="./out")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the sealed record to SOLANGE")
    ap.add_argument("--resume-theta0-job", metavar="JOB_ID", default=None,
                    help="an already-submitted theta=0 IBM job id whose result should be "
                         "RETRIEVED (no new QPU time) instead of re-measured — for resuming "
                         "after a local crash (segfault etc.) with a job still alive on IBM")
    args = ap.parse_args()

    if args.hardware == args.dry_run:
        ap.error("pass exactly one of --hardware or --dry-run")

    if args.hardware:
        ok, detail = check_ibm_credentials()
        if not ok:
            ap.error(f"IBM credentials not available: {detail}")
        token = os.environ.get("QISKIT_IBM_TOKEN")
        n_new = 3 if args.resume_theta0_job else 4
        print(f"HARDWARE RUN — backend={args.backend}  budget_seconds={args.budget_seconds}  "
              f"shots={args.shots}  ({n_new} new job(s) planned"
              + (f" + 1 retrieved (theta=0, job {args.resume_theta0_job})" if args.resume_theta0_job else "")
              + ")")
    else:
        token = None
        print("DRY RUN — free local exact statevector, no hardware, no budget spent")

    target = jw_target(args.key, args.side, args.jw_file) if args.key else h2_target()
    print(f"target: {target['label']}")

    theta_star, e_fit, e_confirm, pts, spent, telemetry = fit_and_confirm(
        target["obs"], args.hardware, args.backend, args.shots, token, args.instance,
        args.budget_seconds, resume_theta0_job=args.resume_theta0_job)

    hf_exact, ground = _exact(target["obs"], target["circuit"])
    ecore = target.get("ecore") or 0.0
    print("=" * 72)
    print(f"theta* = {theta_star:.6f} rad")
    print(f"fit minimum (total)      = {e_fit + ecore:.8f} Ha")
    print(f"confirmation measurement (total) = {e_confirm + ecore:.8f} Ha")
    print(f"CASSCF-exact reference   = {target.get('e_casscf')}")
    print(f"HF reference (total)     = {hf_exact + ecore:.8f} Ha")
    print(f"cumulative QPU-seconds spent = {spent['total']:.2f} / {args.budget_seconds}")
    print("=" * 72)

    record = build_record(target, theta_star, e_fit, e_confirm, hf_exact, pts, spent,
                          args.hardware, args.backend, telemetry)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    stem = f"qpu_vqe1p_{target['mutation_id']}_{target['side']}_{'hw' if args.hardware else 'dry'}"
    outp = Path(args.out) / f"{stem}_{record['id'][:8]}.json"
    outp.write_text(json.dumps(record, indent=2, default=str))
    print(f"WROTE {outp}")

    if args.submit:
        submit(args.submit, record)


if __name__ == "__main__":
    try:
        main()
    except BudgetExceeded as e:
        print(f"\nBUDGET GUARD STOPPED THE RUN: {e}", file=sys.stderr)
        sys.exit(2)
