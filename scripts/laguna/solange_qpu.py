#!/usr/bin/env python3
"""
solange_qpu.py — Phase 3B: real IBM Quantum hardware runs + P1-P9 capture.

PURPOSE (honest scope): prove the MOVING PARTS of the Phase 3B pipeline end to
end on REAL quantum hardware — circuit → IBM QPU → real telemetry → P1-P9
provenance record → LEON notarization in SOLANGE. These are hardware capability
demonstrations, NOT qualified scientific results.

Two targets:
  • H2 (default) — canonical minimal-basis H2 (O'Malley 2016), 2 qubits.
  • --key/--side — one of SOLANGE's OWN 4-qubit CAS(2,2) model-compound
    Hamiltonians from jw_hamiltonians.json (e.g. TP53_C275F native/mutant). This
    runs the dissertation's actual target Hamiltonian on real hardware — still a
    minimal CAS(2,2) active space, NOT the full 44e/88q anchor, and a fixed-state
    <H> measurement with NO on-device VQE optimization (which would burn a scarce
    QPU-time budget for a value already known classically).

Budget-aware (a trial account may have only ~10 minutes of QPU time):
  --dry-run  (DEFAULT) FREE local simulator — validate everything here first.
  --hardware sends ONE job to a real QPU.

USAGE:
  # free validation:
  python solange_qpu.py --key TP53_C275F --side native --dry-run --submit
  # then one real run:
  export QISKIT_IBM_TOKEN=...; export QISKIT_IBM_INSTANCE=...
  python solange_qpu.py --key TP53_C275F --side native --hardware --backend ibm_marrakesh --submit
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

# Canonical minimal-basis (STO-3G) H2 Hamiltonian, JW + 2-qubit reduction, 0.735 Å
# — O'Malley et al., PRX 6:031007 (2016). Terms are (coeff, pauli-label q0q1).
H2_TERMS = [
    (-1.052373245772859,  "II"),
    ( 0.397936358509854,  "IZ"),
    (-0.397936358509854,  "ZI"),
    (-0.011280104256235,  "ZZ"),
    ( 0.180931350397663,  "XX"),
]
_REPO_ROOT = Path(__file__).resolve().parents[2]
_JW_DEFAULT = _REPO_ROOT / "jw_hamiltonians.json"


def check_ibm_credentials():
    """Whether qiskit-ibm-runtime is installed and a token/account is configured.
    No network call. Returns (available, detail)."""
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
    """P9 (ML-decoder lineage). This uses the hardware's default deterministic
    decoding, so P9 is the minimal form — honestly NOT an ML-decoded run."""
    return {"decoder": decoder, "decoder_version": None, "ml_decoded": False}


# ── Target construction ────────────────────────────────────────────────────────
def _pauli_terms_to_op(terms, nq):
    """Build a SparsePauliOp from SOLANGE's {coeff, pauli} term list. pauli is
    'I' (global identity), 'Z0', or space-separated like 'X0 Y1 Z2'."""
    from qiskit.quantum_info import SparsePauliOp
    sparse = []
    for t in terms:
        p, c = t["pauli"], t["coeff"]
        if p in ("I", ""):
            sparse.append(("", [], c))
        else:
            letters, idxs = "", []
            for tok in p.split():
                letters += tok[0]
                idxs.append(int(tok[1:]))
            sparse.append((letters, idxs, c))
    return SparsePauliOp.from_sparse_list(sparse, num_qubits=nq)


def _hf_circuit(nelec, nq):
    """Fixed Hartree-Fock reference: occupy the lowest `nelec` spin-orbitals
    (JW ordering, matching PennyLane hf_state used by the classical runs). No
    parameters, no VQE — one deterministic state prep."""
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(nq)
    for q in range(nelec):
        qc.x(q)
    return qc


def h2_target():
    from qiskit.quantum_info import SparsePauliOp
    obs = SparsePauliOp([p for _, p in H2_TERMS], [c for c, _ in H2_TERMS])
    return {
        "label": "H2 (canonical STO-3G, 0.735 Å — O'Malley 2016)",
        "mutation_id": "H2_smoketest", "compound": "H2", "side": "native",
        "nq": 2, "nelec": 2, "obs": obs, "circuit": _hf_circuit(2, 2),
        "terms": [{"coeff": c, "pauli": p} for c, p in H2_TERMS],
        "ecore": 0.0, "e_casscf": None,
        "encoding": "Jordan-Wigner + 2-qubit reduction (O'Malley 2016)",
        "basis": "STO-3G", "nelecas": 2, "ncas": 2,
    }


def jw_target(key, side, jw_file):
    """Load one of SOLANGE's own CAS(2,2) model-compound Hamiltonians."""
    data = json.loads(Path(jw_file).read_text())
    if key not in data:
        raise SystemExit(f"key '{key}' not in {jw_file}. Available: {', '.join(sorted(data))}")
    if side not in data[key]:
        raise SystemExit(f"side '{side}' not in {key} (have: {', '.join(data[key])})")
    e = data[key][side]
    terms = e["terms"]
    nq = 1 + max((int(tok[1:]) for t in terms for tok in t["pauli"].split()
                  if t["pauli"] not in ("I", "")), default=1)
    nelecas = e.get("nelecas", 2)
    ncas = e.get("ncas", nq // 2)
    return {
        "label": f"{key}/{side} — CAS({nelecas},{ncas}) model compound '{e.get('compound')}'",
        "mutation_id": key, "compound": e.get("compound"), "side": side,
        "nq": nq, "nelec": nelecas, "obs": _pauli_terms_to_op(terms, nq),
        "circuit": _hf_circuit(nelecas, nq), "terms": terms,
        "ecore": e.get("ecore", 0.0), "e_casscf": e.get("e_casscf"),
        "e_active_exact": e.get("e_active_exact"), "e_active_rhf": e.get("e_active_rhf"),
        "encoding": f"Jordan-Wigner (CASSCF({nelecas},{ncas}) -> {len(terms)} Pauli terms)",
        "basis": e.get("basis", "STO-3G"), "nelecas": nelecas, "ncas": ncas,
    }


def _exact(obs, circuit):
    """Exact <H> on the fixed reference state, and exact ground state — computed
    locally so the hardware number has an honest, self-consistent reference."""
    import numpy as np
    from qiskit.quantum_info import Statevector
    hf = float(np.real(Statevector(circuit).expectation_value(obs)))
    ground = float(np.linalg.eigvalsh(obs.to_matrix())[0])
    return hf, ground


def _backend_telemetry(backend):
    """Extract REAL P3/P4 provenance from a live backend — defensively (missing
    field → None). This is the point of Phase 3B: real calibration data, not
    simulated placeholders."""
    tel = {"p3_calibration_epoch": None, "p4_gate_error_rate": None,
           "p4_readout_error_rate": None, "p4_t1_us": None, "p4_t2_us": None,
           "p3_backend_version": None}
    try:
        props = backend.properties()
    except Exception:
        props = None
    if props is not None:
        for key, fn in (("p3_calibration_epoch", lambda: props.last_update_date.isoformat()),
                        ("p4_readout_error_rate", lambda: float(props.readout_error(0))),
                        ("p4_t1_us", lambda: round(float(props.t1(0)) * 1e6, 3)),
                        ("p4_t2_us", lambda: round(float(props.t2(0)) * 1e6, 3))):
            try:
                tel[key] = fn()
            except Exception:
                pass
        try:
            errs = sorted(g.parameters[0].value for g in props.gates
                          if g.gate in ("ecr", "cz", "cx") and g.parameters)
            if errs:
                tel["p4_gate_error_rate"] = float(errs[len(errs) // 2])
        except Exception:
            pass
    try:
        tel["p3_backend_version"] = str(getattr(backend, "version", None) or
                                        getattr(backend, "backend_version", None))
    except Exception:
        pass
    return tel


def measure(target, hardware, backend_name, shots, token, instance):
    """Measure <H> of the target on its fixed reference state. dry-run → local
    simulator (free); hardware → one real QPU job."""
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    obs, qc = target["obs"], target["circuit"]

    if not hardware:
        try:
            from qiskit_aer import AerSimulator
            from qiskit_aer.primitives import EstimatorV2 as AerEstimator
            backend = AerSimulator()
            pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
            isa = pm.run(qc)
            job = AerEstimator().run([(isa, obs.apply_layout(isa.layout))],
                                     precision=1 / math.sqrt(shots))
            energy = float(job.result()[0].data.evs)
            label = "AerSimulator (local dry-run — NOT real hardware)"
        except Exception:
            hf, _ = _exact(obs, qc)      # exact statevector fallback (still free)
            energy, label = hf, "Statevector (local dry-run fallback — NOT real hardware)"
        return energy, label, {"p3_calibration_epoch": None}, {"mode": "dry-run", "shots": shots}

    from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
    kwargs = {}
    if token:    kwargs["token"] = token
    if instance: kwargs["instance"] = instance
    try:
        service = QiskitRuntimeService(channel="ibm_quantum_platform", **kwargs)
    except Exception:
        service = QiskitRuntimeService(channel="ibm_quantum", **kwargs)
    backend = service.backend(backend_name)
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
    isa = pm.run(qc)
    est = EstimatorV2(mode=backend)
    est.options.default_shots = shots
    # IBM's job-queueing API can return transient 5xx / connection errors even when
    # the backend shows operational (the "Error queueing job" / "too many 500" we hit).
    # Retry the SUBMISSION with backoff so a flaky window doesn't kill the run — a
    # failed submit spends NO QPU time, and queue time never counts against budget.
    pub = (isa, obs.apply_layout(isa.layout))
    job, last_err = None, None
    for attempt in range(1, 6):
        try:
            job = est.run([pub]); break
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            print(f"  submit attempt {attempt}/5 failed ({str(e)[:80]}…) — retrying in {wait}s "
                  f"(no QPU time spent)…", file=sys.stderr)
            time.sleep(wait)
    if job is None:
        raise RuntimeError(
            "IBM job submission failed after 5 attempts — almost certainly an IBM-side "
            "queueing outage (us-east), not your code/network. NO QPU time was spent. "
            f"Check quantum.cloud.ibm.com status and retry later. Last error: {last_err}")
    print(f"  submitted job {job.job_id()} to {backend_name} — waiting for result "
          f"(queue time does NOT count against QPU-execution budget) …")
    energy = float(job.result()[0].data.evs)
    tel = _backend_telemetry(backend)
    return (energy, f"{backend_name} (real QPU)", tel,
            {"mode": "hardware", "shots": shots, "job_id": job.job_id(), "backend": backend_name})


def retrieve(job_id, backend_name, token, instance):
    """Fetch the result of an ALREADY-COMPLETED QPU job by id and return it in the
    same shape as measure(). Used when the local script hung on job.result() even
    though IBM finished the job — the result lives on IBM's servers, and retrieving
    it costs NO QPU execution time. No new circuit is run."""
    from qiskit_ibm_runtime import QiskitRuntimeService
    kwargs = {}
    if token:    kwargs["token"] = token
    if instance: kwargs["instance"] = instance
    try:
        service = QiskitRuntimeService(channel="ibm_quantum_platform", **kwargs)
    except Exception:
        service = QiskitRuntimeService(channel="ibm_quantum", **kwargs)
    job = service.job(job_id)
    try:
        status = job.status()
    except Exception:
        status = "?"
    print(f"  retrieving completed job {job_id} (status: {status}) — no QPU time spent …")
    energy = float(job.result()[0].data.evs)
    bname = backend_name
    try:
        b = job.backend()
        bname = b if isinstance(b, str) else getattr(b, "name", backend_name)
    except Exception:
        pass
    tel = {}
    try:
        tel = _backend_telemetry(service.backend(bname))
    except Exception:
        pass
    return (energy, f"{bname} (real QPU)", tel,
            {"mode": "hardware", "shots": None, "job_id": job_id, "backend": bname})


def build_record(target, active_energy, hf_exact_active, backend_label, telemetry, meta):
    """Assemble + seal a P1-P9 record. Energies are reported at the ACTIVE-space
    level and also as totals (ecore + active), matching the classical runs'
    convention. Field names use p1_..p9_ so LEON's existing P8 seal covers them.

    hf_exact_active is the exact classical <H> on the SAME fixed HF state the QPU
    measured — the honest reference for a fixed-state run. The gap between the
    measured value and THIS (not the CASSCF ground state) is the actual hardware
    noise; stored as p7_ref_hf_ha so the ladder can show that honestly instead of
    conflating hardware noise with the intended HF-vs-ground gap."""
    now = datetime.now(timezone.utc).isoformat()
    is_hw = meta.get("mode") == "hardware"
    ecore = target.get("ecore") or 0.0
    total_energy = ecore + active_energy
    total_hf_ref = ecore + hf_exact_active
    rec = {
        "id": str(uuid.uuid4()), "created_at": now,
        "mutation_id": target["mutation_id"],
        "mutation_name": f"{target['mutation_id']} — Phase 3B "
                         f"{'QPU' if is_hw else 'QPU dry-run'} (CAS({target['nelecas']},{target['ncas']}))",
        "side": target["side"],
        "phase": "3B-QPU" if is_hw else "3B-QPU-dryrun",
        "provenance_source": "QPU/IBM" if is_hw else "QPU/local-dryrun",

        "p1_qubit_count": target["nq"],
        "p1_ansatz": "fixed HF reference (no VQE optimization on device)",
        "p1_circuit_hash": hashlib.sha256(json.dumps(
            {"terms": target["terms"], "nelec": target["nelec"]},
            sort_keys=True).encode()).hexdigest(),

        "p2_encoding": target["encoding"], "p2_basis_set": target["basis"],
        "p2_active_electrons": target["nelecas"], "p2_active_orbitals": target["ncas"],
        "p2_model_compound": target["compound"],
        "p2_jw_terms": target["terms"],

        "p3_backend": backend_label,
        "p3_backend_version": telemetry.get("p3_backend_version"),
        "p3_calibration_epoch": telemetry.get("p3_calibration_epoch") or now,

        "p4_gate_error_rate": telemetry.get("p4_gate_error_rate"),
        "p4_readout_error_rate": telemetry.get("p4_readout_error_rate"),
        "p4_t1_us": telemetry.get("p4_t1_us"), "p4_t2_us": telemetry.get("p4_t2_us"),
        "p4_note": ("real IBM QPU calibration data" if is_hw else "dry-run — no device"),

        "p5_shots": meta.get("shots"), "p5_elapsed_s": None,
        "p5_ecore_ha": ecore,
        "p5_active_energy_ha": active_energy,
        "p5_casscf_ref_ha": target.get("e_casscf"),

        "p6_method": "none (raw estimator)", "p6_note": "hardware smoke test — no error mitigation",

        "p7_energy_ha": total_energy,
        "p7_ref_hf_ha": total_hf_ref,   # exact <H> on the SAME fixed HF state — the
                                        # honest reference; |measured - this| = hardware noise
        "p7_method": ("QPU EstimatorV2 <H> on HF ref (total = ecore + active)" if is_hw
                      else "dry-run EstimatorV2 <H> on HF ref"),
        "p7_note": ("Fixed-state expectation, NOT a ground-state search — Δ to the CASSCF "
                    "ground is expected. Hardware noise = |measured - exact HF| (p7_ref_hf_ha). "
                    "CAS(2,2) minimal active space, NOT the full anchor."),

        "p9_applicable": False,
        "p9_note": "deterministic default decoding; not an ML-decoded run",
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


def submit(api, record):
    import urllib.request
    body = {"provenance": record, "jw": {}}
    try:
        req = urllib.request.Request(api.rstrip("/") + "/api/simulate/hpc/submit",
                                     data=json.dumps(body, default=str).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            resp = json.loads(r.read().decode())
        print(f"SUBMITTED → SOLANGE  status={resp.get('status')}  seal_ok={resp.get('seal_ok')}  "
              f"db={resp.get('db_status')}  notary={resp.get('notary')}  run_id={resp.get('run_id')}")
        return resp
    except Exception as e:
        print(f"  SUBMIT FAILED (result still saved locally): {e}", file=sys.stderr)
        return None


def run_one_qpu(key, side, backend, shots, token, instance, jw_file, out_dir, api):
    """Execute ONE real-hardware QPU job end-to-end (used by --agent): build the
    target, measure ⟨H⟩ on hardware, seal the record, save it, and submit to SOLANGE.
    Returns (submit_response_or_None, record)."""
    target = jw_target(key, side, jw_file) if key else h2_target()
    energy, backend_label, telemetry, meta = measure(target, True, backend, shots, token, instance)
    hf_exact, _ground = _exact(target["obs"], target["circuit"])
    record = build_record(target, energy, hf_exact, backend_label, telemetry, meta)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    stem = f"qpu_{target['mutation_id']}_{target['side']}_hw"
    (Path(out_dir) / f"{stem}_{record['id'][:8]}.json").write_text(
        json.dumps(record, indent=2, default=str))
    return submit(api, record), record


def run_agent(api, backend, shots, poll_s, jw_file, instance, out_dir, email, password):
    """QPU pull-agent — the same governed model as the Laguna HPC agent, but for real
    hardware. Runs ONLY in your authenticated session; claims ONLY jobs you queued in
    SOLANGE (job_type=qpu). Each claimed job spends real quantum time on `backend`."""
    import urllib.request
    from solange_hpc import _supabase_login, _ts
    ibm_token = os.environ.get("QISKIT_IBM_TOKEN")   # None ⇒ use the saved IBM account
    jwt = _supabase_login(email, password)
    hdr = {"Authorization": "Bearer " + jwt}

    def _post(path, body):
        req = urllib.request.Request(api.rstrip("/") + path,
            data=json.dumps(body).encode(), method="POST",
            headers={**hdr, "Content-Type": "application/json"})
        return urllib.request.urlopen(req, timeout=30).read()

    def heartbeat():
        try: _post("/api/simulate/hpc/agent/heartbeat", {"agent": "qpu", "agent_id": "qpu"})
        except Exception: pass

    def post_status(did, status, note=None, run_id=None):
        body = {"status": status}
        if note:   body["note"] = str(note)[:400]
        if run_id: body["run_id"] = run_id
        try: _post(f"/api/simulate/hpc/dispatch/{did}/status", body)
        except Exception as e: print(f"[{_ts()}] [qpu-agent] status post failed: {e}", file=sys.stderr)

    print(f"[{_ts()}] [qpu-agent] up · backend={backend} · polling {api} every {poll_s}s · Ctrl+C to stop")
    print(f"[{_ts()}] [qpu-agent] ⚠ each claimed job spends REAL quantum time on {backend}.")
    while True:
        heartbeat()
        try:
            req = urllib.request.Request(
                api.rstrip("/") + "/api/simulate/hpc/dispatch/next?job_type=qpu", headers=hdr)
            with urllib.request.urlopen(req, timeout=30) as r:
                job = json.loads(r.read().decode()).get("job")
        except Exception as e:
            print(f"[{_ts()}] [qpu-agent] poll error: {e}", file=sys.stderr); time.sleep(poll_s); continue
        if not job:
            time.sleep(poll_s); continue
        did = job["id"]; key = job.get("key"); side = job.get("side", "native")
        print(f"[{_ts()}] [qpu-agent] job {did[:8]} · {key}/{side} → REAL QPU on {backend}")
        try:
            resp, _record = run_one_qpu(key, side, backend, shots, ibm_token, instance, jw_file, out_dir, api)
            stored = bool(resp) and resp.get("db_status") in ("stored", "stored_no_payload")
            ok = bool(stored and resp.get("seal_ok"))
            note = "ok" if ok else ("verified but not stored" if resp else "submit failed")
            post_status(did, "done" if ok else "failed", note=note,
                        run_id=(resp.get("run_id") if resp else None))
            print(f"[{_ts()}] [qpu-agent] job {did[:8]} → {'DONE' if ok else 'FAILED'}  "
                  f"db={resp.get('db_status') if resp else '—'}")
        except KeyboardInterrupt:
            raise
        except (Exception, SystemExit) as e:
            # jw_target() raises SystemExit for a clean CLI error message on a
            # single direct run (--key ... without --agent) — correct there, but
            # SystemExit does NOT inherit from Exception, so a bare 'except
            # Exception' here would let it escape uncaught and silently kill the
            # whole agent process (looks like it "just stopped", no traceback,
            # nothing left polling). Catch it explicitly so one bad key just fails
            # that job and the agent keeps polling for the next one.
            note = str(e) or repr(e)
            post_status(did, "failed", note=note)
            print(f"[{_ts()}] [qpu-agent] job {did[:8]} error: {note}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="SOLANGE Phase 3B — real QPU hardware smoke test.")
    ap.add_argument("--key", help="SOLANGE target key from jw_hamiltonians.json "
                                   "(e.g. TP53_C275F). Omit to run the default H2.")
    ap.add_argument("--side", default="native", choices=["native", "mutant"])
    ap.add_argument("--jw-file", default=str(_JW_DEFAULT))
    ap.add_argument("--hardware", action="store_true",
                    help="send ONE job to a real QPU (spends QPU time). Omit for a free dry-run.")
    ap.add_argument("--dry-run", action="store_true",
                    help="explicit free local-simulator run (already the default without --hardware).")
    ap.add_argument("--backend", default="ibm_marrakesh",
                    help="QPU name (with --hardware). ibm_marrakesh = Heron r2, open trial instance.")
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--instance", default=os.environ.get("QISKIT_IBM_INSTANCE"))
    ap.add_argument("--out", default="./out")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the sealed record to SOLANGE (LEON notarizes it).")
    ap.add_argument("--check-credentials", action="store_true")
    ap.add_argument("--retrieve", metavar="JOB_ID",
                    help="fetch the result of an already-COMPLETED QPU job by id (e.g. when "
                         "the local run hung on job.result()). Costs NO QPU time. Requires "
                         "--key/--side/--backend to rebuild and submit the record.")
    # ── QPU pull-agent (same governed model as the Laguna HPC agent) ──────────
    ap.add_argument("--agent", action="store_true",
                    help="run as the QPU pull-agent: poll SOLANGE for queued job_type=qpu "
                         "jobs and run each on --backend. Each claimed job spends REAL "
                         "quantum time. Runs only in your authenticated session.")
    ap.add_argument("--api", default="https://qcaihpc-simulation-api.onrender.com",
                    help="SOLANGE backend base URL (agent mode)")
    ap.add_argument("--poll", type=int, default=15, help="agent poll interval seconds")
    ap.add_argument("--email", default="guest@solange.bio", help="SOLANGE login (agent mode)")
    ap.add_argument("--password", default="Solange2026", help="SOLANGE password (agent mode)")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    token = os.environ.get("QISKIT_IBM_TOKEN")

    print("=" * 72)
    print("SOLANGE Phase 3B — real QPU hardware smoke test (moving-parts proof, not chemistry)")
    available, detail = check_ibm_credentials()
    print(f"IBM Quantum credentials: {'AVAILABLE' if available else 'NOT AVAILABLE'} — {detail}")
    if args.check_credentials:
        print("=" * 72); sys.exit(0 if available else 1)
    if (args.hardware or args.retrieve) and not available:
        print("-" * 72)
        print("REFUSING hardware/retrieve without credentials — no fabricated result. "
              "Set QISKIT_IBM_TOKEN and retry.")
        print("=" * 72); sys.exit(1)
    if args.retrieve and not args.key:
        ap.error("--retrieve requires --key/--side to rebuild the P1-P9 record")

    if args.agent:
        if not available:
            print("-" * 72)
            print("REFUSING to start the QPU agent without credentials — it would only fail "
                  "every claimed job. Set QISKIT_IBM_TOKEN (or save an account) and retry.")
            print("=" * 72); sys.exit(1)
        print(f"Mode: QPU PULL-AGENT · backend={args.backend}")
        print("=" * 72)
        run_agent(args.api, args.backend, args.shots, args.poll, args.jw_file,
                  args.instance, args.out, args.email, args.password)
        return

    target = jw_target(args.key, args.side, args.jw_file) if args.key else h2_target()
    _mode = ('RETRIEVE (fetch completed job — no QPU time)' if args.retrieve
             else 'REAL HARDWARE · backend='+args.backend if args.hardware
             else 'DRY-RUN (free local simulator)')
    print(f"Mode: {_mode}")
    print(f"Target: {target['label']}  ·  {target['nq']} qubits")

    if args.retrieve:
        energy, backend_label, telemetry, meta = retrieve(
            args.retrieve, args.backend, token, args.instance)
    else:
        energy, backend_label, telemetry, meta = measure(
            target, args.hardware, args.backend, args.shots, token, args.instance)
    hf_exact, ground = _exact(target["obs"], target["circuit"])
    record = build_record(target, energy, hf_exact, backend_label, telemetry, meta)

    ecore = target.get("ecore") or 0.0
    print("-" * 72)
    print(f"backend            : {backend_label}")
    print(f"measured <H> active: {energy:.6f} Ha   (on the fixed HF reference state)")
    print(f"exact <H> active HF: {hf_exact:.6f} Ha   (what a faithful pipeline approaches)")
    print(f"exact ground active: {ground:.6f} Ha   (reference only; no optimization was done)")
    if ecore:
        print(f"total energy (P7)  : {ecore + energy:.6f} Ha   (ecore {ecore:.6f} + active)")
    print(f"Δ(measured-HF)     : {(energy - hf_exact)*1000:.2f} mHa   (hardware noise; expected, not an error)")
    if telemetry.get("p3_calibration_epoch"):
        print(f"REAL P3 epoch      : {telemetry['p3_calibration_epoch']}")
        print(f"REAL P4 budget     : gate_err={telemetry.get('p4_gate_error_rate')} "
              f"readout_err={telemetry.get('p4_readout_error_rate')} "
              f"T1={telemetry.get('p4_t1_us')}us T2={telemetry.get('p4_t2_us')}us")
    print(f"P8 seal            : {record['p8_hash'][:16]}…")

    stem = f"qpu_{target['mutation_id']}_{target['side']}_" + \
           ("hw" if meta.get("mode") == "hardware" else "dryrun")
    p = Path(args.out) / f"{stem}_{record['id'][:8]}.json"
    p.write_text(json.dumps(record, indent=2, default=str))
    print(f"WROTE {p}")
    if args.submit:
        submit(args.submit, record)
    print("=" * 72)
    if not args.hardware:
        print("Dry-run clean? Re-run with --hardware (and --submit) for ONE real QPU job.")


if __name__ == "__main__":
    main()
