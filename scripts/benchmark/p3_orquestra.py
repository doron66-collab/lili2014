#!/usr/bin/env python3
"""Pathway P3 of the Comparative Benchmark Protocol (BM-2026-001): Orquestra.

The protocol required this pathway's availability to be confirmed by obtaining
and running it, not inferred from Zapata Computing's wind-down in October 2024.
It was obtained: orquestra-core 0.12.0 installs from PyPI and its quantum, vqa
and opt subpackages all import. The pathway is therefore assessable, and is
assessed here rather than written off.

Two costs of it are real and are recorded rather than smoothed over.

It performs no quantum chemistry, so the Hamiltonian is built by
export_hamiltonian.py and carried in as a file. And it cannot share a Python
environment with the Qiskit pathway - Orquestra pins numpy below 2, Qiskit 2
requires 2 or above, and installing one breaks the other at runtime with an
error that never mentions numpy. Each pathway needs its own interpreter, which
is itself a setup step.

Run (from the repo root, in the isolated environment):
    /tmp/venv_orquestra/bin/python scripts/benchmark/p3_orquestra.py
"""
import json
import sys
from pathlib import Path

import numpy as np

MAXITER = 5000
ANSATZ_REPS = 3          # matched to P2 so the two ansaetze have equal depth
CHEM_ACC_MHA = 1.6


def load_hamiltonian():
    src = Path("out/benchmark/reference_hamiltonian.json")
    if not src.exists():
        sys.exit("run scripts/benchmark/export_hamiltonian.py first (main environment)")
    return json.loads(src.read_text())


def build_ansatz(circuits, n_qubits, params):
    """EfficientSU2-equivalent, assembled from primitives.

    Orquestra ships QAOA and singlet-UCCSD ansaetze but no generic
    hardware-efficient one, so the circuit P2 got from a library call is built
    by hand here. That asymmetry is a finding about the pathway and is counted.
    """
    Circuit, RY, RZ, CNOT, X = (circuits.Circuit, circuits.RY, circuits.RZ,
                                circuits.CNOT, circuits.X)
    ops = [X(0), X(1)]                                   # Hartree-Fock reference
    k = 0
    for _ in range(ANSATZ_REPS):
        for q in range(n_qubits):
            ops.append(RY(params[k])(q)); k += 1
            ops.append(RZ(params[k])(q)); k += 1
        for a in range(n_qubits):                        # full entanglement, as in P2
            for b in range(a + 1, n_qubits):
                ops.append(CNOT(a, b))
    for q in range(n_qubits):
        ops.append(RY(params[k])(q)); k += 1
        ops.append(RZ(params[k])(q)); k += 1
    # n_qubits must be declared: Orquestra sizes a circuit from the qubits its
    # gates touch, and a parameter set that leaves one idle would silently
    # produce a wavefunction of the wrong dimension.
    return Circuit(ops, n_qubits=n_qubits)


def main():
    from orquestra.quantum import circuits
    from orquestra.quantum.operators import PauliSum, PauliTerm, get_sparse_operator
    from orquestra.quantum.runners import SymbolicSimulator
    from scipy.optimize import minimize

    data = load_hamiltonian()
    n_qubits = data["n_qubits"]
    e_casci, e_hf = data["e_casci"], data["e_hf"]

    # Pauli strings -> Orquestra's own operator type.
    #
    # The index is NOT reversed, and getting this wrong is silent. The first
    # version mapped qiskit's string position i to qubit n-1-i, on the reasoning
    # that qiskit prints qubit 0 rightmost. The resulting operator had exactly
    # the right spectrum - its lowest eigenvalue matched CASCI to eight decimals,
    # because relabelling qubits is a permutation and a permutation does not move
    # eigenvalues - so a spectrum check passed it. What it broke was the STATE:
    # the Hartree-Fock preparation X(0)X(1) then pointed at a basis state of the
    # relabelled operator carrying 0.459 Ha instead of -1.117, and VQE spent its
    # whole budget climbing back to a reference it should have started from.
    #
    # The check that catches this is below and now runs on every execution:
    # evaluate the prepared HF state against the assembled operator and require
    # it to reproduce the mean-field energy. An operator can be right and a
    # convention still wrong, and only a state-level test can tell them apart.
    terms = []
    for t in data["terms"]:
        s = t["pauli"]
        label = "*".join(f"{c}{i}" for i, c in enumerate(s) if c != "I")
        terms.append(PauliTerm(label, t["coeff"]) if label else PauliTerm("I0", t["coeff"]))
    H = PauliSum(terms)

    n_params = 2 * n_qubits * (ANSATZ_REPS + 1)
    runner = SymbolicSimulator()
    # Built once. The first version of this called it inside the loop, which is
    # both slow and pointless - the Hamiltonian does not change with the
    # parameters.
    H_sparse = get_sparse_operator(H)

    # Convention gate. Prepare Hartree-Fock and require the assembled operator to
    # return the mean-field energy PySCF computed. Everything downstream is
    # meaningless if this fails, so it aborts rather than warns.
    hf_circ = circuits.Circuit([circuits.X(0), circuits.X(1)], n_qubits=n_qubits)
    psi_hf = np.asarray(runner.get_wavefunction(hf_circ).amplitudes, dtype=complex)
    e_hf_check = float(np.real(np.conj(psi_hf) @ (H_sparse @ psi_hf)))
    if abs(e_hf_check - e_hf) > 1e-8:
        sys.exit(f"ABORT: qubit convention mismatch. <HF|H|HF> = {e_hf_check:.8f} Ha "
                 f"but PySCF gives {e_hf:.8f} Ha. The operator's spectrum can still be "
                 f"correct while the basis ordering is not.")
    print(f"  convention check  <HF|H|HF> = {e_hf_check:.8f} Ha == PySCF HF  OK")

    calls = {"n": 0}

    def energy(theta):
        calls["n"] += 1
        circ = build_ansatz(circuits, n_qubits, list(theta))
        psi = np.asarray(runner.get_wavefunction(circ).amplitudes, dtype=complex)
        return float(np.real(np.conj(psi) @ (H_sparse @ psi)))

    rng = np.random.default_rng(0)
    x0 = rng.normal(scale=0.1, size=n_params)
    res = minimize(energy, x0, method="COBYLA", options={"maxiter": MAXITER})
    e_vqe = float(res.fun)
    err = abs(e_vqe - e_casci) * 1000

    print("=" * 66)
    print("P3 - Orquestra, reference task BM-2026-001 §3")
    print("=" * 66)
    print(f"  availability      CONFIRMED — orquestra-core 0.12.0 installed and ran")
    print(f"  qubits            {n_qubits}   Pauli terms: {len(data['terms'])}")
    print(f"  ansatz            hand-built SU2, reps={ANSATZ_REPS}, {n_params} parameters")
    print(f"  simulator calls   {calls['n']}")
    print(f"  Hartree-Fock      {e_hf:.8f} Ha")
    print(f"  VQE energy        {e_vqe:.8f} Ha")
    print(f"  CASCI reference   {e_casci:.8f} Ha")
    print(f"  error vs exact    {err:.4f} mHa  "
          f"{'within chemical accuracy' if err < CHEM_ACC_MHA else '*** ABOVE 1.6 mHa ***'}")

    print()
    print("  Part 11 elements produced WITHOUT work beyond the counted steps:")
    part11 = [
        ("E1", "attributable",      False, "no identity is captured by the SDK"),
        ("E2", "timestamped",       False, "workflow runs can be timestamped by the platform "
                                           "service, which is what ceased operating; the local "
                                           "SDK writes none"),
        ("E3", "tamper-evident",    False, "results are plain values with no integrity check"),
        ("E4", "complete inputs",   True,  "operator and circuit are fixed in this file, the "
                                           "Hamiltonian in the file it reads"),
        ("E5", "device state",      False, "the simulator reports none"),
        ("E6", "re-verifiable",     False, "nothing a third party can check without trusting "
                                           "the operator"),
        ("E7", "append-only trail", False, "no audit record is written locally"),
    ]
    for eid, name, ok, why in part11:
        print(f"    {eid} {name:<18} {'yes' if ok else 'NO ':<4} - {why}")
    satisfied = sum(1 for *_, ok, _ in part11 if ok)
    print(f"    -> {satisfied} of 7 satisfied")

    out = Path("out/benchmark"); out.mkdir(parents=True, exist_ok=True)
    (out / "p3_orquestra_result.json").write_text(json.dumps({
        "pathway": "P3 Orquestra",
        "availability": "confirmed — orquestra-core 0.12.0, isolated environment required",
        "n_qubits": n_qubits, "pauli_terms": len(data["terms"]),
        "ansatz_reps": ANSATZ_REPS, "simulator_calls": calls["n"],
        "e_hf": e_hf, "e_vqe": e_vqe, "e_casci": e_casci, "error_mha": err,
        "part11_satisfied": satisfied,
        "part11": {e: {"name": n, "satisfied": o, "note": w} for e, n, o, w in part11},
    }, indent=2))
    print(f"\n  wrote {out/'p3_orquestra_result.json'}")
    print("=" * 66)
    return 0


if __name__ == "__main__":
    sys.exit(main())
