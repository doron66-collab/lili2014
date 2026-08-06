#!/usr/bin/env python3
"""Pathway P2 of the Comparative Benchmark Protocol (BM-2026-001): Qiskit alone.

The reference task carried out with IBM's SDK and no orchestration layer, as a
competent researcher would actually write it. The chemistry still comes from
PySCF - Qiskit does not do quantum chemistry on its own, and pretending
otherwise would understate this pathway's step count rather than overstate it.

On the Runtime primitive. The protocol fixes a statevector backend and excludes
real hardware, because queue time is not attributable to any pathway. Qiskit's
local StatevectorEstimator implements the same primitive interface that Runtime
exposes for hardware, so the code below is what a Runtime submission looks like
minus the credential and the queue. That substitution is recorded here rather
than buried: it REMOVES steps from this pathway (no account setup, no backend
selection, no job polling), so the count it produces is a lower bound on what
the hardware route would cost.

Run:
    python scripts/benchmark/p2_qiskit.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

BOND_ANGSTROM = 0.7414
BASIS = "sto-3g"
NCAS = 2
NELECAS = 2
MAXITER = 5000
ANSATZ_REPS = 3          # raised from 1; see the note at the ansatz


def active_space_matrix():
    """Same active space and same qubit Hamiltonian as P4, so the two pathways
    are solving one problem rather than two similar ones."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from p4_manual import build_active_space, jordan_wigner_hamiltonian
    h1e, h2e, ecore, e_casci, e_hf = build_active_space()
    return jordan_wigner_hamiltonian(h1e, h2e, ecore), e_casci, e_hf


def main():
    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter
    from qiskit.quantum_info import SparsePauliOp
    from qiskit.primitives import StatevectorEstimator
    from qiskit.circuit.library import efficient_su2

    H, e_casci, e_hf = active_space_matrix()
    n_qubits = int(np.log2(H.shape[0]))

    # Dense matrix -> Pauli sum. Qiskit's own conversion, so the operator the
    # estimator sees is built by the pathway under test rather than by us.
    op = SparsePauliOp.from_operator(H)

    # Hardware-efficient ansatz. Qiskit's own EfficientSU2 rather than one
    # hand-rolled here, so the circuit under test is the library's.
    #
    # REPS IS NOT ARBITRARY. A single entangling layer - the first thing written
    # here - cannot express the doubly-excited configuration the correlated
    # ground state needs, and VQE converged twenty millihartree high, at the
    # Hartree-Fock energy, from every starting point tried. Nothing in the
    # pathway reported a problem: the number returned looked entirely
    # reasonable. Depth was raised until the energy reached chemical accuracy,
    # and the iterating is itself part of what this pathway costs.
    hf = QuantumCircuit(n_qubits)
    hf.x(0); hf.x(1)                                  # Hartree-Fock reference
    ansatz = efficient_su2(n_qubits, reps=ANSATZ_REPS, entanglement="full")
    qc = hf.compose(ansatz)
    theta = list(qc.parameters)

    estimator = StatevectorEstimator()
    calls = {"n": 0}

    def energy(params):
        calls["n"] += 1
        job = estimator.run([(qc, op, list(params))])
        return float(job.result()[0].data.evs.real)

    # NOT zeros. At theta = 0 this ansatz reproduces the Hartree-Fock state
    # exactly, and by its symmetry that point is stationary: COBYLA evaluates,
    # finds no improving direction, and returns the HF energy. The first run of
    # this script did exactly that - it reported -1.11668438 Ha, twenty
    # millihartree above the true ground state, with nothing anywhere in the
    # pathway flagging the answer as wrong. Recorded here because it is a
    # finding about the pathway, not an accident of this file: a bare SDK will
    # return a plausible number from a broken setup and say nothing.
    rng = np.random.default_rng(0)
    x0 = rng.normal(scale=0.1, size=len(theta))
    res = minimize(energy, x0, method="COBYLA", options={"maxiter": MAXITER})
    e_vqe = float(res.fun)

    print("=" * 66)
    print("P2 - Qiskit alone, reference task BM-2026-001 §3")
    print("=" * 66)
    print(f"  qubits            {n_qubits}   Pauli terms: {len(op)}")
    print(f"  ansatz            EfficientSU2, reps={ANSATZ_REPS}, {len(theta)} parameters")
    print(f"  estimator calls   {calls['n']}")
    print(f"  Hartree-Fock      {e_hf:.8f} Ha")
    print(f"  VQE energy        {e_vqe:.8f} Ha")
    print(f"  CASCI reference   {e_casci:.8f} Ha")
    err = abs(e_vqe - e_casci) * 1000
    print(f"  error vs exact    {err:.4f} mHa  {'within chemical accuracy' if err < 1.6 else '*** ABOVE 1.6 mHa ***'}")

    print()
    print("  Part 11 elements produced WITHOUT work beyond the counted steps:")
    part11 = [
        ("E1", "attributable",      False, "no identity is captured by the SDK"),
        ("E2", "timestamped",       False, "no timestamp is written; the operator would add one"),
        ("E3", "tamper-evident",    False, "results are plain values with no integrity check"),
        ("E4", "complete inputs",   True,  "circuit, operator and parameters are all in this file"),
        ("E5", "device state",      False, "a simulator reports none, and on hardware the operator "
                                           "must fetch calibration separately"),
        ("E6", "re-verifiable",     False, "nothing a third party can check without trusting the operator"),
        ("E7", "append-only trail", False, "no audit record; Runtime job history is IBM's, not the "
                                           "researcher's, and is not an append-only record of the science"),
    ]
    for eid, name, ok, why in part11:
        print(f"    {eid} {name:<18} {'yes' if ok else 'NO ':<4} - {why}")
    satisfied = sum(1 for *_, ok, _ in part11 if ok)
    print(f"    -> {satisfied} of 7 satisfied")

    out = Path("out/benchmark"); out.mkdir(parents=True, exist_ok=True)
    (out / "p2_qiskit_result.json").write_text(json.dumps({
        "pathway": "P2 Qiskit alone",
        "n_qubits": n_qubits, "pauli_terms": len(op), "ansatz_reps": ANSATZ_REPS,
        "estimator_calls": calls["n"],
        "e_hf": e_hf, "e_vqe": e_vqe, "e_casci": e_casci, "error_mha": err,
        "part11_satisfied": satisfied,
        "part11": {e: {"name": n, "satisfied": o, "note": w} for e, n, o, w in part11},
    }, indent=2))
    print(f"\n  wrote {out/'p2_qiskit_result.json'}")
    print("=" * 66)
    return 0


if __name__ == "__main__":
    sys.exit(main())
