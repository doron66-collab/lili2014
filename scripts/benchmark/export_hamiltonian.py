#!/usr/bin/env python3
"""Export the reference task's qubit Hamiltonian as Pauli strings.

Orquestra performs no quantum chemistry, so the Hamiltonian has to be built
elsewhere and carried in. That hand-off is a real cost of the P3 pathway and is
counted as a step, not hidden inside it. It also cannot share a Python
environment with the Qiskit pathway - Orquestra pins numpy below 2 and Qiskit 2
requires 2 or above - so the transfer is through a file rather than an import.
"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p4_manual import build_active_space, jordan_wigner_hamiltonian
from qiskit.quantum_info import SparsePauliOp

h1e, h2e, ecore, e_casci, e_hf = build_active_space()
H = jordan_wigner_hamiltonian(h1e, h2e, ecore)
op = SparsePauliOp.from_operator(H)
terms = [{"pauli": str(p), "coeff": float(c.real)}
         for p, c in zip(op.paulis, op.coeffs) if abs(c) > 1e-12]

out = Path("out/benchmark"); out.mkdir(parents=True, exist_ok=True)
(out / "reference_hamiltonian.json").write_text(json.dumps(
    {"n_qubits": op.num_qubits, "e_casci": e_casci, "e_hf": e_hf, "terms": terms}, indent=2))
print(f"wrote {out/'reference_hamiltonian.json'}: {len(terms)} Pauli terms, "
      f"{op.num_qubits} qubits, CASCI = {e_casci:.8f} Ha")
