#!/usr/bin/env python3
"""Pathway P4 of the Comparative Benchmark Protocol (BM-2026-001): fully manual.

This is the baseline a researcher without any platform actually occupies. It
computes the reference task with nothing but a quantum-chemistry library and
linear algebra: build the molecule, reduce to the active space, map to qubits,
and evaluate the energy.

It is written to be HONEST about the baseline rather than to make it look bad.
Everything a competent researcher would in fact automate is automated here -
the whole file is one command invocation. What it cannot automate away is the
record: whatever this script prints, nobody can later show it was not edited,
nor who ran it, nor on what. That absence is the measurement, and §3 below
counts it explicitly rather than leaving it to be inferred.

Run:
    python scripts/benchmark/p4_manual.py
"""
import json
import sys
from pathlib import Path

import numpy as np

# The reference task, fixed by the protocol (BM-2026-001 §3).
BOND_ANGSTROM = 0.7414          # H2 at equilibrium
BASIS = "sto-3g"
NCAS = 2                        # CAS(2e,2o) -> 4 qubits under Jordan-Wigner
NELECAS = 2


def build_active_space():
    """Molecule -> active-space integrals. The chemistry half of the task."""
    from pyscf import gto, scf, mcscf, ao2mo

    mol = gto.M(atom=f"H 0 0 0; H 0 0 {BOND_ANGSTROM}", basis=BASIS,
                spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol).run()
    mc = mcscf.CASCI(mf, NCAS, NELECAS)
    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), NCAS)
    e_exact = mc.kernel()[0]
    return h1e, np.asarray(h2e), float(ecore), float(e_exact), float(mf.e_tot)


def jordan_wigner_hamiltonian(h1e, h2e, ecore):
    """Second-quantised Hamiltonian -> a Pauli operator on 2*NCAS qubits.

    Written out rather than imported, because importing it is exactly the
    convenience a manual pathway does not have. Spin-orbital p is qubit p, with
    even indices alpha and odd beta.
    """
    n_so = 2 * NCAS
    I2 = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    Sp = np.array([[0, 1], [0, 0]], dtype=complex)      # |0><1|, the annihilator
    Sm = Sp.conj().T

    def kron(ops):
        out = np.array([[1]], dtype=complex)
        for o in ops:
            out = np.kron(out, o)
        return out

    def annihilate(p):
        # Jordan-Wigner: Z on every qubit before p enforces fermionic sign.
        return kron([Z] * p + [Sp] + [I2] * (n_so - p - 1))

    a = [annihilate(p) for p in range(n_so)]
    ad = [x.conj().T for x in a]

    dim = 2 ** n_so
    H = np.zeros((dim, dim), dtype=complex)
    H += ecore * np.eye(dim, dtype=complex)

    # One-body: spatial orbital i,j carried by both spins.
    for i in range(NCAS):
        for j in range(NCAS):
            for s in (0, 1):
                H += h1e[i, j] * (ad[2 * i + s] @ a[2 * j + s])

    # Two-body, chemist notation (ij|kl).
    for i in range(NCAS):
        for j in range(NCAS):
            for k in range(NCAS):
                for l in range(NCAS):
                    v = 0.5 * h2e[i, j, k, l]
                    if v == 0:
                        continue
                    for s in (0, 1):
                        for t in (0, 1):
                            H += v * (ad[2 * i + s] @ ad[2 * k + t]
                                      @ a[2 * l + t] @ a[2 * j + s])
    return H


def main():
    h1e, h2e, ecore, e_casci, e_hf = build_active_space()
    H = jordan_wigner_hamiltonian(h1e, h2e, ecore)

    evals = np.linalg.eigvalsh(H)
    e_ground = float(evals[0].real)

    print("=" * 66)
    print("P4 - fully manual pathway, reference task BM-2026-001 §3")
    print("=" * 66)
    print(f"  system            H2 @ {BOND_ANGSTROM} A / {BASIS}")
    print(f"  active space      CAS({NELECAS},{NCAS}) -> {2*NCAS} qubits (Jordan-Wigner)")
    print(f"  Hartree-Fock      {e_hf:.8f} Ha")
    print(f"  qubit-H ground    {e_ground:.8f} Ha")
    print(f"  CASCI reference   {e_casci:.8f} Ha")
    agree = abs(e_ground - e_casci) * 1000
    print(f"  agreement         {agree:.6f} mHa  {'OK' if agree < 1e-3 else '*** MISMATCH ***'}")

    # The metric this pathway exists to expose. Each answer is a fact about what
    # the script above does, not a judgement about the pathway.
    print()
    print("  Part 11 elements produced WITHOUT work beyond the counted steps:")
    part11 = [
        ("E1", "attributable",      False, "no identity is captured; the shell user is whoever ran it"),
        ("E2", "timestamped",       False, "nothing records when this ran unless the operator notes it"),
        ("E3", "tamper-evident",    False, "stdout can be edited afterwards and nothing detects it"),
        ("E4", "complete inputs",   True,  "geometry, basis and active space are fixed in this file"),
        ("E5", "device state",      False, "no backend state is captured"),
        ("E6", "re-verifiable",     False, "a third party must trust the operator's copy of the output"),
        ("E7", "append-only trail", False, "no audit record of any kind is written"),
    ]
    for eid, name, ok, why in part11:
        print(f"    {eid} {name:<18} {'yes' if ok else 'NO ':<4} - {why}")
    satisfied = sum(1 for *_, ok, _ in [(e, n, o, w) for e, n, o, w in part11] if ok)
    print(f"    -> {satisfied} of 7 satisfied")

    out = Path("out/benchmark"); out.mkdir(parents=True, exist_ok=True)
    (out / "p4_manual_result.json").write_text(json.dumps({
        "pathway": "P4 fully manual",
        "e_hf": e_hf, "e_ground_qubit_hamiltonian": e_ground, "e_casci": e_casci,
        "agreement_mha": agree,
        "part11_satisfied": satisfied,
        "part11": {e: {"name": n, "satisfied": o, "note": w} for e, n, o, w in part11},
    }, indent=2))
    print(f"\n  wrote {out/'p4_manual_result.json'}")
    print("=" * 66)
    return 0 if agree < 1e-3 else 1


if __name__ == "__main__":
    sys.exit(main())
