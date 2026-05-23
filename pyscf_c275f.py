"""
PySCF CASSCF(2,2) ground-state computation for TP53 C275F mutation site.

Molecular model:
  Phenylalanine sidechain proxy: toluene (PhCH3)
  Represents the Phe275 aromatic π system introduced by the C275F substitution.
  Standard model compound for amino acid sidechain QM/MM partitioning.
  Reference: Gao J & Truhlar DG (2002) Annu Rev Phys Chem 53:467–505

Active space:
  CASSCF(2,2) — 2 electrons in 2 orbitals (π/π* HOMO–LUMO pair)
  This is the Phase 3A proxy active space.
  Jordan-Wigner encoding: 4 spin-orbitals = 4 qubits (IBM Sherbrooke hardware run)

Geometry:
  Idealized from Engh & Huber (1991) Acta Cryst A47:392–400
  Cβ–Cγ (sp3→aromatic): 1.512 Å
  Aromatic C–C: 1.396 Å (benzene), ring angles 120°
  C–H: 1.083 Å (aromatic), 1.090 Å (sp3)

Basis:
  STO-3G (minimal Slater-type orbital, 3 Gaussians)
  Appropriate for Phase 3A proof-of-concept; Phase 3B will use cc-pVDZ or 6-31G*

References:
  Cho Y et al. (1994) Science 265:346–355 — TP53 crystal structure (1TUP)
  Joerger AC et al. (2006) PNAS 103:15056–15061 — C275F structure (2BIM)
  Joerger AC & Fersht AR (2008) Annu Rev Biochem 77:557–579 — TP53 review
  Sun Q et al. (2020) J Chem Phys 153:024109 — PySCF 2.0 paper
"""

import numpy as np

# ─── Geometry ─────────────────────────────────────────────────────────────────
# Toluene: methyl group (Cβ proxy) + benzene ring (Phe275 sidechain π system)
# Coordinate frame: ring in xz-plane, Cβ below ring along z-axis
# All positions in Angstroms.
#
# Atom layout:
#   [0]  Cβ   — sp3 carbon (Cα–Cβ attachment point, methyl proxy)
#   [1]  Cγ   — ring C1 (ipso, bonded to Cβ)
#   [2]  Cδ1  — ring C2 (ortho)
#   [3]  Cε1  — ring C3 (meta)
#   [4]  Cζ   — ring C4 (para)
#   [5]  Cε2  — ring C5 (meta)
#   [6]  Cδ2  — ring C6 (ortho)
#   [7-9]  Hβ (3 methyl H's)
#   [10-14] Hδ1, Hε1, Hζ, Hε2, Hδ2 (5 ring H's)

TOLUENE_ATOMS = [
    # Cβ (methyl carbon, Cα attachment proxy)
    ('C',  ( 0.000,  0.000,  0.000)),
    # Ring carbons (benzene, C–C = 1.396 Å, circumradius = 1.396 Å)
    # Cγ = ipso, placed 1.512 Å above Cβ along z
    ('C',  ( 0.000,  0.000,  1.512)),   # Cγ  (ipso)
    ('C',  ( 1.209,  0.000,  2.210)),   # Cδ1 (ortho, +x side)
    ('C',  ( 1.209,  0.000,  3.602)),   # Cε1 (meta,  +x side)
    ('C',  ( 0.000,  0.000,  4.300)),   # Cζ  (para)
    ('C',  (-1.209,  0.000,  3.602)),   # Cε2 (meta,  -x side)
    ('C',  (-1.209,  0.000,  2.210)),   # Cδ2 (ortho, -x side)
    # Methyl H's (tetrahedral around Cβ, staggered w.r.t. ring)
    ('H',  ( 0.000,  1.026, -0.363)),   # Hβ1
    ('H',  ( 0.889, -0.513, -0.363)),   # Hβ2
    ('H',  (-0.889, -0.513, -0.363)),   # Hβ3
    # Ring H's (C–H = 1.083 Å, pointing radially outward from ring)
    ('H',  ( 2.156,  0.000,  1.672)),   # Hδ1
    ('H',  ( 2.156,  0.000,  4.140)),   # Hε1
    ('H',  ( 0.000,  0.000,  5.383)),   # Hζ
    ('H',  (-2.156,  0.000,  4.140)),   # Hε2
    ('H',  (-2.156,  0.000,  1.672)),   # Hδ2
]

def build_mol_string(atoms):
    lines = []
    for sym, (x, y, z) in atoms:
        lines.append(f'{sym}  {x:.4f}  {y:.4f}  {z:.4f}')
    return '\n'.join(lines)


def run_casscf_phase3a(verbose=3):
    """
    Run RHF → CASSCF(2,2) for the Phe275 sidechain model.

    Returns dict with energies and integrals for VQE Pauli decomposition.
    """
    try:
        from pyscf import gto, scf, mcscf
    except ImportError:
        print("\n  PySCF not installed. Install with: pip install pyscf")
        print("  Returning placeholder values used in Phase 3A proxy run.\n")
        return _placeholder_result()

    print("\n" + "="*60)
    print("  TP53 C275F — Phase 3A PySCF CASSCF(2,2)")
    print("  Model: Phe275 sidechain (toluene proxy, STO-3G)")
    print("="*60)

    # Build molecule
    mol = gto.Mole()
    mol.atom = build_mol_string(TOLUENE_ATOMS)
    mol.basis = 'sto-3g'
    mol.charge = 0
    mol.spin = 0       # closed-shell singlet
    mol.verbose = verbose
    mol.build()

    print(f"\n  Molecule: toluene (C7H8)")
    print(f"  Atoms: {mol.natm}, Electrons: {mol.nelectron}, Basis functions: {mol.nao}")

    # Step 1: RHF reference wavefunction
    print("\n  [1/3] Running RHF...")
    mf = scf.RHF(mol)
    mf.max_cycle = 200
    mf.conv_tol = 1e-10
    e_rhf = mf.kernel()
    print(f"  RHF converged: {mf.converged}")
    print(f"  RHF energy:    {e_rhf:.10f} Ha")

    # Step 2: CASSCF(2,2) — 2 electrons in 2 π orbitals (HOMO/LUMO)
    # For toluene, HOMO and LUMO are the degenerate π1/π1* pair of the ring.
    # PySCF selects the active orbitals automatically from MO ordering.
    print("\n  [2/3] Running CASSCF(2,2)...")
    mc = mcscf.CASSCF(mf, ncas=2, nelecas=2)
    mc.conv_tol = 1e-9
    mc.conv_tol_grad = 1e-5
    mc.max_cycle_macro = 100
    e_casscf = mc.kernel()[0]
    print(f"  CASSCF converged: {mc.converged}")
    print(f"  CASSCF energy:    {e_casscf:.10f} Ha")

    correlation_energy = e_casscf - e_rhf

    # Step 3: Extract active-space integrals for VQE Pauli decomposition
    print("\n  [3/3] Extracting active-space Hamiltonian integrals...")
    h1e, ecore = mc.get_h1eff()
    h2e = mc.get_h2eff()

    # Natural orbital occupancies (diagnostic for active space quality)
    natocc = mc.mo_occ

    print(f"\n{'='*60}")
    print(f"  RESULTS — TP53 C275F Phase 3A CASSCF(2,2)")
    print(f"{'='*60}")
    print(f"  Model compound:     Phe275 sidechain (toluene, C7H8)")
    print(f"  Basis set:          STO-3G (minimal)")
    print(f"  Active space:       CAS(2e, 2o) — Phe π/π* HOMO–LUMO")
    print(f"  JW encoding:        4 spin-orbitals → 4 qubits")
    print(f"")
    print(f"  RHF energy:         {e_rhf:+.8f} Ha")
    print(f"  CASSCF energy:      {e_casscf:+.8f} Ha")
    print(f"  Correlation energy: {correlation_energy:+.8f} Ha")
    print(f"")
    print(f"  Active-space 1e integrals (h1e), shape: {h1e.shape}")
    print(f"  h1e =")
    print(f"    [[{h1e[0,0]:+.8f}  {h1e[0,1]:+.8f}]")
    print(f"     [{h1e[1,0]:+.8f}  {h1e[1,1]:+.8f}]]")
    print(f"")
    print(f"  Core energy (frozen electrons): {ecore:+.8f} Ha")
    print(f"")
    # h2e from get_h2eff() is stored compressed (n_pair × n_pair).
    # Restore to full 4D tensor for readable indexing.
    from pyscf import ao2mo
    h2e_full = ao2mo.restore(1, h2e, mc.ncas)   # shape (ncas,ncas,ncas,ncas)
    print(f"  Active-space 2e integrals (h2e), shape: {h2e_full.shape} (full 4D)")
    print(f"  h2e[0,0,0,0] = {h2e_full[0,0,0,0]:+.8f} Ha  (Coulomb integral J_aa)")
    print(f"  h2e[0,0,1,1] = {h2e_full[0,0,1,1]:+.8f} Ha  (Coulomb integral J_ab)")
    print(f"  h2e[0,1,1,0] = {h2e_full[0,1,1,0]:+.8f} Ha  (Exchange integral K_ab)")
    print(f"{'='*60}")

    result = {
        'model':              'toluene (Phe275 sidechain proxy)',
        'basis':              'STO-3G',
        'active_space':       'CAS(2e,2o)',
        'jw_qubits':          4,
        'e_rhf_Ha':           e_rhf,
        'e_casscf_Ha':        e_casscf,
        'correlation_Ha':     correlation_energy,
        'ecore_Ha':           float(ecore),
        'h1e':                h1e.tolist(),
        'h2e_J_aa':           float(h2e_full[0, 0, 0, 0]),
        'h2e_J_ab':           float(h2e_full[0, 0, 1, 1]),
        'h2e_K_ab':           float(h2e_full[0, 1, 1, 0]),
        'converged':          bool(mc.converged),
        'pauli_note': (
            'Decompose h1e + h2e into Pauli operators via Jordan-Wigner transform '
            'using openfermion or qiskit-nature for VQE circuit construction.'
        ),
    }
    return result


def _placeholder_result():
    """Placeholder values matching the proxy Hamiltonian used in Phase 3A runs."""
    return {
        'model':          'toluene (Phe275 sidechain proxy)',
        'basis':          'STO-3G',
        'active_space':   'CAS(2e,2o)',
        'jw_qubits':      4,
        'e_rhf_Ha':       None,
        'e_casscf_Ha':    None,
        'correlation_Ha': None,
        'converged':      False,
        'note':           'PySCF not available — install with: pip install pyscf',
    }


def pauli_decomposition_note():
    print("""
  ─── VQE PAULI DECOMPOSITION (next step after this script) ───────────────
  The integrals h1e and h2e define the molecular Hamiltonian H in the
  active orbital basis. To run VQE on IBM hardware:

  1. Jordan-Wigner transform (4 qubits for 2e/2o CAS):
       from openfermion import jordan_wigner, InteractionOperator
       # or: from qiskit_nature.second_q.mappers import JordanWignerMapper

  2. Pauli decomposition yields a sum of Pauli strings:
       H = Σ_i c_i · P_i   where P_i ∈ {I,X,Y,Z}^⊗4

  3. VQE ansatz: hardware-efficient Ry–CNOT–Ry (Phase 3A proxy)
       or UCCSD (unitary coupled-cluster, Phase 3B)

  4. Phase 3A run: IBM Sherbrooke (127-qubit Eagle r3)
       backend = service.backend('ibm_sherbrooke')
       vqe = VQE(ansatz, optimizer=COBYLA(), estimator=Estimator())
       result = vqe.compute_minimum_eigenvalue(qubit_op)

  5. Phase 3B: Same pipeline, 88-qubit active space on IBM Heron r3
       (pending access authorization — see SOLANGE Phase 3B submission pkg)
  ──────────────────────────────────────────────────────────────────────────
""")


# ─── Active space electron count (matches casscf_analysis.py methodology) ─────
def print_c275f_active_space_summary():
    print("""
  TP53 C275F — CASSCF Active Space Summary
  ─────────────────────────────────────────
  Native C275:
    Valence electrons (5 Å shell): see casscf_analysis.py
    Key active orbitals: Cys Sγ lone pair (nS) + S–H σ* + backbone π
    H-bond donor: Sγ → R248 Nε (3.40 Å) — stabilizes L3 loop

  Mutant C275F:
    Valence electrons (5 Å shell): same count as native (ΔE = 0)
    Key active orbitals: Phe aromatic π / π* (HOMO–LUMO)
    No H-bond donor capacity → R248 displaced → LOF
    Phase 3A CAS(2e,2o): HOMO/LUMO of Phe π system (this script)
    Phase 3B CAS(44e,44o): full DBD functional interface (88 qubits)

  Jordan-Wigner encoding:
    Phase 3A:  2e × 2 spin-orbitals/orbital = 4 qubits
    Phase 3B: 44e × 2 spin-orbitals/orbital = 88 qubits

  PDB references:
    Native: 2AC0 (2.05 Å), 1TUP (2.35 Å, with DNA)
    Mutant: 2BIM (1.80 Å) — Joerger et al. (2006) PNAS 103:15056
""")


if __name__ == '__main__':
    print_c275f_active_space_summary()
    result = run_casscf_phase3a(verbose=3)
    pauli_decomposition_note()

    import json
    out_file = 'c275f_casscf_phase3a.json'
    with open(out_file, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n  Integrals written to: {out_file}")
    print("  Feed h1e, h2e, ecore into openfermion/qiskit-nature for VQE.\n")
