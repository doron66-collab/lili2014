"""
pyscf_mutations.py — PySCF CASSCF(2,2) for all NSCLC mutations on the SOLANGE platform.

For each mutation, we compute the ground-state energy of a small model compound
that represents the chemically active sidechain of the native and mutant residue.
The CASSCF(2,2) active space targets the HOMO/LUMO pair — the two orbitals most
responsible for the change in electronic structure at the mutation site.

Model compounds (sidechain proxies):
  Cys  → methanethiol   (CH3SH)        active space: S lone pair / S–H σ*
  Phe  → toluene        (C7H8)         active space: aromatic π / π* (HOMO–LUMO)
  Tyr  → p-cresol       (C7H8O)        active space: phenol π / π* (HOMO–LUMO)
  Gly  → formamide      (HCONH2)       active space: backbone C=O π / π*
  Arg  → guanidine      (CH5N3)        active space: guanidinium π / π*
  Gln  → acetamide      (C2H5NO)       active space: amide C=O π / π*
  Asp  → acetic acid    (C2H4O2)       active space: carboxylate C=O π / π*
  Asn  → acetamide      (C2H5NO)       active space: amide C=O π / π*
  Leu  → isobutane      (C4H10)        active space: C–H σ / σ* (aliphatic, no π)

Basis: STO-3G (minimal — appropriate for Phase 3A proof of concept)
Method: RHF → CASSCF(2,2) → active-space integral extraction

Geometries: idealized from Engh & Huber (1991) Acta Cryst A47:392–400.
            Bond lengths and angles rounded to 3 decimal places.

References:
  Engh RA & Huber R (1991) Acta Cryst A47:392 — standard amino acid geometry
  Veryazov V et al. (2011) Int J Quantum Chem 111:3329 — active space selection
  Sun Q et al. (2020) J Chem Phys 153:024109 — PySCF 2.0
"""

import numpy as np
import json
from datetime import datetime

# ── Model compound geometries (Angstroms) ─────────────────────────────────────

GEOM = {}

# ── Toluene (C7H8) — Phe sidechain proxy ──────────────────────────────────────
# Ring in xz-plane. Cβ (methyl) below ring, ring carbons at 1.396 Å C–C.
# Cβ–Cγ (sp3→aryl): 1.512 Å. Aromatic C–H: 1.083 Å. sp3 C–H: 1.090 Å.
GEOM['toluene'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # Cβ (methyl)
    ('C', ( 0.000,  0.000,  1.512)),  # Cγ (ipso)
    ('C', ( 1.209,  0.000,  2.210)),  # Cδ1 (ortho)
    ('C', ( 1.209,  0.000,  3.602)),  # Cε1 (meta)
    ('C', ( 0.000,  0.000,  4.300)),  # Cζ  (para)
    ('C', (-1.209,  0.000,  3.602)),  # Cε2 (meta)
    ('C', (-1.209,  0.000,  2.210)),  # Cδ2 (ortho)
    ('H', ( 1.026,  0.000, -0.363)),  # Hβ1
    ('H', (-0.513,  0.889, -0.363)),  # Hβ2
    ('H', (-0.513, -0.889, -0.363)),  # Hβ3
    ('H', ( 2.156,  0.000,  1.672)),  # Hδ1
    ('H', ( 2.156,  0.000,  4.140)),  # Hε1
    ('H', ( 0.000,  0.000,  5.383)),  # Hζ
    ('H', (-2.156,  0.000,  4.140)),  # Hε2
    ('H', (-2.156,  0.000,  1.672)),  # Hδ2
]

# ── Methanethiol (CH3SH) — Cys sidechain proxy ────────────────────────────────
# C–S: 1.819 Å. S–H: 1.339 Å. H–S–C: 96.5°. C–H: 1.090 Å.
GEOM['methanethiol'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # Cβ
    ('S', ( 0.000,  0.000,  1.819)),  # Sγ
    ('H', ( 1.330,  0.000,  1.971)),  # Hγ (thiol H; H–S–C = 96.5°)
    ('H', ( 1.026,  0.000, -0.363)),  # Hβ1
    ('H', (-0.513,  0.889, -0.363)),  # Hβ2
    ('H', (-0.513, -0.889, -0.363)),  # Hβ3
]

# ── p-Cresol (C7H8O, 4-methylphenol) — Tyr sidechain proxy ───────────────────
# Toluene ring + phenolic OH at para. C–O: 1.364 Å. O–H: 0.963 Å. C–O–H: 109.5°.
GEOM['p_cresol'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # Cβ (methyl)
    ('C', ( 0.000,  0.000,  1.512)),  # Cγ (ipso)
    ('C', ( 1.209,  0.000,  2.210)),  # Cδ1
    ('C', ( 1.209,  0.000,  3.602)),  # Cε1
    ('C', ( 0.000,  0.000,  4.300)),  # Cζ (para — has OH)
    ('C', (-1.209,  0.000,  3.602)),  # Cε2
    ('C', (-1.209,  0.000,  2.210)),  # Cδ2
    ('H', ( 1.026,  0.000, -0.363)),  # Hβ1
    ('H', (-0.513,  0.889, -0.363)),  # Hβ2
    ('H', (-0.513, -0.889, -0.363)),  # Hβ3
    ('H', ( 2.156,  0.000,  1.672)),  # Hδ1
    ('H', ( 2.156,  0.000,  4.140)),  # Hε1
    ('H', (-2.156,  0.000,  4.140)),  # Hε2
    ('H', (-2.156,  0.000,  1.672)),  # Hδ2
    ('O', ( 0.000,  0.000,  5.664)),  # phenolic O
    ('H', ( 0.908,  0.000,  5.985)),  # phenolic H (C–O–H = 109.5°)
]

# ── Formamide (HCONH2) — Gly backbone proxy ───────────────────────────────────
# Smallest amide — represents the peptide backbone of Gly (no sidechain).
# Planar molecule. C=O: 1.193 Å. C–N: 1.352 Å. C–H: 1.102 Å. N–H: 1.002 Å.
# O–C–N: 124.2°. H–C–N: 112.7°.
GEOM['formamide'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # carbonyl C
    ('O', ( 0.000,  1.193,  0.000)),  # C=O (double bond, up)
    ('N', ( 1.121, -0.756,  0.000)),  # amide N (O–C–N = 124.2°)
    ('H', (-0.925, -0.600,  0.000)),  # H on C (H–C–N = 112.7°, opposite side from O)
    ('H', ( 2.014, -0.300,  0.000)),  # H on N (H–N–C = 121°)
    ('H', ( 1.034, -1.754,  0.000)),  # H on N (H–N–H = 122°)
]

# ── Acetic acid (CH3COOH) — Asp sidechain proxy ───────────────────────────────
# Carboxylate model. C=O: 1.214 Å. C–O: 1.364 Å. O–H: 0.972 Å. Carboxyl planar.
GEOM['acetic_acid'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # methyl C
    ('C', ( 1.520,  0.000,  0.000)),  # carbonyl C
    ('O', ( 2.127,  1.051,  0.000)),  # C=O (double bond)
    ('O', ( 2.139, -1.215,  0.000)),  # C–OH (single bond)
    ('H', ( 1.440, -1.890,  0.000)),  # O–H (C–O–H ≈ 107°)
    ('H', (-0.390,  1.026,  0.000)),  # methyl H1
    ('H', (-0.390, -0.513,  0.889)),  # methyl H2
    ('H', (-0.390, -0.513, -0.889)),  # methyl H3
]

# ── Acetamide (CH3CONH2) — Asn/Gln sidechain proxy ───────────────────────────
# Amide model. C=O: 1.220 Å. C–N: 1.380 Å (partial double bond). Amide group planar.
GEOM['acetamide'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # methyl C
    ('C', ( 1.522,  0.000,  0.000)),  # carbonyl C
    ('O', ( 2.151,  1.046,  0.000)),  # C=O
    ('N', ( 2.105, -1.250,  0.000)),  # amide N
    ('H', ( 3.100, -1.355,  0.000)),  # N–H1 (H–N–C = 121°)
    ('H', ( 1.546, -2.079,  0.000)),  # N–H2
    ('H', (-0.390,  1.026,  0.000)),  # methyl H1
    ('H', (-0.390, -0.513,  0.889)),  # methyl H2
    ('H', (-0.390, -0.513, -0.889)),  # methyl H3
]

# ── Guanidine (CH5N3) — Arg sidechain proxy ───────────────────────────────────
# Guanidinium π system: all three N–C bonds have partial double-bond character.
# Planar (D3h-like). C–N: 1.352 Å. N–H: 1.010 Å. N–C–N: 120°.
GEOM['guanidine'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # central C
    ('N', ( 1.352,  0.000,  0.000)),  # N1 (=NH, imine — 1 H only)
    ('N', (-0.676,  1.171,  0.000)),  # N2 (–NH2)
    ('N', (-0.676, -1.171,  0.000)),  # N3 (–NH2)
    ('H', ( 1.857,  0.875,  0.000)),  # N1–H  (H–N–C = 120°, in plane)
    ('H', (-0.200,  2.050,  0.000)),  # N2–H1
    ('H', (-1.660,  1.110,  0.000)),  # N2–H2
    ('H', (-0.200, -2.050,  0.000)),  # N3–H1
    ('H', (-1.660, -1.110,  0.000)),  # N3–H2
]
# Formula: CH5N3  →  6 + 5 + 21 = 32 electrons (even, spin=0 ✓)

# ── Isobutane (C4H10) — Leu sidechain proxy ───────────────────────────────────
# Aliphatic — no π system. CAS(2,2) targets C–H σ/σ* of the tertiary C–H bond.
# C–C: 1.522 Å. C–H: 1.090 Å. Tetrahedral geometry.
# NOTE: This is the physically least meaningful active space in this set.
#       Leu has no π electrons. The correlation energy will be near zero.
GEOM['isobutane'] = [
    ('C', ( 0.000,  0.000,  0.000)),  # central CH (tertiary)
    ('C', ( 0.000,  0.000,  1.522)),  # methyl 1
    ('C', ( 1.436,  0.000, -0.507)),  # methyl 2
    ('C', (-0.718,  1.243, -0.507)),  # methyl 3
    ('H', ( 0.000,  0.000, -1.090)),  # H on central C
    ('H', ( 1.028,  0.000,  1.885)),  # methyl 1 H's (tetrahedral)
    ('H', (-0.514,  0.890,  1.885)),
    ('H', (-0.514, -0.890,  1.885)),
    ('H', ( 2.121,  0.000,  0.342)),  # methyl 2 H's
    ('H', ( 1.607,  0.891, -1.113)),
    ('H', ( 1.607, -0.891, -1.113)),
    ('H', (-1.050,  1.819,  0.285)),  # methyl 3 H's
    ('H', ( 0.308,  0.916, -0.870)),
    ('H', (-1.745,  1.526, -0.870)),
]

# ── Mutation → model compound mapping ─────────────────────────────────────────
MUTATION_MODELS = {
    'TP53_Y220C':  {
        'label':  'TP53 Y220C (p.Tyr220Cys)',
        'gene':   'TP53',
        'native': {'compound': 'p_cresol',    'residue': 'Tyr220', 'note': 'Phenol π/π* HOMO–LUMO'},
        'mutant': {'compound': 'methanethiol','residue': 'Cys220', 'note': 'S lone pair / S–H σ*'},
        'pdb_native': '2OCJ', 'pdb_mutant': '2VUK',
        'ref': 'Joerger AC et al. (2006) PNAS 103:15056',
    },
    'TP53_C275F':  {
        'label':  'TP53 C275F (p.Cys275Phe)',
        'gene':   'TP53',
        'native': {'compound': 'methanethiol','residue': 'Cys275', 'note': 'S lone pair / S–H σ*'},
        'mutant': {'compound': 'toluene',     'residue': 'Phe275', 'note': 'Phe π/π* HOMO–LUMO'},
        'pdb_native': '2AC0', 'pdb_mutant': '2BIM',
        'ref': 'Joerger AC et al. (2006) PNAS 103:15056',
    },
    'KEAP1_G333C': {
        'label':  'KEAP1 G333C (p.Gly333Cys)',
        'gene':   'KEAP1',
        'native': {'compound': 'formamide',   'residue': 'Gly333', 'note': 'Backbone C=O π/π* (no sidechain)'},
        'mutant': {'compound': 'methanethiol','residue': 'Cys333', 'note': 'S lone pair / S–H σ*'},
        'pdb_native': '1U6D', 'pdb_mutant': None,
        'ref': 'Lo SC et al. (2006) EMBO J 25:3605',
    },
    'KEAP1_R320Q': {
        'label':  'KEAP1 R320Q (p.Arg320Gln)',
        'gene':   'KEAP1',
        'native': {'compound': 'guanidine',  'residue': 'Arg320', 'note': 'Guanidinium π system'},
        'mutant': {'compound': 'acetamide',  'residue': 'Gln320', 'note': 'Amide C=O π/π*'},
        'pdb_native': '1U6D', 'pdb_mutant': None,
        'ref': 'Tong KI et al. (2006) Biochemistry 45:6845',
    },
    'STK11_F354L': {
        'label':  'STK11 F354L (p.Phe354Leu)',
        'gene':   'STK11',
        'native': {'compound': 'toluene',    'residue': 'Phe354', 'note': 'Phe π/π* HOMO–LUMO (R-spine RS4)'},
        'mutant': {'compound': 'isobutane',  'residue': 'Leu354', 'note': 'C–H σ/σ* (aliphatic — no π)'},
        'pdb_native': '2WTK', 'pdb_mutant': None,
        'ref': 'Zeqiraj E et al. (2009) Science 326:1707',
    },
    'STK11_D194N': {
        'label':  'STK11 D194N (p.Asp194Asn)',
        'gene':   'STK11',
        'native': {'compound': 'acetic_acid','residue': 'Asp194', 'note': 'Carboxylate C=O π/π* (DFG Asp)'},
        'mutant': {'compound': 'acetamide',  'residue': 'Asn194', 'note': 'Amide C=O π/π*'},
        'pdb_native': '2WTK', 'pdb_mutant': None,
        'ref': 'Zeqiraj E et al. (2009) Science 326:1707',
    },
}


def geom_to_string(atom_list):
    return '\n'.join(f'{sym}  {x:.4f}  {y:.4f}  {z:.4f}' for sym, (x, y, z) in atom_list)


def run_casscf(name, atom_list, label='', verbose=1):
    """
    Run RHF → CASSCF(2,2) for a given molecule.
    Returns dict with energies and active-space integrals.
    """
    try:
        from pyscf import gto, scf, mcscf, ao2mo
    except ImportError:
        return {'error': 'PySCF not installed. Run: pip install pyscf'}

    mol = gto.Mole()
    mol.atom   = geom_to_string(atom_list)
    mol.basis  = 'sto-3g'
    mol.charge = 0
    mol.spin   = 0
    mol.verbose = verbose
    try:
        mol.build()
    except Exception as e:
        return {'error': f'mol.build() failed: {e}'}

    # RHF
    mf = scf.RHF(mol)
    mf.max_cycle = 300
    mf.conv_tol  = 1e-10
    try:
        e_rhf = mf.kernel()
    except Exception as e:
        return {'error': f'RHF failed: {e}'}

    if not mf.converged:
        return {'error': 'RHF did not converge'}

    # CASSCF(2,2)
    mc = mcscf.CASSCF(mf, ncas=2, nelecas=2)
    mc.conv_tol      = 1e-9
    mc.conv_tol_grad = 1e-5
    mc.max_cycle_macro = 150
    try:
        e_casscf = mc.kernel()[0]
    except Exception as e:
        return {'error': f'CASSCF failed: {e}'}

    # Active-space integrals
    h1e, ecore = mc.get_h1eff()
    h2e_compressed = mc.get_h2eff()
    h2e = ao2mo.restore(1, h2e_compressed, mc.ncas)

    return {
        'molecule':       name,
        'label':          label,
        'formula':        mol.atom,
        'n_atoms':        mol.natm,
        'n_electrons':    mol.nelectron,
        'n_basis':        mol.nao,
        'basis':          'STO-3G',
        'active_space':   'CAS(2e, 2o)',
        'jw_qubits':      4,
        'rhf_converged':  bool(mf.converged),
        'casscf_converged': bool(mc.converged),
        'e_rhf_Ha':       float(e_rhf),
        'e_casscf_Ha':    float(e_casscf),
        'correlation_Ha': float(e_casscf - e_rhf),
        'ecore_Ha':       float(ecore),
        'h1e':            h1e.tolist(),
        'h2e_J_aa':       float(h2e[0, 0, 0, 0]),
        'h2e_J_ab':       float(h2e[0, 0, 1, 1]),
        'h2e_K_ab':       float(h2e[0, 1, 1, 0]),
    }


def run_all_mutations(verbose=1):
    results = {}
    computed_cache = {}  # avoid re-computing same molecule twice

    print("\n" + "="*70)
    print("  SOLANGE — PySCF CASSCF(2,2) for all NSCLC mutations")
    print("  Basis: STO-3G   Active space: CAS(2e,2o)   JW qubits: 4")
    print("="*70)

    for mut_id, cfg in MUTATION_MODELS.items():
        print(f"\n{'─'*70}")
        print(f"  {cfg['label']}")
        print(f"  Ref: {cfg['ref']}")
        print(f"{'─'*70}")

        mut_result = {'mutation_id': mut_id, 'label': cfg['label'],
                      'gene': cfg['gene'], 'pdb_native': cfg['pdb_native'],
                      'pdb_mutant': cfg['pdb_mutant'], 'ref': cfg['ref']}

        for side in ('native', 'mutant'):
            compound = cfg[side]['compound']
            residue  = cfg[side]['residue']
            note     = cfg[side]['note']
            cache_key = compound

            print(f"\n  [{side.upper()}] {residue} — model: {compound}")
            print(f"  Active space: {note}")

            if cache_key in computed_cache:
                print(f"  ↳ Using cached result for {compound}")
                r = dict(computed_cache[cache_key])
                r['label'] = note
            else:
                r = run_casscf(compound, GEOM[compound], label=note, verbose=verbose)
                computed_cache[cache_key] = r

            if 'error' in r:
                print(f"  ERROR: {r['error']}")
            else:
                print(f"  RHF energy:     {r['e_rhf_Ha']:+.8f} Ha   (converged: {r['rhf_converged']})")
                print(f"  CASSCF energy:  {r['e_casscf_Ha']:+.8f} Ha   (converged: {r['casscf_converged']})")
                print(f"  Correlation:    {r['correlation_Ha']:+.8f} Ha")
                print(f"  h1e diagonal:   [{r['h1e'][0][0]:+.6f},  {r['h1e'][1][1]:+.6f}]")
                print(f"  J_ab = {r['h2e_J_ab']:+.6f} Ha    K_ab = {r['h2e_K_ab']:+.6f} Ha")

            mut_result[side] = {
                'residue':  residue,
                'compound': compound,
                'note':     note,
                'casscf':   r,
            }

        # ΔE native→mutant (CASSCF level)
        if ('error' not in mut_result.get('native', {}).get('casscf', {'error': ''})):
            if ('error' not in mut_result.get('mutant', {}).get('casscf', {'error': ''})):
                dE = (mut_result['mutant']['casscf']['e_casscf_Ha']
                    - mut_result['native']['casscf']['e_casscf_Ha'])
                mut_result['delta_E_casscf_Ha'] = dE
                print(f"\n  ΔE (mutant − native, CASSCF): {dE:+.8f} Ha")
                print(f"  Note: ΔE reflects change in model-compound electronic structure,")
                print(f"  not binding affinity. Interpret relative to active space character.")

        results[mut_id] = mut_result

    return results


def print_summary(results):
    print("\n\n" + "="*70)
    print("  SUMMARY — CASSCF(2,2) energies, STO-3G, all mutations")
    print("="*70)
    print(f"  {'Mutation':<18} {'Side':<8} {'Model':<14} {'CASSCF (Ha)':>16} {'Corr (Ha)':>13}")
    print(f"  {'-'*68}")
    for mut_id, r in results.items():
        for side in ('native', 'mutant'):
            if side not in r:
                continue
            casscf = r[side].get('casscf', {})
            if 'error' in casscf:
                print(f"  {mut_id:<18} {side:<8} {'ERROR':<14}")
                continue
            print(f"  {mut_id:<18} {side:<8} {r[side]['compound']:<14} "
                  f"{casscf['e_casscf_Ha']:>+16.8f} {casscf['correlation_Ha']:>+13.8f}")
        if 'delta_E_casscf_Ha' in r:
            print(f"  {'':>18} {'ΔE':>8} {'':14} {r['delta_E_casscf_Ha']:>+16.8f}")
        print()
    print("="*70)
    print("""
  INTERPRETATION NOTES:
  1. Energies are for model compounds (sidechain proxies), not full residues.
     Absolute values are NOT comparable across different molecules.
  2. Within one mutation (native vs mutant), ΔE reflects the change in
     electronic structure of the active orbital pair.
  3. Correlation energy = CASSCF − RHF. Larger |correlation| means stronger
     multi-reference character (more quantum, less classical).
  4. Leu (isobutane) has no π system — its CAS(2,2) targets C–H σ/σ*,
     which is near-zero correlation. This is physically correct.
  5. For VQE: decompose h1e + h2e (per mutation) into Pauli operators via
     Jordan-Wigner transform using openfermion or qiskit-nature.
""")


if __name__ == '__main__':
    results = run_all_mutations(verbose=0)
    print_summary(results)

    # Save full results
    out = {
        'title':        'SOLANGE — PySCF CASSCF(2,2) for all NSCLC mutations',
        'generated_at': datetime.utcnow().isoformat() + 'Z',
        'basis':        'STO-3G',
        'method':       'RHF → CASSCF(2,2)',
        'active_space': 'CAS(2e, 2o) — HOMO/LUMO pair',
        'jw_qubits':    4,
        'mutations':    results,
    }
    with open('all_mutations_casscf.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  Full results saved to: all_mutations_casscf.json\n")
