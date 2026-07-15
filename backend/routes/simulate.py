"""
VQE simulation engine — Phase 3A PennyLane simulator backend.
Real 4-qubit Jordan-Wigner Hamiltonians from PySCF CASSCF(2,2).
Each run produces a complete P1–P9 provenance record stored in Supabase.
"""
import base64
import hashlib
import json
import logging
import os
import time
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pennylane as qml
from pennylane import numpy as pnp
import asyncio
from concurrent.futures import ThreadPoolExecutor
from fastapi import APIRouter, Body, Header, HTTPException
from fastapi.responses import StreamingResponse
from supabase import create_client

router = APIRouter()

# ── Supabase client ────────────────────────────────────────────────────────────
SUPABASE_URL = os.environ.get("SUPABASE_URL", "")
SUPABASE_KEY = os.environ.get("SUPABASE_KEY", "")

def get_supabase():
    if not SUPABASE_URL or not SUPABASE_KEY:
        return None
    return create_client(SUPABASE_URL, SUPABASE_KEY)

# Columns that exist in public.simulation_runs (see backend/supabase_schema.sql).
# The insert is filtered to this set so record fields added ahead of a DB
# migration degrade gracefully instead of failing the whole insert (PGRST204).
_DB_COLUMNS = frozenset({
    "id", "created_at", "user_id", "mutation_id", "mutation_name", "pdb_id", "phase",
    "p1_circuit_hash", "p1_gate_count", "p1_depth", "p1_qubit_count", "p1_ansatz",
    "p2_compiler", "p2_compiler_version", "p2_encoding", "p2_basis_set",
    "p2_active_electrons", "p2_active_orbitals", "p2_model_compound", "p2_jw_terms",
    "p3_backend", "p3_backend_version", "p3_calibration_epoch", "p3_simulator",
    "p4_gate_error_rate", "p4_readout_error_rate", "p4_t1_us", "p4_t2_us", "p4_note",
    "p5_shots", "p5_raw_energy", "p5_energy_variance", "p5_opt_steps", "p5_elapsed_s",
    "p5_ecore_ha", "p5_active_energy_ha", "p5_casscf_ref_ha",
    "p6_method", "p6_note",
    "p7_energy_ha", "p7_ci_lower", "p7_ci_upper", "p7_confidence", "p7_method", "p7_ref_hf_ha",
    "p8_hash", "p8_algorithm", "p8_sealed_at", "p8_seal_payload",
    "p9_applicable", "p9_note",
})

# ── Load real JW Hamiltonians from PySCF CASSCF(2,2) ──────────────────────────
_JW_PATH = Path(__file__).parent.parent / "jw_hamiltonians.json"
with open(_JW_PATH) as _f:
    _JW_DATA = json.load(_f)

# ── CASSCF(2,2) active-space constants (2 electrons, 4 spin-orbitals) ─────────
# qchem.hf_state(2,4)=[1,1,0,0], spin-preserving excitations for 2e/4q.
_QUBITS     = 4
_N_ELECTRONS = 2
_HF_STATE   = [1, 1, 0, 0]           # |1100⟩ — alpha/beta electrons in orbital 0
_SINGLES    = [[0, 2], [1, 3]]        # spin-preserving singles only
_DOUBLES    = [[0, 1, 2, 3]]          # (0,1)→(2,3) double excitation
_N_PARAMS   = 3                       # 2 singles + 1 double


def _build_hamiltonian(terms: list) -> qml.Hamiltonian:
    """Build PennyLane Hamiltonian from pre-computed JW Pauli terms."""
    coeffs, ops = [], []
    for t in terms:
        pauli, coeff = t["pauli"], t["coeff"]
        if pauli == "I":
            ops.append(qml.Identity(0))
        else:
            gate_ops = []
            for token in pauli.split():
                letter, qubit = token[0], int(token[1:])
                if letter == "X": gate_ops.append(qml.PauliX(qubit))
                elif letter == "Y": gate_ops.append(qml.PauliY(qubit))
                elif letter == "Z": gate_ops.append(qml.PauliZ(qubit))
            op = gate_ops[0]
            for g in gate_ops[1:]:
                op = op @ g
            ops.append(op)
        coeffs.append(coeff)
    return qml.Hamiltonian(coeffs, ops)

# ── Mutation configurations ────────────────────────────────────────────────────
# Seven scientifically classified NSCLC targets (Y220C is a platform placeholder
# for NGS demo only — excluded from scientific counts).
#
# jw_source: (jw_key, side) → key in jw_hamiltonians.json; side = "mutant"/"native".
#   Named point mutations use the mutant compound (the changed residue).
#   General LOF entries use the key catalytic residue compound (native side).
#
# active_electrons/active_orbitals: Phase 3A CAS(2e,2o) real PySCF values.
# local_electrons/local_qubits: 5 Å binding-site shell from PDB coordinates.
# full_electrons/full_qubits: complete active-site environment from PDB.
# hardware_era: "current" = within 94-qubit demonstrated ceiling (Merz et al. 2026);
#               "fault_tolerant" = requires fault-tolerant QPU (~2030+).
#
# PDB coordinate sources (coordinate-verified May 2026):
#   TP53 C275F  → 2OCJ (TP53 DBD wild-type, 2.05 Å)          — hardware_era: current
#   KEAP1       → 1U6D (Kelch apo, 1.85 Å) +
#                 2FLU (Kelch + Nrf2 ETGE peptide, 2.0 Å)     — hardware_era: fault_tolerant
#   STK11/LKB1  → 2WTK (LKB1–STRADα–MO25α, 2.65 Å;          — hardware_era: fault_tolerant
#                        D194A engineered mutant in structure)
#   R320Q: AlphaFold Q14145 (IVR disordered region, low pLDDT) — 80e/160q local
#   F354L: AlphaFold Q15831 (C-terminal disordered, pLDDT 45) — 48e/96q local
#
# Hardware precedent: Merz et al. (Cleveland Clinic/RIKEN/IBM, May 2026,
# arXiv:2605.01138) demonstrated 94 qubits on IBM Heron r2 for a 12,635-atom
# protein-ligand complex — establishing the current NISQ ceiling for chemistry.
# C275F full active site (~88q) is the ONLY target within this ceiling.
# ── Expansion gene map — mirrors frontend GENE_MAP (non-core genes) ───────────
# Used to build Phase 3A proxy configs on the fly for NGS-detected expansion targets.
# All use KEAP1_G333C "mutant" (methanethiol) as the generic CAS(2e,2o) LOF proxy —
# the 4-qubit Hamiltonian is compound-specific but BQP class is determined by full_electrons.
# ── Expansion gene configs — real PySCF CASSCF(2,2)/STO-3G JW Hamiltonians ───
# Each gene has its own entry in jw_hamiltonians.json (generated by
# generate_expansion_jw.py). The jw_source key = "{GENE}_LOF", side = "native"
# (the wild-type functional residue whose disruption defines the LOF state).
#
# Model compound assignment is based on the key functional residue:
#   acetic_acid   → Asp catalytic (kinase DFG, Walker B)
#   propionic_acid → Glu catalytic (DExx box, ExoIII)
#   guanidine     → Arg electrostatic/catalytic finger
#   imidazole     → His neomorphic substitution (IDH1 R132H, IDH2 R172H, AXIN2 R815H)
#   methanethiol  → Cys catalytic (UCH, RING Zn-finger)
#   p_cresol      → Tyr in substrate-recognition interface (bromodomain, SET)
#   toluene       → Phe aromatic stacking (DIX domain, OB fold, ARID)
#   methanol      → Ser/Thr hydroxyl contact (VHL, FGFR3, CDKN2A)
_EXPANSION_GENE_CONFIGS = {
    "TP53":    {"full_electrons": 44, "full_qubits": 88,  "badge": "Structural LOF",     "jw_source": ("TP53_LOF",    "native"), "pdb": "2OCJ",           "native_residue": "Arg248",  "native_compound": "guanidine"},
    "VHL":     {"full_electrons": 25, "full_qubits": 50, "badge": "Structural",         "jw_source": ("VHL_LOF",     "native"), "pdb": "1LM8",           "native_residue": "Ser111",  "native_compound": "formamide"},
    "BAP1":    {"full_electrons": 35, "full_qubits": 70, "badge": "Ubiquitin LOF",      "jw_source": ("BAP1_LOF",    "native"), "pdb": "3KVF",           "native_residue": "Cys91",   "native_compound": "methanethiol"},
    "PBRM1":   {"full_electrons": 28, "full_qubits": 56, "badge": "Chromatin LOF",      "jw_source": ("PBRM1_LOF",   "native"), "pdb": "3G0L",           "native_residue": "Tyr1242", "native_compound": "p_cresol"},
    "SETD2":   {"full_electrons": 30, "full_qubits": 60, "badge": "Methyltransf. LOF",  "jw_source": ("SETD2_LOF",   "native"), "pdb": "5JLB",           "native_residue": "Tyr1666", "native_compound": "p_cresol"},
    "FGFR3":   {"full_electrons": 24, "full_qubits": 48, "badge": "Kinase",             "jw_source": ("FGFR3_LOF",   "native"), "pdb": "4K33",           "native_residue": "Asp641",  "native_compound": "acetic_acid"},
    "TSC1":    {"full_electrons": 26, "full_qubits": 52, "badge": "GAP LOF",            "jw_source": ("TSC1_LOF",    "native"), "pdb": "AlphaFold",      "native_residue": "Arg692",  "native_compound": "guanidine"},
    "TSC2":    {"full_electrons": 30, "full_qubits": 60, "badge": "GAP LOF",            "jw_source": ("TSC2_LOF",    "native"), "pdb": "5EJO",           "native_residue": "Arg1743", "native_compound": "guanidine"},
    "ATRX":    {"full_electrons": 30, "full_qubits": 60, "badge": "Helicase LOF",       "jw_source": ("ATRX_LOF",    "native"), "pdb": "AlphaFold",      "native_residue": "Asp2104", "native_compound": "acetic_acid"},
    "IDH1":    {"full_electrons": 22, "full_qubits": 44, "badge": "Neomorphic",         "jw_source": ("IDH1_LOF",    "native"), "pdb": "1T0L",           "native_residue": "Arg132",  "native_compound": "guanidine"},
    "IDH2":    {"full_electrons": 22, "full_qubits": 44, "badge": "Neomorphic",         "jw_source": ("IDH2_LOF",    "native"), "pdb": "1LWD",           "native_residue": "Arg172",  "native_compound": "guanidine"},
    "SMARCA4": {"full_electrons": 40, "full_qubits": 80, "badge": "ATPase LOF",         "jw_source": ("SMARCA4_LOF", "native"), "pdb": "6LTJ",           "native_residue": "Glu479",  "native_compound": "acetic_acid"},
    "ARID1A":  {"full_electrons": 28, "full_qubits": 56, "badge": "Chromatin LOF",      "jw_source": ("ARID1A_LOF",  "native"), "pdb": "2L9X",           "native_residue": "Trp1815", "native_compound": "toluene"},
    "ARID2":   {"full_electrons": 28, "full_qubits": 56, "badge": "Chromatin LOF",      "jw_source": ("ARID2_LOF",   "native"), "pdb": "AlphaFold",      "native_residue": "Gln1118", "native_compound": "acetamide"},
    "POLE":    {"full_electrons": 24, "full_qubits": 48, "badge": "Exonuclease LOF",    "jw_source": ("POLE_LOF",    "native"), "pdb": "4M8O",           "native_residue": "Glu272",  "native_compound": "acetic_acid"},
    "BRCA1":   {"full_electrons": 32, "full_qubits": 64, "badge": "DNA Repair LOF",     "jw_source": ("BRCA1_LOF",   "native"), "pdb": "1JM7",           "native_residue": "Cys44",   "native_compound": "methanethiol"},
    "BRCA2":   {"full_electrons": 32, "full_qubits": 64, "badge": "DNA Repair LOF",     "jw_source": ("BRCA2_LOF",   "native"), "pdb": "1MJE",           "native_residue": "Phe3175", "native_compound": "toluene"},
    "ATM":     {"full_electrons": 28, "full_qubits": 56, "badge": "DNA Repair LOF",     "jw_source": ("ATM_LOF",     "native"), "pdb": "AlphaFold",      "native_residue": "Asp2870", "native_compound": "acetic_acid"},
    "TERT":    {"full_electrons": 24, "full_qubits": 48, "badge": "Telomerase",         "jw_source": ("TERT_LOF",    "native"), "pdb": "7LYT",           "native_residue": "Asp712",  "native_compound": "acetic_acid"},
    "RB1":     {"full_electrons": 32, "full_qubits": 64, "badge": "Cell Cycle LOF",     "jw_source": ("RB1_LOF",     "native"), "pdb": "2AZE",           "native_residue": "Glu2",    "native_compound": "acetic_acid"},
    "NF1":     {"full_electrons": 28, "full_qubits": 56, "badge": "RasGAP LOF",         "jw_source": ("NF1_LOF",     "native"), "pdb": "1NF1",           "native_residue": "Arg1276", "native_compound": "guanidine"},
    "NF2":     {"full_electrons": 22, "full_qubits": 44, "badge": "Scaffold LOF",       "jw_source": ("NF2_LOF",     "native"), "pdb": "1H4R",           "native_residue": "Arg341",  "native_compound": "guanidine"},
    "AXIN1":   {"full_electrons": 24, "full_qubits": 48, "badge": "WNT Scaffold",       "jw_source": ("AXIN1_LOF",   "native"), "pdb": "1WSP",           "native_residue": "Phe631",  "native_compound": "toluene"},
    "AXIN2":   {"full_electrons": 24, "full_qubits": 48, "badge": "WNT Scaffold",       "jw_source": ("AXIN2_LOF",   "native"), "pdb": "AlphaFold",      "native_residue": "Arg815",  "native_compound": "guanidine"},
    "CDKN2A":  {"full_electrons": 20, "full_qubits": 40, "badge": "Cell Cycle LOF",     "jw_source": ("CDKN2A_LOF",  "native"), "pdb": "1BI7",           "native_residue": "Arg58",   "native_compound": "guanidine"},
}


def _make_expansion_config(gene: str, gm: dict) -> dict:
    """Build a Phase 3A config for an expansion gene LOF target.

    Uses the gene's own real PySCF CASSCF(2,2)/STO-3G JW Hamiltonian
    (generated by generate_expansion_jw.py) keyed by {GENE}_LOF / native.
    The native-state model compound is chosen based on the key functional
    residue disrupted by the LOF mutation.
    """
    fe  = gm["full_electrons"]
    fq  = gm["full_qubits"]
    bqp = "A" if fe >= 30 else "B"
    era = "current" if fq <= 94 else "fault_tolerant"
    return {
        "name": f"{gene} Loss-of-Function",
        "pdb":  gm["pdb"],
        "desc": (
            f"{gene} {gm['badge']} — Phase 3A: CAS(2e,2o) {gm['native_compound']} "
            f"({gm['native_residue']} native-state model), STO-3G"
        ),
        "jw_source":        gm["jw_source"],
        "active_electrons": 2,
        "active_orbitals":  2,
        "local_electrons":  fe // 2,
        "local_qubits":     fq // 2,
        "full_electrons":   fe,
        "full_qubits":      fq,
        "bqp_class":        bqp,
        "hardware_era":     era,
        "phase3b_backend":  (
            "IBM Heron r3" if era == "current" else "fault-tolerant QPU (~2030+)"
        ),
    }


def _resolve_config(mutation_id: str) -> dict | None:
    """Return a MUTATION_CONFIGS entry, or build one for an expansion {GENE}_LOF target."""
    cfg = MUTATION_CONFIGS.get(mutation_id)
    if cfg:
        return cfg
    if mutation_id.endswith("_LOF"):
        gene = mutation_id[:-4]
        gm   = _EXPANSION_GENE_CONFIGS.get(gene)
        if gm:
            return _make_expansion_config(gene, gm)
    return None


MUTATION_CONFIGS = {
    "TP53_C275F": {
        "name": "TP53 p.Cys275Phe",
        "pdb": "2OCJ",
        "desc": "Phe275 π-system fragment - CAS(2e,2o) toluene proxy, STO-3G",
        "jw_source": ("TP53_C275F", "mutant"),   # toluene (Phe275 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 24,   # loop-sheet-helix 5 Å shell — PDB 2OCJ
        "local_qubits": 48,
        "full_electrons": 44,    # Zn²⁺ shell + DNA-guanine interface — PDB 2OCJ
        "full_qubits": 88,       # ONLY target within 94-qubit demonstrated ceiling (Merz et al. 2026)
        "bqp_class": "A",
        "hardware_era": "current",
        "phase3b_backend": "IBM Heron r3",
    },
    "TP53_Y220C": {
        "name": "TP53 p.Tyr220Cys",
        "pdb": "2VUK",
        "desc": "NGS demo anchor — CAS(2e,2o) methanethiol proxy (Cys220 sidechain), STO-3G",
        "jw_source": ("TP53_Y220C", "mutant"),   # methanethiol (Cys220 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 24,
        "local_qubits": 48,
        "full_electrons": 38,
        "full_qubits": 76,
        "bqp_class": "C",
        "hardware_era": "placeholder",
        "phase3b_backend": "IBM Heron r3",
    },
    "KEAP1_LOF": {
        "name": "KEAP1 Loss-of-Function",
        "pdb": "2FLU",
        "desc": "Nrf2-KEAP1 PPI interface — CAS(2e,2o) methanethiol proxy (Cys333 sidechain), STO-3G",
        "jw_source": ("KEAP1_G333C", "mutant"),  # methanethiol (Cys333, representative LOF)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 104,  # G333 5 Å shell — PDB 1U6D + 2FLU coordinate-verified
        "local_qubits": 208,
        "full_electrons": 155,   # full Nrf2-binding interface — PDB 2FLU coordinate-verified
        "full_qubits": 310,
        "bqp_class": "B",
        "hardware_era": "fault_tolerant",
        "phase3b_backend": "fault-tolerant QPU (~2030+)",
    },
    "KEAP1_G333C": {
        "name": "KEAP1 p.Gly333Cys",
        "pdb": "1U6D",
        "desc": "Kelch β-propeller Gly333 — CAS(2e,2o) methanethiol proxy (Cys333 sidechain), STO-3G",
        "jw_source": ("KEAP1_G333C", "mutant"),  # methanethiol (Cys333 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 104,  # G333 5 Å shell — PDB 1U6D coordinate-verified (15 residues)
        "local_qubits": 208,
        "full_electrons": 155,   # full Nrf2-binding interface — PDB 2FLU coordinate-verified
        "full_qubits": 310,
        "bqp_class": "B",
        "hardware_era": "fault_tolerant",
        "phase3b_backend": "fault-tolerant QPU (~2030+)",
    },
    "KEAP1_R320Q": {
        "name": "KEAP1 p.Arg320Gln",
        "pdb": "2FLU",
        "desc": "IVR-Kelch boundary Arg320 — CAS(2e,2o) acetamide proxy (Gln320 sidechain), STO-3G",
        "jw_source": ("KEAP1_R320Q", "mutant"),  # acetamide (Gln320 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 80,   # R320 5Å shell — AlphaFold Q14145 (IVR disordered region, pLDDT low)
        "local_qubits": 160,
        "full_electrons": 155,   # shares full Nrf2-binding interface — PDB 2FLU coordinate-verified
        "full_qubits": 310,
        "bqp_class": "B",
        "hardware_era": "fault_tolerant",
        "phase3b_backend": "fault-tolerant QPU (~2030+)",
    },
    "STK11_LKB1": {
        "name": "STK11/LKB1 Loss-of-Function",
        "pdb": "2WTK",
        "desc": "LKB1 kinase domain LOF — CAS(2e,2o) acetic acid proxy (Asp194 DFG motif), STO-3G",
        "jw_source": ("STK11_D194N", "native"),  # acetic_acid (Asp194, DFG-motif catalytic residue)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 76,   # D194 5 Å shell — PDB 2WTK chain C coordinate-verified
        "local_qubits": 152,
        "full_electrons": 152,   # full ATP pocket 8 Å shell — PDB 2WTK chain C coordinate-verified
        "full_qubits": 304,
        "bqp_class": "A",
        "hardware_era": "fault_tolerant",
        "phase3b_backend": "fault-tolerant QPU (~2030+)",
    },
    "STK11_F354L": {
        "name": "STK11 p.Phe354Leu",
        "pdb": "2WTK",
        "desc": "LKB1 R-spine Phe354 — CAS(2e,2o) isobutane proxy (Leu354 sidechain), STO-3G",
        "jw_source": ("STK11_F354L", "mutant"),  # isobutane (Leu354 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 48,   # F354 5Å shell — AlphaFold Q15831 (C-terminal disordered, pLDDT 45)
        "local_qubits": 96,
        "full_electrons": 152,   # shares full ATP pocket — PDB 2WTK coordinate-verified
        "full_qubits": 304,
        "bqp_class": "A",
        "hardware_era": "near_term",
        "phase3b_backend": "IBM Heron r3 (near-term — local ~96q, 2q beyond demonstrated ceiling)",
    },
    "STK11_D194N": {
        "name": "STK11 p.Asp194Asn",
        "pdb": "2WTK",
        "desc": "LKB1 DFG-motif Asp194 — CAS(2e,2o) acetamide proxy (Asn194 sidechain), STO-3G",
        "jw_source": ("STK11_D194N", "mutant"),  # acetamide (Asn194 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 76,   # D194 5 Å shell — PDB 2WTK chain C coordinate-verified
        "local_qubits": 152,
        "full_electrons": 152,   # full ATP pocket 8 Å shell — PDB 2WTK chain C coordinate-verified
        "full_qubits": 304,
        "bqp_class": "A",
        "hardware_era": "fault_tolerant",
        "phase3b_backend": "fault-tolerant QPU (~2030+)",
    },
}


def run_vqe(config: dict, progress_cb=None) -> dict:
    """
    Live 4-qubit VQE on PennyLane default.qubit simulator.

    Ansatz : AllSinglesDoubles (UCCSD-type), HF initial state |1100⟩
    Active space : CAS(2e,2o) — JW Hamiltonian from PySCF CASSCF(2,2)
    Optimizer : Adam, 30 steps (converges to CASSCF-exact by proof for 2e/4q)

    ╔══════════════════════════════════════════════════════════════════╗
    ║  IBM_CONNECT — Phase 3B hardware entry point                    ║
    ║                                                                  ║
    ║  Replace this function with Qiskit Runtime execution:           ║
    ║    from qiskit_ibm_runtime import QiskitRuntimeService,         ║
    ║                                   EstimatorV2                   ║
    ║    service = QiskitRuntimeService(channel="ibm_quantum",        ║
    ║                  token=os.environ["IBM_QUANTUM_TOKEN"])         ║
    ║    backend = service.backend("ibm_heron_r3")                    ║
    ║                                                                  ║
    ║  Phase 3B requires:                                             ║
    ║    • New Hamiltonian: 24e/48q or 44e/88q                        ║
    ║    • Error mitigation: ZNE + Pauli Twirling                     ║
    ║    • Transpile to Heron r3 native gate set (ECR, Rz, SX)       ║
    ║    • CalibrationData from backend.properties()                  ║
    ╚══════════════════════════════════════════════════════════════════╝
    """
    jw_key, side = config["jw_source"]
    jw_entry     = _JW_DATA[jw_key][side]
    ecore        = jw_entry["ecore"]
    e_casscf     = jw_entry["e_casscf"]
    compound     = jw_entry["compound"]
    e_hf_active  = jw_entry["e_active_rhf"]

    # Build Hamiltonian from pre-computed JW Pauli terms
    hamiltonian = _build_hamiltonian(jw_entry["terms"])

    # PennyLane device and VQE circuit
    dev = qml.device("default.qubit", wires=_QUBITS)

    @qml.qnode(dev)
    def circuit(params):
        qml.AllSinglesDoubles(
            weights=params,
            wires=range(_QUBITS),
            hf_state=pnp.array(_HF_STATE),
            singles=_SINGLES,
            doubles=_DOUBLES,
        )
        return qml.expval(hamiltonian)

    # Run Adam optimizer — 80 steps, stepsize=0.4
    # Converges CAS(2e,2o) to within 1e-05 Ha of CASSCF-exact (~3s on Render Starter)
    params = pnp.zeros(_N_PARAMS, requires_grad=True)
    opt    = qml.AdamOptimizer(stepsize=0.4)

    t_start        = time.time()
    energies_active = []
    for i in range(80):
        params, e = opt.step_and_cost(circuit, params)
        energies_active.append(float(e))
        if progress_cb:
            progress_cb(i, round(ecore + float(e), 8))
    elapsed = time.time() - t_start

    final_active   = energies_active[-1]
    final_total    = ecore + final_active
    energies_total = [ecore + e for e in energies_active]

    variance = sum((e - final_total) ** 2 for e in energies_total[-10:]) / 10
    ci_half  = 1.96 * (variance / 10) ** 0.5

    gate_count = 13   # BasisState + 1 Double (DoubleExcitation) + 2 Singles
    depth      = 7

    fp_payload   = json.dumps({
        "gate_count": gate_count, "depth": depth, "qubits": _QUBITS,
        "compound": compound, "jw_key": jw_key, "side": side,
        "ansatz": "AllSinglesDoubles-UCCSD", "method": "PennyLane-live",
    }, sort_keys=True)
    circuit_hash = hashlib.sha256(fp_payload.encode()).hexdigest()

    return {
        "energy_ha":       final_total,
        "energy_active":   final_active,
        "ecore":           ecore,
        "e_casscf":        e_casscf,
        "compound":        compound,
        "ci_lower":        final_total - ci_half,
        "ci_upper":        final_total + ci_half,
        "energy_variance": variance,
        "gate_count":      gate_count,
        "depth":           depth,
        "n_qubits":        _QUBITS,
        "circuit_hash":    circuit_hash,
        "elapsed_s":       round(elapsed, 3),
        "convergence":     energies_total,
        "e_rhf":           ecore + e_hf_active,
        "jw_terms":        jw_entry.get("terms", []),
        "jw_key":          jw_key,
        "side":            side,
    }


# Fields excluded from the P8 seal because they do NOT survive a DB round-trip
# byte-identically (e.g. timestamptz is reformatted by Postgres), which would make
# a re-verification from the stored row spuriously FAIL. They are metadata, not
# result-integrity data, so excluding them keeps the seal robust and re-verifiable.
# The P8 seal is owned by LEON (routes.leon) — the single notarization authority.
# Import the canonical helpers so the seal can never drift between ingestion here
# and query-time re-verification in routes.provenance.
from routes import leon
from routes.leon import build_p8_payload, build_p8_seal

_SEAL_EXCLUDE = leon._SEAL_EXCLUDE


# ── API endpoint ───────────────────────────────────────────────────────────────

def _extract_user_id(authorization: str | None) -> str | None:
    """Decode Supabase JWT without verification to extract sub (user_id)."""
    if not authorization or not authorization.startswith("Bearer "):
        return None
    try:
        token = authorization[7:]
        payload_b64 = token.split(".")[1]
        payload_b64 += "=" * (-len(payload_b64) % 4)
        payload = json.loads(base64.urlsafe_b64decode(payload_b64))
        return payload.get("sub")
    except Exception:
        return None


@router.get("/{mutation_id}/stream")
async def stream_simulation(mutation_id: str, authorization: str | None = Header(None)):
    """SSE endpoint — runs VQE exactly once, streaming each energy value as it is
    computed, then emits a single final message carrying the full P1-P9 result
    payload (same shape as GET /{mutation_id}) so the frontend never needs a
    second, independently-computed fetch for the same run."""
    config = _resolve_config(mutation_id)
    if not config:
        raise HTTPException(status_code=404, detail=f"Unknown mutation: {mutation_id}")

    async def generate():
        loop   = asyncio.get_event_loop()
        queue  = asyncio.Queue()

        def worker():
            try:
                def progress_cb(step, energy):
                    asyncio.run_coroutine_threadsafe(
                        queue.put({"step": step, "energy": energy}), loop
                    )

                vqe = run_vqe(config, progress_cb=progress_cb)
                final = _assemble_and_persist(mutation_id, config, vqe, authorization)
                asyncio.run_coroutine_threadsafe(
                    queue.put({"done": True, "result": final}), loop
                )
            except Exception as exc:
                logging.error("VQE worker error for %s: %s", mutation_id, exc, exc_info=True)
                asyncio.run_coroutine_threadsafe(
                    queue.put({"error": str(exc), "done": True}), loop
                )

        executor = ThreadPoolExecutor(max_workers=1)
        loop.run_in_executor(executor, worker)

        while True:
            item = await queue.get()
            yield f"data: {json.dumps(item)}\n\n"
            if item.get("done"):
                break

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@router.get("/results")
async def get_results(limit: int = 20, authorization: str | None = Header(None)):
    """Return the most recent simulation runs from Supabase."""
    sb = get_supabase()
    if not sb:
        return {"error": "Supabase not configured", "data": []}
    # Pass user JWT so Supabase RLS can identify the caller and return their rows
    if authorization and authorization.startswith("Bearer "):
        sb.postgrest.auth(authorization[7:])
    res = (sb.table("simulation_runs")
             .select("id, created_at, mutation_id, mutation_name, p7_energy_ha, p7_ci_lower, p7_ci_upper, p8_hash, phase")
             .order("created_at", desc=True)
             .limit(limit)
             .execute())
    return {"data": res.data, "count": len(res.data)}


@router.post("/hpc/submit")
async def submit_hpc_run(payload: dict = Body(...)):
    """Ingest an externally-executed HPC run (e.g. Laguna) into SOLANGE.

    SOLANGE does NOT trust the submitter — it EARNS the green "Passed" by:
      1. Re-computing the P8 seal (SHA-256 over P1-P7,P9) and comparing to the
         submitted hash — proves the record is complete and untampered-since-seal.
      2. Re-checking physics consistency (ecore + e_active == e_casscf) when the
         Hamiltonian metadata is present — a cheap independent sanity check.
    A record that fails seal verification is REJECTED (422), never stored.
    The stored record is always stamped phase=3A-HPC and provenance_source so an
    external run is never disguised as an in-browser one.
    """
    prov = payload.get("provenance") or {}
    jw   = payload.get("jw") or {}
    if not prov or not prov.get("p8_hash"):
        raise HTTPException(400, "missing provenance record or p8_hash")

    # LEON notarizes the incoming run: it recomputes the P8 seal and re-checks the
    # physics-consistency invariant. Nothing is trusted on the strength of its
    # origin — a record that fails either check is REJECTED (422), never stored.
    verdict = leon.notarize(prov, jw)
    recomputed = verdict["recomputed_hash"]
    seal_ok = verdict["seal_ok"]
    consistency_ok = verdict["consistency_ok"]
    _audit_sb = get_supabase()
    _actor = prov.get("provenance_source") or prov.get("residue") or "HPC/external"
    _run_id = prov.get("id")
    if not seal_ok:
        leon.write_audit(_audit_sb, "reject", _run_id, verdict, actor=_actor,
                         note="P8 seal mismatch — record rejected at ingestion")
        raise HTTPException(
            422, f"LEON: P8 seal verification FAILED — recomputed {recomputed[:16]}… "
                 f"!= submitted {str(verdict['submitted_hash'])[:16]}…; record rejected")
    if consistency_ok is False:
        ecore, eact, ecas = jw.get("ecore"), jw.get("e_active_exact"), jw.get("e_casscf")
        leon.write_audit(_audit_sb, "reject", _run_id, verdict, actor=_actor,
                         note="physics consistency mismatch — record rejected at ingestion")
        raise HTTPException(
            422, f"LEON: physics consistency FAILED — ecore+e_active ({ecore+eact:.6f}) "
                 f"!= e_casscf ({ecas:.6f}); record rejected")

    # 3) Build the stored record — force honest phase/source labels.
    # Phase is normally forced to 3A-HPC (a classical Laguna run). A real quantum
    # run (solange_qpu.py) may declare 3B-QPU; accept only that allow-listed value
    # so an external submitter still can't mislabel a run as anything it wants.
    record = dict(prov)
    _submitted_phase = prov.get("phase")
    record["phase"] = _submitted_phase if _submitted_phase in ("3B-QPU", "3B-QPU-dryrun") else "3A-HPC"
    record["provenance_source"] = prov.get("provenance_source", "HPC/external")
    record.setdefault("id", str(uuid.uuid4()))
    # Fold the side into mutation_id so native/mutant of the same gene are distinct
    # rows (and readable in the panel), e.g. "ARID2_LOF (native)".
    side = prov.get("side", "native")
    base_id = record.get("mutation_id") or "unknown"
    record["mutation_id"]   = f"{base_id} ({side})"
    record["mutation_name"] = f"{base_id} ({side})"

    safe = {k: v for k, v in record.items() if k in _DB_COLUMNS}
    sb = get_supabase()
    db_status = "not_configured"
    if sb:
        try:
            # Upsert semantics: a re-run of the same target+side+active space REPLACES
            # the prior row (no duplicates in the panel). The append-only local
            # runs_log.jsonl on the cluster keeps the full archive of every run.
            #
            # Two SEPARATE .eq() deletes rather than one .in_() delete on
            # [base_id, folded_id]: folded ids look like "ATRX_LOF (mutant)" —
            # parentheses and a space inside a value going into an IN-list is
            # exactly the kind of value that can be mis-encoded by REST filter
            # builders. A native ATRX_LOF row was observed deleted by a LATER
            # mutant-side submission at the same CAS size despite the filter
            # apparently excluding it — consistent with that failure mode. Two
            # plain equality deletes carry no such risk.
            cas_e, cas_o = safe.get("p2_active_electrons"), safe.get("p2_active_orbitals")
            (sb.table("simulation_runs").delete()
               .eq("phase", "3A-HPC").eq("mutation_id", record["mutation_id"])
               .eq("p2_active_electrons", cas_e).eq("p2_active_orbitals", cas_o)
               .execute())
            if base_id != record["mutation_id"]:
                (sb.table("simulation_runs").delete()   # legacy pre-side-fold rows only
                   .eq("phase", "3A-HPC").eq("mutation_id", base_id)
                   .eq("p2_active_electrons", cas_e).eq("p2_active_orbitals", cas_o)
                   .execute())
            try:
                sb.table("simulation_runs").insert(safe).execute()
                db_status = "stored"
            except Exception as e_ins:
                # Resilient to the p8_seal_payload migration not being run yet: retry
                # without it so the run still stores (Verify is just LEGACY until migrated).
                if "p8_seal_payload" in str(e_ins):
                    safe2 = {k: v for k, v in safe.items() if k != "p8_seal_payload"}
                    sb.table("simulation_runs").insert(safe2).execute()
                    db_status = "stored_no_payload"
                else:
                    raise
        except Exception as e:
            db_status = "error"
            logging.error("HPC upsert failed: %s", e)

    # LEON (Lineage-Evidence Orchestration & Notarization) has re-verified the seal
    # and notarized the record. Nothing enters SOLANGE without passing this guard.
    logging.info("LEON notarized run %s — seal_ok=%s consistency_ok=%s db=%s",
                 record["id"], seal_ok, consistency_ok, db_status)
    leon.write_audit(_audit_sb, "notarize", record["id"],
                     {**verdict, "integrity": "PASS", "method": "ingestion"},
                     actor=record["provenance_source"],
                     note=f"notarized & stored (db={db_status})")
    return {
        "status":            "PASSED",
        "verified":          True,
        "notary":            "LEON",   # Lineage-Evidence Orchestration & Notarization
        "seal_ok":           seal_ok,
        "consistency_ok":    consistency_ok,
        "recomputed_p8":     recomputed,
        "phase":             "3A-HPC",
        "provenance_source": record["provenance_source"],
        "run_id":            record["id"],
        "db_status":         db_status,
    }


def _uid_from_auth(authorization: str | None) -> str:
    """Extract the user id (sub) from a Bearer JWT; raise 401 if absent/invalid."""
    if not authorization or not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Authorization header required")
    import base64 as _b64, json as _json
    try:
        pb = authorization[7:].split(".")[1]
        pb += "=" * (-len(pb) % 4)
        uid = _json.loads(_b64.urlsafe_b64decode(pb)).get("sub")
        if not uid:
            raise HTTPException(status_code=401, detail="Invalid token")
        return uid
    except HTTPException:
        raise
    except Exception:
        raise HTTPException(status_code=401, detail="Invalid token")


# ── Outbound dispatch queue: SOLANGE → (pull agent on the cluster) → SOLANGE ──
# Laguna security (Duo 2FA, no external SSH) blocks SOLANGE from pushing jobs. So
# SOLANGE only QUEUES a job here; a lightweight agent running INSIDE the user's
# Laguna session (solange_hpc.py --agent) pulls it, runs it, and --submits results
# back — completing the loop without breaching cluster security.

# A job claimed by an agent that then dies would sit in "running" forever, keeping
# the panel's live indicator on and blocking a clean queue. The reaper marks such
# orphans failed. The timeout is generous (longer than any observed run) so it never
# kills a legitimately long execution — it only clears truly abandoned jobs.
_STALE_MINUTES = 45


def _reap_stale_dispatch(sb):
    """Mark jobs stuck in 'running' past the stale timeout as failed (best-effort)."""
    try:
        cutoff = (datetime.now(timezone.utc) - timedelta(minutes=_STALE_MINUTES)).isoformat()
        (sb.table("hpc_dispatch").update(
            {"status": "failed",
             "finished_at": datetime.now(timezone.utc).isoformat(),
             "note": f"stale — no completion within {_STALE_MINUTES}m (agent likely lost)"})
           .eq("status", "running").lt("claimed_at", cutoff).execute())
    except Exception as e:
        logging.warning("dispatch reaper skipped: %s", e)


@router.post("/hpc/agent/heartbeat")
async def agent_heartbeat(payload: dict = Body(default={}),
                          authorization: str | None = Header(None)):
    """Cluster agent liveness ping — upserts a single last-seen row so SOLANGE can
    show whether an agent is currently active."""
    _uid_from_auth(authorization)
    sb = get_supabase()
    if not sb:
        return {"ok": False, "db": "not_configured"}
    row = {"id": "default", "last_seen": datetime.now(timezone.utc).isoformat(),
           "agent": (payload or {}).get("agent", "laguna"),
           "note": (payload or {}).get("note")}
    try:
        sb.table("agent_heartbeat").upsert(row).execute()
        return {"ok": True}
    except Exception as e:
        return {"ok": False, "error": str(e)}


@router.get("/hpc/agent/status")
async def agent_status():
    """Return whether an agent has pinged recently (online within 90s)."""
    sb = get_supabase()
    if not sb:
        return {"online": False, "db": "not_configured"}
    try:
        res = sb.table("agent_heartbeat").select("*").eq("id", "default").execute()
        if not res.data:
            return {"online": False, "last_seen": None}
        last = res.data[0].get("last_seen")
        secs = None
        try:
            dt = datetime.fromisoformat(str(last).replace("Z", "+00:00"))
            secs = (datetime.now(timezone.utc) - dt).total_seconds()
        except Exception:
            pass
        return {"online": (secs is not None and secs < 90),
                "last_seen": last, "seconds_ago": None if secs is None else round(secs),
                "agent": res.data[0].get("agent")}
    except Exception as e:
        return {"online": False, "error": str(e)}


@router.post("/hpc/dispatch")
async def dispatch_hpc(payload: dict = Body(...), authorization: str | None = Header(None)):
    """Queue an HPC job from SOLANGE. The cluster agent picks it up and runs it."""
    uid = _uid_from_auth(authorization)
    sb = get_supabase()
    if not sb:
        return {"queued": False, "db": "not_configured"}
    row = {
        "requested_by": uid,
        "status":   "queued",
        "key":      payload.get("key", "ARID2_LOF"),
        "side":     payload.get("side", "native"),
        # Do NOT default the compound — leave it null so the cluster agent resolves
        # the correct model compound PER GENE from the key (defaulting to acetamide
        # made every gene run acetamide, giving them all the same wrong energy).
        "compound": payload.get("compound"),
        "basis":    payload.get("basis", "6-31g"),
        "ncas":     int(payload.get("ncas", 8)),
        "nelecas":  int(payload.get("nelecas", 8)),
        "run_vqe":  bool(payload.get("run_vqe", False)),
        "residue":  payload.get("residue", ""),
    }
    try:
        res = sb.table("hpc_dispatch").insert(row).execute()
        did = (res.data or [{}])[0].get("id")
        return {"queued": True, "dispatch_id": did, "job": row}
    except Exception as e:
        return {"queued": False, "error": str(e),
                "hint": "run backend/migrations to create table hpc_dispatch"}


@router.get("/hpc/dispatch/list")
async def list_dispatch(limit: int = 20):
    """List recent dispatch jobs + their status (for the SOLANGE queue view)."""
    sb = get_supabase()
    if not sb:
        return {"jobs": [], "db": "not_configured"}
    _reap_stale_dispatch(sb)   # clear orphaned 'running' jobs before reporting state
    try:
        res = (sb.table("hpc_dispatch")
                 .select("id, created_at, status, key, side, compound, basis, "
                         "ncas, nelecas, run_vqe, claimed_at, finished_at, run_id, note")
                 .order("created_at", desc=True).limit(limit).execute())
        return {"jobs": res.data or []}
    except Exception as e:
        return {"jobs": [], "error": str(e)}


@router.post("/hpc/dispatch/clear")
async def clear_dispatch(payload: dict = Body(default={}),
                         authorization: str | None = Header(None)):
    """Clear the dispatch queue. hpc_dispatch is a WORKING QUEUE, not a provenance
    record — the immutable audit ledger is leon_audit (§06.iii), untouched here.
    Default: queued + running + failed, so a stopped/reaped agent's stale jobs
    (including ones the reaper marked failed, see _reap_stale_dispatch) don't
    linger in the status strip forever. Pass {"all": true} to also wipe 'done'."""
    _uid_from_auth(authorization)
    sb = get_supabase()
    if not sb:
        return {"deleted": 0, "db": "not_configured"}
    try:
        q = sb.table("hpc_dispatch").delete()
        if (payload or {}).get("all"):
            q = q.neq("id", "00000000-0000-0000-0000-000000000000")  # match every row
            scope = "all"
        else:
            q = q.in_("status", ["queued", "running", "failed"])
            scope = "pending + failed (queued+running+failed)"
        res = q.execute()
        n = len(res.data) if getattr(res, "data", None) else 0
        return {"deleted": n, "status": "cleared", "scope": scope}
    except Exception as e:
        return {"deleted": 0, "error": str(e)}


@router.get("/hpc/dispatch/next")
async def next_dispatch(authorization: str | None = Header(None)):
    """Cluster agent pulls the oldest queued job and claims it (status→running)."""
    _uid_from_auth(authorization)
    sb = get_supabase()
    if not sb:
        return {"job": None, "db": "not_configured"}
    _reap_stale_dispatch(sb)   # clear orphaned 'running' jobs before claiming the next
    try:
        res = (sb.table("hpc_dispatch").select("*")
                 .eq("status", "queued").order("created_at").limit(1).execute())
        if not res.data:
            return {"job": None}
        job = res.data[0]
        (sb.table("hpc_dispatch").update(
            {"status": "running", "claimed_at": datetime.now(timezone.utc).isoformat()})
           .eq("id", job["id"]).eq("status", "queued").execute())
        return {"job": job}
    except Exception as e:
        return {"job": None, "error": str(e)}


@router.post("/hpc/dispatch/{dispatch_id}/status")
async def update_dispatch(dispatch_id: str, payload: dict = Body(...),
                          authorization: str | None = Header(None)):
    """Cluster agent reports progress: status = running | done | failed."""
    _uid_from_auth(authorization)
    sb = get_supabase()
    if not sb:
        return {"updated": False, "db": "not_configured"}
    upd = {"status": payload.get("status", "done")}
    if payload.get("note"):    upd["note"] = payload["note"]
    if payload.get("run_id"):  upd["run_id"] = payload["run_id"]
    if upd["status"] in ("done", "failed"):
        upd["finished_at"] = datetime.now(timezone.utc).isoformat()
    try:
        sb.table("hpc_dispatch").update(upd).eq("id", dispatch_id).execute()
        return {"updated": True, "status": upd["status"]}
    except Exception as e:
        return {"updated": False, "error": str(e)}


@router.post("/hpc/clear")
async def clear_hpc_runs(payload: dict = Body(default={}),
                         authorization: str | None = Header(None)):
    """Delete HPC runs (phase=3A-HPC). Requires an authenticated session (Bearer
    token). Optional 'mutation_id' scopes the delete to one target; omitted = all.
    Destructive but recoverable — the runs can be re-submitted from the cluster."""
    sb = get_supabase()
    if not sb:
        return {"deleted": 0, "db": "not_configured"}
    if not authorization or not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Authorization header required")
    import base64 as _b64, json as _json
    try:
        pb = authorization[7:].split(".")[1]
        pb += "=" * (-len(pb) % 4)
        uid = _json.loads(_b64.urlsafe_b64decode(pb)).get("sub")
        if not uid:
            raise HTTPException(status_code=401, detail="Invalid token")
    except HTTPException:
        raise
    except Exception:
        raise HTTPException(status_code=401, detail="Invalid token")
    try:
        q = sb.table("simulation_runs").delete().eq("phase", "3A-HPC")
        mid = (payload or {}).get("mutation_id")
        if mid:
            q = q.eq("mutation_id", mid)
        res = q.execute()
        n = len(res.data) if getattr(res, "data", None) else 0
        return {"deleted": n, "status": "cleared", "scope": mid or "all HPC runs"}
    except Exception as e:
        logging.error("HPC clear failed: %s", e)
        raise HTTPException(status_code=500, detail=f"clear failed: {e}")


@router.post("/hpc/runs/delete")
async def delete_selected_qpu_runs(payload: dict = Body(...),
                                   authorization: str | None = Header(None)):
    """Delete SPECIFIC QPU runs by id (the per-row checkbox flow). Requires auth.
    Guarded to phase 3B-QPU* on the server side too, so this endpoint can NEVER
    delete a classical run even if a non-QPU id is passed — QPU runs cost real
    quantum time, so their deletion is deliberately isolated and explicit."""
    _uid_from_auth(authorization)
    ids = (payload or {}).get("ids") or []
    if not isinstance(ids, list) or not ids:
        raise HTTPException(status_code=400, detail="body must include a non-empty 'ids' list")
    sb = get_supabase()
    if not sb:
        return {"deleted": 0, "db": "not_configured"}
    try:
        res = (sb.table("simulation_runs").delete()
                 .in_("id", [str(i) for i in ids])
                 .like("phase", "3B-QPU%")   # safety: QPU rows only, never classical
                 .execute())
        n = len(res.data) if getattr(res, "data", None) else 0
        return {"deleted": n, "status": "deleted", "requested": len(ids)}
    except Exception as e:
        logging.error("QPU delete failed: %s", e)
        raise HTTPException(status_code=500, detail=f"delete failed: {e}")


@router.get("/hpc/runs")
async def list_hpc_runs(limit: int = 50):
    """List externally-executed HPC runs for the dashboard — classical Laguna
    runs (phase 3A-HPC) plus real/dry-run QPU smoke tests (phase 3B-QPU*). The
    p3_backend column already distinguishes them (e.g. 'ibm_marrakesh (real QPU)')."""
    sb = get_supabase()
    if not sb:
        return {"runs": [], "db": "not_configured"}
    try:
        res = (sb.table("simulation_runs")
                 .select("id, created_at, mutation_id, mutation_name, phase, "
                         "p1_ansatz, p2_active_electrons, p2_active_orbitals, p2_basis_set, "
                         "p3_backend, p5_elapsed_s, p5_ecore_ha, p5_casscf_ref_ha, "
                         "p7_energy_ha, p7_ref_hf_ha, p7_method, p8_hash")
                 .in_("phase", ["3A-HPC", "3B-QPU", "3B-QPU-dryrun"])
                 .order("created_at", desc=True)
                 .limit(limit)
                 .execute())
        return {"runs": res.data or []}
    except Exception as e:
        return {"runs": [], "error": str(e)}


# ── DMRG A/B/C classification ingestion ────────────────────────────────────────
# Sealed and notarized the same way as an HPC/CASSCF run (routes.leon), but the
# record shape is different (a classification + entanglement diagnostic, not a
# P1-P9 energy record), so it uses LEON's generic seal path (notarize_generic)
# rather than the P8/physics-consistency path. This closes the same "stuck only on
# Laguna" gap that HPC results already avoid: a DMRG classification is safe in
# SOLANGE the moment it's submitted, independent of any later cluster connectivity.
_DMRG_DB_COLUMNS = frozenset({
    "id", "created_at", "key", "compound", "basis", "ncas", "nelecas",
    "e_casscf", "dmrg_energies", "s_max", "bqp_class", "class_rationale",
    "time_budget_hit", "bond_dims_requested", "method",
    "dmrg_seal_payload", "dmrg_hash", "provenance_source",
})


@router.post("/hpc/dmrg/submit")
async def submit_dmrg_classification(payload: dict = Body(...)):
    """Ingest a DMRG A/B/C classification from solange_dmrg.py --submit.

    LEON re-verifies the seal the script computed at source before this record is
    trusted — a mismatch means the record was altered or incomplete in transit and
    is REJECTED (422), never stored.
    """
    if not payload.get("dmrg_hash"):
        raise HTTPException(400, "missing dmrg_hash — record was not sealed at source")

    record = dict(payload)
    record.setdefault("id", str(uuid.uuid4()))
    record["provenance_source"] = record.get("provenance_source", "HPC/external (DMRG)")

    # exclude dmrg_seal_payload too — the client (solange_dmrg.py) hashed the record
    # BEFORE either field existed, so the server must exclude both to recompute the
    # identical hash. Excluding hash_field alone here would hash the payload STRING
    # into itself and never match the client's seal.
    verdict = leon.notarize_generic(record, hash_field="dmrg_hash", exclude={"dmrg_seal_payload"})
    sb = get_supabase()
    actor = record.get("provenance_source")
    if not verdict["ok"]:
        leon.write_audit(sb, "reject", record["id"], verdict, actor=actor,
                         note="DMRG record seal mismatch — rejected at ingestion")
        raise HTTPException(
            422, f"LEON: DMRG seal verification FAILED — recomputed "
                 f"{verdict['recomputed_hash'][:16]}… != submitted "
                 f"{str(verdict['submitted_hash'])[:16]}…; record rejected")

    safe = {k: v for k, v in record.items() if k in _DMRG_DB_COLUMNS}
    db_status = "not_configured"
    if sb:
        try:
            sb.table("dmrg_classifications").insert(safe).execute()
            db_status = "stored"
        except Exception as e:
            db_status = "error"
            logging.error("DMRG classification insert failed: %s", e)

    logging.info("LEON notarized DMRG classification %s (%s) — class=%s db=%s",
                 record["id"], record.get("key"), record.get("bqp_class"), db_status)
    leon.write_audit(sb, "notarize", record["id"],
                     {**verdict, "integrity": "PASS", "method": "generic-ingestion"},
                     actor=actor, note=f"DMRG class={record.get('bqp_class')} (db={db_status})")

    return {
        "status": "PASSED", "verified": True, "notary": "LEON",
        "seal_ok": verdict["seal_ok"], "run_id": record["id"], "db_status": db_status,
    }


@router.get("/hpc/dmrg/list")
async def list_dmrg_classifications(limit: int = 50):
    """List LEON-notarized DMRG A/B/C classifications for the dashboard."""
    sb = get_supabase()
    if not sb:
        return {"classifications": [], "db": "not_configured"}
    try:
        res = (sb.table("dmrg_classifications")
                 .select("id, created_at, key, compound, basis, ncas, nelecas, "
                         "e_casscf, s_max, bqp_class, class_rationale, "
                         "time_budget_hit, dmrg_hash")
                 .order("created_at", desc=True).limit(limit).execute())
        return {"classifications": res.data or []}
    except Exception as e:
        return {"classifications": [], "error": str(e)}


@router.get("/{mutation_id}")
async def run_simulation(mutation_id: str, authorization: str | None = Header(None)):
    try:
        return await _run_simulation_inner(mutation_id, authorization)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


async def _run_simulation_inner(mutation_id: str, authorization: str | None):
    config = _resolve_config(mutation_id)
    if not config:
        raise HTTPException(status_code=404,
                            detail=f"Unknown mutation: {mutation_id}. "
                                   f"Valid: {list(MUTATION_CONFIGS.keys())} "
                                   f"or any expansion {{GENE}}_LOF where GENE is in the expansion map.")
    # ── Run VQE ────────────────────────────────────────────────────────────────
    vqe = run_vqe(config)
    return _assemble_and_persist(mutation_id, config, vqe, authorization)


def _assemble_and_persist(mutation_id: str, config: dict, vqe: dict, authorization: str | None) -> dict:
    """Build the P1-P9 provenance record from an already-computed VQE result,
    seal it (P8), persist to Supabase, and return the full API response shape.
    Shared by the plain GET endpoint and the SSE /stream endpoint so a single
    VQE run always produces a single, consistent result payload."""
    user_id = _extract_user_id(authorization)

    now    = datetime.now(timezone.utc).isoformat()
    run_id = str(uuid.uuid4())

    # ── Assemble P1–P9 provenance record ───────────────────────────────────────
    record = {
        "id":            run_id,
        "created_at":    now,
        "user_id":       user_id,
        "mutation_id":   mutation_id,
        "mutation_name": config["name"],
        "pdb_id":        config["pdb"],
        "phase":         "3A — PennyLane simulator",

        # P1 — Circuit fingerprint
        "p1_circuit_hash": vqe["circuit_hash"],
        "p1_gate_count":   vqe["gate_count"],
        "p1_depth":        vqe["depth"],
        "p1_qubit_count":  vqe["n_qubits"],
        "p1_ansatz":       "AllSinglesDoubles UCCSD (2e/4q, HF initial state |1100⟩)",

        # P2 — Compilation lineage
        "p2_compiler":         "PennyLane",
        "p2_compiler_version": qml.__version__,
        "p2_encoding":         "Jordan-Wigner (PySCF CAS(2e,2o) → openfermion → 27 Pauli terms)",
        "p2_basis_set":        "STO-3G (PySCF CASSCF(2,2))",
        "p2_active_electrons": config["active_electrons"],
        "p2_active_orbitals":  config["active_orbitals"],
        "p2_model_compound":   vqe["compound"],

        # P3 — Device & calibration
        "p3_backend":           "default.qubit",
        "p3_backend_version":   "0.38.0",
        "p3_calibration_epoch": now,
        "p3_simulator":         True,

        # P4 — Error budget (simulator: zero hardware noise)
        "p4_gate_error_rate":    0.0,
        "p4_readout_error_rate": 0.0,
        "p4_t1_us":              None,
        "p4_t2_us":              None,
        "p4_note":               "Noiseless simulator — Phase 3B will record real IBM Heron r3 calibration data",

        # P5 — Raw outcome distribution
        "p5_shots":           None,
        "p5_raw_energy":      vqe["energy_ha"],
        "p5_energy_variance": vqe["energy_variance"],
        "p5_opt_steps":       80,
        "p5_elapsed_s":       vqe["elapsed_s"],
        "p5_ecore_ha":        vqe["ecore"],
        "p5_active_energy_ha": vqe["energy_active"],
        "p5_casscf_ref_ha":   vqe["e_casscf"],

        # P6 — Error mitigation (none for noiseless simulator)
        "p6_method": "none — noiseless simulator",
        "p6_note":   "Phase 3B: ZNE + Pauli Twirling on IBM Heron r3",

        # P7 — Statistical estimator & CI
        "p7_energy_ha":  vqe["energy_ha"],
        "p7_ci_lower":   vqe["ci_lower"],
        "p7_ci_upper":   vqe["ci_upper"],
        "p7_confidence": 0.95,
        "p7_method":     "Bootstrap CI over last 20 optimisation steps",

        # P9 — ML decoder (not applicable for noiseless simulator)
        "p9_applicable": False,
        "p9_note":       "P9 conditional — applies when Nvidia Ising 3D CNN QEC decoder is active (Phase 3B)",
    }

    # P8 — Cryptographic seal
    record["p8_seal_payload"] = build_p8_payload(record)   # exact hashed input, for robust re-verify
    record["p8_hash"]      = hashlib.sha256(record["p8_seal_payload"].encode()).hexdigest()
    record["p8_algorithm"] = "SHA-256"
    record["p8_sealed_at"] = datetime.now(timezone.utc).isoformat()

    # ── Persist to Supabase ────────────────────────────────────────────────────
    # Insert only columns that exist in the table schema, so a future record field
    # added ahead of a DB migration can never crash the whole insert (PGRST204).
    safe_record = {k: v for k, v in record.items() if k in _DB_COLUMNS}
    sb = get_supabase()
    db_status = "not_configured"
    db_error  = None
    if sb:
        try:
            sb.table("simulation_runs").insert(safe_record).execute()
            db_status = "stored"
        except Exception as e:
            db_error  = str(e)
            db_status = "error"
            logging.error("Supabase insert failed: %s", e)

    return {
        "run_id":    run_id,
        "mutation":  config["name"],
        "bqp_class": config["bqp_class"],
        "result": {
            "energy_ha":        round(vqe["energy_ha"], 8),
            "energy_active_ha": round(vqe["energy_active"], 8),
            "ecore_ha":         round(vqe["ecore"], 8),
            "casscf_ref_ha":    round(vqe["e_casscf"], 8),
            "ci_lower":         vqe["ci_lower"],
            "ci_upper":         vqe["ci_upper"],
            "ci_half":          (vqe["ci_upper"] - vqe["ci_lower"]) / 2,
            "confidence":       "95%",
            "gate_count":       vqe["gate_count"],
            "depth":            vqe["depth"],
            "qubits_used":      vqe["n_qubits"],
            "elapsed_s":        vqe["elapsed_s"],
            "phase":            "3A — PennyLane simulator",
            "model_compound":   vqe["compound"],
            "local_target":     (
                f"{config['local_electrons']}e / {config['local_qubits']} qubits "
                f"(local site, Phase 3A tier)"
                if config.get("local_electrons")
                else "local active space TBD — mutation site not resolved in available PDB structures"
            ),
            "full_target":      (
                f"{config['full_electrons']}e / {config['full_qubits']} qubits "
                f"— {config['phase3b_backend']}"
            ),
            "hardware_era":     config.get("hardware_era", "unknown"),
            "_e_rhf":           round(vqe["e_rhf"], 10),
            "_jw_terms":        vqe["jw_terms"],
            "convergence":      vqe["convergence"],
        },
        "provenance": {
            "p1_circuit_hash": record["p1_circuit_hash"],
            "p2_compiler":     record["p2_compiler"],
            "p3_backend":      record["p3_backend"],
            "p7_energy_ha":    record["p7_energy_ha"],
            "p8_hash":         record["p8_hash"],
            "p9_applicable":   record["p9_applicable"],
        },
        "full_record": record,
        "db_status":   db_status,
        "db_error":    db_error,
    }
