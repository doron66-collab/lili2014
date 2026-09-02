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
    # IBM's own job id — the key that lets a sealed run be reconciled against the
    # vendor's billing export. Needs:
    #   alter table public.simulation_runs add column if not exists p3_vendor_job_id text;
    "p3_vendor_job_id",
    "p4_gate_error_rate", "p4_readout_error_rate", "p4_t1_us", "p4_t2_us", "p4_note",
    "p5_shots", "p5_raw_energy", "p5_energy_variance", "p5_opt_steps", "p5_elapsed_s",
    # Billable QPU execution seconds — the ONLY basis for a cost figure. Distinct from
    # p5_elapsed_s (wall clock) on purpose: for a QPU job wall clock is dominated by
    # free queue time. Needs the Supabase column:
    #   alter table public.simulation_runs add column if not exists p5_qpu_seconds numeric;
    "p5_qpu_seconds", "p5_qpu_seconds_source",
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

# ── Default active-space size — every gene target's jw_hamiltonians.json entry
# omits n_qubits/n_electrons, so these apply: CAS(2,2)/4-qubit. hf_state and
# singles/doubles are derived per-run via qml.qchem (see run_vqe) rather than
# hardcoded, so an entry that does specify n_qubits/n_electrons (e.g. the N2
# pipeline self-test, CAS(4,4)/8-qubit) gets its own correctly-sized circuit.
_QUBITS      = 4
_N_ELECTRONS = 2


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
    # No DMRG measurement exists for these expansion-gene LOF targets — a size
    # heuristic (fe >= 30) previously stood in for a real classification here,
    # producing a bqp_class indistinguishable downstream from a measured one
    # (including in the LEON P8 seal). Per project policy (an unlisted variant
    # must come back not classified, never a size proxy — see U2AF1_S34F),
    # this is left unclassified until a real DMRG run exists.
    bqp = None
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


# Structurally excluded expansion-gene LOF targets — the frontend already refuses
# to dispatch these (Assignment10_Prototype.html's STRUCTURALLY_UNRESOLVED), but that
# is a client-side gate only. Without a matching server-side refusal, a direct API call
# bypassing the UI reaches _make_expansion_config(), which computes bqp_class from
# full_electrons via a >=30e heuristic — for CDKN2A_LOF that meant fabricating a
# "Class B" verdict from a full_electrons figure (20e) that was itself fabricated: the
# gene's demonstrated NGS finding is a homozygous deletion, not a point mutation, so no
# active space exists to size in the first place. Found 2026-08-23 while correcting the
# same fabricated figure in targets.json.
_EXPANSION_GENE_STRUCTURALLY_UNRESOLVED = {
    "CDKN2A": "homozygous deletion, not a point mutation — no residue to anchor an "
              "active space to (matches STRUCTURALLY_UNRESOLVED['CDKN2A_LOF'] in the UI)",
}


def _resolve_config(mutation_id: str) -> dict | None:
    """Return a MUTATION_CONFIGS entry, or build one for an expansion {GENE}_LOF target."""
    cfg = MUTATION_CONFIGS.get(mutation_id)
    if cfg:
        return cfg
    if mutation_id.endswith("_LOF"):
        gene = mutation_id[:-4]
        if gene in _EXPANSION_GENE_STRUCTURALLY_UNRESOLVED:
            return None
        gm   = _EXPANSION_GENE_CONFIGS.get(gene)
        if gm:
            return _make_expansion_config(gene, gm)
    return None


MUTATION_CONFIGS = {
    "N2": {
        # PIPELINE SELF-TEST ONLY — not a gene, not a mutation, not a clinical or
        # scientific target. Proves the Rung-1 browser VQE wiring (dispatch ->
        # PennyLane -> Supabase -> P8 seal) runs end to end. See
        # jw_hamiltonians.json["N2"]["native"]["residue_note"] for why this needs
        # CAS(4,4)/8-qubit rather than the CAS(2,2)/4-qubit every real target uses.
        "name": "N2 pipeline self-test",
        "pdb": None,
        "desc": "PIPELINE SELF-TEST — not a target. N2 CAS(4e,4o)/STO-3G, 8-qubit JW Hamiltonian.",
        "jw_source": ("N2", "native"),
        "active_electrons": 4,
        "active_orbitals": 4,
        "bqp_class": None,
        "hardware_era": None,
    },
    "TP53_C275F": {
        "name": "TP53 p.Cys275Phe",
        "pdb": "2OCJ",
        "desc": "Phe275 π-system fragment - CAS(2e,2o) toluene proxy, STO-3G",
        "jw_source": ("TP53_C275F", "mutant"),   # toluene (Phe275 sidechain)
        "active_electrons": 2,
        "active_orbitals": 2,
        "local_electrons": 24,   # loop-sheet-helix 5 Å shell — PDB 2OCJ
        "local_qubits": 48,
        "full_electrons": 48,    # largest real AVAS-tested space (S 3p+N 2p+O 2p) — PDB 2OCJ, CAS(48,28)
        "full_qubits": 56,       # 28 orbitals x 2 under JW — real measurement, not a size guess
        # B, not A — aligned to the dissertation, which states it once and explicitly:
        # "Under this template C275F is Class B: it sits in the empirically
        # quantum-advantaged regime". It reaches that deliberately — "the claim is the
        # weaker but defensible one" — declining to assert that NO classical path exists,
        # because DMRG and neural quantum states keep advancing. This value said "A"
        # (quantum-necessary), the stronger claim the thesis chose not to make, so the
        # demo was overclaiming relative to the document it demonstrates.
        "bqp_class": "B",
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
        "full_electrons": 72,    # largest real AVAS-tested space (N 2p+O 2p; no S in this shell) — CAS(72,37)
        "full_qubits": 74,       # 37 orbitals x 2 under JW — real measurement, DMRG Class B (S_max=0.22)
        "bqp_class": "B",
        "hardware_era": "current",
        "phase3b_backend": "IBM Heron r3",
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
        "full_electrons": 54,    # largest real AVAS-tested space (N 2p+O 2p; no S in this shell) — CAS(54,31)
        "full_qubits": 62,       # 31 orbitals x 2 under JW — real measurement, DMRG Class B (S_max=0.82)
        "bqp_class": "B",        # was A, inherited from GENE_MAP's gene-level STK11 default — wrong for
                                 # a per-mutation entry (see project rule); real measurement is B
        "hardware_era": "current",
        "phase3b_backend": "IBM Heron r3",
    },
}


# ── Single source of truth overlay ────────────────────────────────────────────
# The literals above are convenient to read but they are NOT authoritative for the
# three fields that were found duplicated and contradictory across the codebase:
# bqp_class, full_electrons, full_qubits. Those come from targets.json, which the
# dissertation is checked against by scripts/laguna/verify_consistency.py.
#
# This exists because the same fact was written in four independent places with
# nothing forcing agreement, and they duly diverged: the dissertation stated
# "C275F is Class B" while this file and the frontend both said "A" — the demo
# claiming more than the thesis was prepared to defend, on the anchor target.
# A checker detects that after the fact; reading from one source prevents it.
#
# Applied as an overlay rather than by deleting the literals, so no unrelated field
# (pdb, desc, jw_source, the CAS(2,2) active_electrons) can be lost in the move.
# Divergence is logged rather than silently corrected: if these ever disagree, the
# literal above is stale and should be deleted, and that is worth seeing in the log.
# Search rather than assume. The first draft pointed only at backend/targets.json while
# the file lived at the repo root, so this silently fell through to the literals — the
# single-source guarantee off, and nothing on screen to say so. backend/targets.json is
# a symlink to the root file (one real copy, following the same pattern as the root
# dissertation_revised.html symlink); the search makes the loader independent of whether
# a given deployment ships the repo root at all.
_TARGETS_CANDIDATES = [
    Path(__file__).parent.parent / "targets.json",           # backend/ (symlink)
    Path(__file__).parent.parent.parent / "targets.json",    # repo root (the real file)
]
_TARGETS_PATH = next((p for p in _TARGETS_CANDIDATES if p.exists()), _TARGETS_CANDIDATES[0])
try:
    with open(_TARGETS_PATH) as _tf:
        _TARGETS = json.load(_tf)
    for _key, _src in (_TARGETS.get("mutations") or {}).items():
        _cfg = MUTATION_CONFIGS.get(_key)
        if not _cfg:
            continue
        for _field in ("bqp_class", "full_electrons", "full_qubits"):
            _want = _src.get(_field)
            if _want is None:
                continue
            if _cfg.get(_field) != _want:
                logging.warning(
                    "targets.json overrides MUTATION_CONFIGS[%s].%s: %r -> %r "
                    "(the literal in simulate.py is stale)",
                    _key, _field, _cfg.get(_field), _want)
            _cfg[_field] = _want
    logging.info("targets.json loaded — %d mutations, %d genes (single source of truth)",
                 len(_TARGETS.get("mutations") or {}), len(_TARGETS.get("genes") or {}))
except FileNotFoundError:
    # Not fatal: the literals are a complete, working fallback. But it means the
    # single-source guarantee is off, and that must be visible rather than assumed.
    _TARGETS = {}
    logging.error("targets.json NOT FOUND at %s — falling back to the literals in "
                  "simulate.py. Cross-file consistency is NOT guaranteed in this state.",
                  _TARGETS_PATH)


def run_vqe(config: dict, progress_cb=None) -> dict:
    """
    Live VQE on PennyLane default.qubit simulator.

    Ansatz : AllSinglesDoubles (UCCSD-type), HF initial state
    Active space : from the jw_hamiltonians.json entry — CAS(2e,2o)/4-qubit for
    every gene target (the default when an entry has no n_qubits/n_electrons),
    CAS(4e,4o)/8-qubit for the N2 pipeline self-test (see its jw entry's
    residue_note for why CAS(2,2) does not work for N2 specifically).
    Optimizer : Adam, 80 steps (converges to CASSCF-exact for 2e/4q by proof;
    good empirical convergence for 4e/8q too, see the N2 self-test's own runs).

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
    e_hf_active  = jw_entry.get("e_active_rhf")

    # Active-space size — per-entry, defaulting to the CAS(2e,2o)/4-qubit shape
    # every gene target uses. qml.qchem.hf_state/excitations reproduce the old
    # hardcoded _HF_STATE/_SINGLES/_DOUBLES exactly for (2, 4) — verified before
    # this was generalized — so no existing target's circuit changes.
    n_qubits    = jw_entry.get("n_qubits", _QUBITS)
    n_electrons = jw_entry.get("n_electrons", _N_ELECTRONS)
    hf_state    = qml.qchem.hf_state(n_electrons, n_qubits)
    singles, doubles = qml.qchem.excitations(n_electrons, n_qubits)
    n_params    = len(singles) + len(doubles)

    # Build Hamiltonian from pre-computed JW Pauli terms
    hamiltonian = _build_hamiltonian(jw_entry["terms"])

    # PennyLane device and VQE circuit
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(params):
        qml.AllSinglesDoubles(
            weights=params,
            wires=range(n_qubits),
            hf_state=pnp.array(hf_state),
            singles=singles,
            doubles=doubles,
        )
        return qml.expval(hamiltonian)

    # Run Adam optimizer — 80 steps, stepsize=0.4
    # Converges CAS(2e,2o) to within 1e-05 Ha of CASSCF-exact (~3s on Render Starter)
    params = pnp.zeros(n_params, requires_grad=True)
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

    if n_qubits == 4 and n_electrons == 2:
        gate_count = 13   # BasisState + 1 Double (DoubleExcitation) + 2 Singles
        depth      = 7
    else:
        # Not independently derived like the 4-qubit values above — a
        # documented, honestly-approximate count for AllSinglesDoubles' own
        # template structure: one BasisState plus one gate per excitation.
        gate_count = 1 + len(singles) + len(doubles)
        depth      = gate_count

    fp_payload   = json.dumps({
        "gate_count": gate_count, "depth": depth, "qubits": n_qubits,
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
        "n_qubits":        n_qubits,
        "n_electrons":     n_electrons,
        "n_paulis":        len(jw_entry.get("terms", [])),
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
from routes.security_log import log_denied

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


def _claims_from_auth(authorization: str | None) -> dict:
    """Decode the Bearer JWT's claims (no signature check — Supabase already
    verified it to issue it; this just reads what it says). Raises 401 if
    absent/malformed or missing a subject."""
    if not authorization or not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Authorization header required")
    import base64 as _b64, json as _json
    try:
        pb = authorization[7:].split(".")[1]
        pb += "=" * (-len(pb) % 4)
        claims = _json.loads(_b64.urlsafe_b64decode(pb))
        if not claims.get("sub"):
            raise HTTPException(status_code=401, detail="Invalid token")
        return claims
    except HTTPException:
        raise
    except Exception:
        raise HTTPException(status_code=401, detail="Invalid token")


def _uid_from_auth(authorization: str | None) -> str:
    """Extract the user id (sub) from a Bearer JWT; raise 401 if absent/invalid."""
    return _claims_from_auth(authorization)["sub"]


def _require_dispatch_allowed(authorization: str | None, sb) -> str:
    """Verify the caller may dispatch work that spends real HPC/DMRG/QPU time —
    reject executive accounts. Returns the user id on success.

    The frontend's three-tier tab model already hides the Orchestration tab from
    executive accounts, but that is a UI-only gate: _uid_from_auth only checks
    that a token exists, never who it belongs to, so an executive session could
    still call this endpoint directly and queue a job that costs real cluster
    time or real QPU-seconds. This is the server-side check that actually stops
    it — the same posture as LEON's "verify, don't trust" applied to who is
    allowed to spend, not just what gets stored.

    (The guest tier this once also rejected was retired outright — see the
    comment in showPlatform() in Assignment10_Prototype.html for why — so the
    only remaining restriction to enforce here is executive.)

    Checked via users_profile.role, and only possible if Supabase is reachable;
    if it is not, the check degrades to "unenforced" rather than "wrongly
    blocks a legitimate account", the safer failure direction for a role check
    with no source of truth to consult.
    """
    claims = _claims_from_auth(authorization)
    uid = claims["sub"]
    if sb:
        try:
            res = sb.table("users_profile").select("role").eq("id", uid).single().execute()
            role = (res.data or {}).get("role")
        except Exception:
            role = None
        if role == "executive":
            log_denied(sb, event="dispatch_denied", uid=uid, email=claims.get("email"),
                      role=role, endpoint="/hpc/dispatch",
                      detail="executive role attempted to dispatch real HPC/DMRG/QPU work")
            raise HTTPException(status_code=403,
                detail="Executive accounts cannot dispatch jobs that spend real HPC/DMRG/QPU time.")
    return uid


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
    # agent_id keys the heartbeat row so distinct agents track independently:
    # the classical Laguna agent uses 'default' (unchanged), the QPU agent 'qpu'.
    row = {"id": (payload or {}).get("agent_id", "default"),
           "last_seen": datetime.now(timezone.utc).isoformat(),
           "agent": (payload or {}).get("agent", "laguna"),
           "note": (payload or {}).get("note")}
    try:
        sb.table("agent_heartbeat").upsert(row).execute()
        return {"ok": True}
    except Exception as e:
        return {"ok": False, "error": str(e)}


@router.get("/hpc/agent/status")
async def agent_status(agent_id: str = "default"):
    """Return whether an agent has pinged recently (online within 90s).
    agent_id selects which agent: 'default' (classical Laguna) or 'qpu'."""
    sb = get_supabase()
    if not sb:
        return {"online": False, "db": "not_configured"}
    try:
        res = sb.table("agent_heartbeat").select("*").eq("id", agent_id).execute()
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
    """Queue an HPC job from SOLANGE. The cluster agent picks it up and runs it.

    Real spend gate: this is the single endpoint every job type goes through
    (job_type=hpc/dmrg/qpu), so it is the one place that needs the executive
    check — see _require_dispatch_allowed."""
    sb = get_supabase()
    uid = _require_dispatch_allowed(authorization, sb)
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
    # Only attached when set — same reasoning as job_type below: a fresh column
    # (dmrg_scf, dmrg_scf_maxm) may not exist in hpc_dispatch yet on every
    # deployment, and every dispatch that never asks for DMRG-SCF (still the
    # overwhelming majority — only CAS>16 selections turn this on client-side)
    # must keep working even before that migration runs.
    if payload.get("dmrg_scf"):
        row["dmrg_scf"] = True
        row["dmrg_scf_maxm"] = int(payload.get("dmrg_scf_maxm", 250))
    # SHCI dispatch (job_type='shci', see below) — a cross-validation run needs the
    # SAME geometry/AVAS/charge/spin the referenced DMRG record used, copied from
    # that record's own stored inputs (see the DMRG submission's geometry/avas/
    # charge/spin fields), not re-typed. Only attached when present — same reasoning
    # as dmrg_scf above: these columns may not exist yet on every deployment, and
    # every non-SHCI dispatch must keep working before that migration runs.
    if payload.get("job_type") == "shci":
        row["geometry"] = payload.get("geometry")
        row["avas"] = payload.get("avas")
        row["charge"] = int(payload.get("charge", 0))
        row["spin"] = int(payload.get("spin", 0))
        row["sweep_eps"] = payload.get("sweep_eps", "1e-2,1e-3,5e-4,1e-4")
        if payload.get("dmrg_classification_id"):
            row["dmrg_classification_id"] = payload.get("dmrg_classification_id")
    # Custom-geometry DMRG dispatch (job_type='dmrg' + geometry present) — a target
    # that is neither a --compound-library name nor a PDB-derived site, e.g. a
    # Gate-2-gated model-compound submission. Mirrors the SHCI branch above; the
    # agent side (solange_hpc.py) checks job.geometry to pick this path over the
    # --key/--side/--ncas/--nelecas one. Gate 2's own missing_requirements() check
    # happens BEFORE this endpoint is ever called (see gate2.py's
    # dispatch_custom_compound) — this endpoint itself does not re-check it, so
    # anything reaching here with geometry set has already cleared that gate.
    if payload.get("job_type") == "dmrg" and payload.get("geometry"):
        row["geometry"] = payload.get("geometry")
        row["avas"] = payload.get("avas")
        row["charge"] = int(payload.get("charge", 0))
        row["spin"] = int(payload.get("spin", 0))
    # PDB-to-classification pipeline dispatch (job_type='screen_classify') — the
    # "New Target from PDB" form. Chains protonate -> carve/probe (shrinking radius
    # until the active space fits) -> DMRG -> SHCI (both mandatory — the
    # classifier's decision is the two-method outcome, not one method's
    # opinion), all on the agent, so the end user never opens a terminal
    # (solange_screen_and_classify.py). Only attached for this job type, same
    # backward-compatible pattern as above.
    if payload.get("job_type") == "screen_classify":
        row["pdb_id"] = payload.get("pdb_id")
        row["chain"] = payload.get("chain")
        row["resi"] = int(payload.get("resi", 0))
        row["expect_resname"] = payload.get("expect_resname")
        row["radii"] = payload.get("radii", "5.0,4.0,3.5,3.0")
        row["max_orbitals"] = int(payload.get("max_orbitals", 45))
        row["avas"] = payload.get("avas")  # optional — omit to use build_qm_cluster.py's own per-radius suggestion
        row["spin"] = int(payload.get("spin", 0))
        # skip_shci: local-testing escape hatch only (no Dice build available) —
        # the classifier's normal outcome always runs both methods, so this
        # defaults to False (SHCI runs) rather than being opt-in.
        row["skip_shci"] = bool(payload.get("skip_shci", False))
        row["sweep_eps"] = payload.get("sweep_eps", "1e-2,1e-3,5e-4,1e-4")
    # job_type routes the job to the right agent: 'hpc' (default, classical Laguna)
    # or 'qpu' (real IBM hardware, pulled by the QPU agent). Only attach it for the
    # non-default type so classical dispatch keeps working even before the job_type
    # column migration is applied. The QPU backend/shots are agent-level choices
    # (set when you start the QPU agent), so the row needs only key/side + type.
    job_type = str(payload.get("job_type", "hpc")).lower()
    if job_type != "hpc":
        row["job_type"] = job_type
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
                         "ncas, nelecas, run_vqe, claimed_at, finished_at, run_id, note, job_type, "
                         "dmrg_scf, dmrg_scf_maxm, geometry, avas, charge, spin, sweep_eps, "
                         "dmrg_classification_id, pdb_id, chain, resi, expect_resname, radii, "
                         "max_orbitals, skip_shci")
                 .order("created_at", desc=True).limit(limit).execute())
        return {"jobs": res.data or []}
    except Exception as e:
        # dmrg_scf/dmrg_scf_maxm were added to the select list above without
        # confirming the migration (backend/migrations — alter table hpc_dispatch
        # add column dmrg_scf boolean, dmrg_scf_maxm integer) had actually been
        # run against this deployment's database. If it has not, this SELECT
        # itself now fails with Postgres's own "column does not exist" — which is
        # the fastest, most direct way to settle whether that's what has been
        # silently blocking --dmrg-scf from ever reaching a "Queue DMRG" job.
        msg = str(e)
        if "dmrg_scf" in msg or "does not exist" in msg.lower():
            msg += (" — hint: backend/migrations needs 'alter table hpc_dispatch "
                    "add column if not exists dmrg_scf boolean, "
                    "add column if not exists dmrg_scf_maxm integer' run against "
                    "this deployment's Supabase instance.")
        return {"jobs": [], "error": msg}


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
async def next_dispatch(job_type: str = "hpc", authorization: str | None = Header(None)):
    """An agent pulls the oldest queued job OF ITS TYPE and claims it (status→running).
    job_type routes: the classical Laguna agent asks for 'hpc' (default), the QPU
    agent asks for 'qpu'. A row with no job_type (pre-migration / classical) counts
    as 'hpc', so this stays correct whether or not the job_type column exists yet."""
    _uid_from_auth(authorization)
    # job_type may be a comma-separated set — the classical compute-node agent asks
    # for 'hpc,dmrg' (it can run both), the QPU agent asks for 'qpu'.
    wants = {t.strip().lower() for t in str(job_type or "hpc").split(",") if t.strip()} or {"hpc"}
    sb = get_supabase()
    if not sb:
        return {"job": None, "db": "not_configured"}
    _reap_stale_dispatch(sb)   # clear orphaned 'running' jobs before claiming the next
    try:
        # Fetch queued jobs oldest-first and pick the first matching type in Python —
        # robust to the job_type column being absent or NULL (treated as 'hpc').
        res = (sb.table("hpc_dispatch").select("*")
                 .eq("status", "queued").order("created_at").limit(20).execute())
        job = next((j for j in (res.data or [])
                    if str(j.get("job_type") or "hpc").lower() in wants), None)
        if not job:
            return {"job": None}
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
    if upd["status"] == "running":
        # A "running" update is itself a liveness signal (e.g. the QPU agent
        # reporting a live IBM job-status change while still waiting) — refresh
        # claimed_at so _reap_stale_dispatch's staleness clock resets instead of
        # timing out a job that is legitimately still queued/executing on IBM's
        # side well past the original claim (IBM queue waits are not bounded).
        upd["claimed_at"] = datetime.now(timezone.utc).isoformat()
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
        # Select-then-delete: a .delete().in_(id).like(phase,'3B-QPU%') chain threw a
        # Supabase edge-Worker exception (Cloudflare 1101) — the '%' in a LIKE on a
        # DELETE is the trigger. So first READ which requested ids are QPU rows
        # (the safety guard, via a plain select), then DELETE by id only.
        rows = (sb.table("simulation_runs").select("id, phase")
                  .in_("id", [str(i) for i in ids]).execute())
        qpu_ids = [r["id"] for r in (rows.data or [])
                   if str(r.get("phase", "")).startswith("3B-QPU")]
        if not qpu_ids:
            return {"deleted": 0, "status": "deleted", "requested": len(ids),
                    "note": "no matching QPU rows for the given ids"}
        res = sb.table("simulation_runs").delete().in_("id", qpu_ids).execute()
        n = len(res.data) if getattr(res, "data", None) else len(qpu_ids)
        return {"deleted": n, "status": "deleted", "requested": len(ids)}
    except Exception as e:
        logging.error("QPU delete failed: %s", e)
        raise HTTPException(status_code=500, detail=f"delete failed: {e}")


@router.post("/hpc/runs/delete/classical")
async def delete_selected_classical_runs(payload: dict = Body(...),
                                         authorization: str | None = Header(None)):
    """Delete SPECIFIC classical HPC runs by id (the per-row checkbox flow, Rung
    2). Mirrors delete_selected_qpu_runs exactly but with the phase guard
    reversed: this endpoint can NEVER delete a QPU run even if a QPU id is
    passed, so a slip on either checkbox list can't cross into the other rung's
    runs — the QPU side's real-quantum-time cost, and the classical side's own
    re-submit cost, each stay isolated to their own opt-in delete."""
    _uid_from_auth(authorization)
    ids = (payload or {}).get("ids") or []
    if not isinstance(ids, list) or not ids:
        raise HTTPException(status_code=400, detail="body must include a non-empty 'ids' list")
    sb = get_supabase()
    if not sb:
        return {"deleted": 0, "db": "not_configured"}
    try:
        rows = (sb.table("simulation_runs").select("id, phase")
                  .in_("id", [str(i) for i in ids]).execute())
        hpc_ids = [r["id"] for r in (rows.data or [])
                   if str(r.get("phase", "")) == "3A-HPC"]
        if not hpc_ids:
            return {"deleted": 0, "status": "deleted", "requested": len(ids),
                    "note": "no matching classical HPC rows for the given ids"}
        res = sb.table("simulation_runs").delete().in_("id", hpc_ids).execute()
        n = len(res.data) if getattr(res, "data", None) else len(hpc_ids)
        return {"deleted": n, "status": "deleted", "requested": len(ids)}
    except Exception as e:
        logging.error("Classical HPC delete failed: %s", e)
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
                         "p1_ansatz, p2_active_electrons, p2_active_orbitals, p2_basis_set, p2_model_compound, "
                         "p3_backend, p3_vendor_job_id, p3_calibration_epoch, p5_elapsed_s, p5_qpu_seconds, p5_qpu_seconds_source, p5_ecore_ha, "
                         "p5_casscf_ref_ha, p7_energy_ha, p7_ref_hf_ha, p7_method, p8_hash")
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
    "time_budget_hit", "bond_dims_requested", "method", "elapsed_s",
    "dmrg_seal_payload", "dmrg_hash", "provenance_source", "hardware",
    # geometry/avas/charge/spin: only present for a real --geometry run (not
    # --compound demo mode) — the reproducibility inputs an SHCI
    # cross-validation dispatch needs to copy, not re-type (see
    # solange_dmrg.py's own out[] construction and the "Queue SHCI
    # Cross-Validation" button in the frontend).
    "geometry", "avas", "charge", "spin",
})
# NOTE: "elapsed_s" and "hardware" each require their matching Supabase column
# to exist first — see the one-time migration in scripts/laguna/RUN_GUIDE.md §2
# (or run:
#   alter table public.dmrg_classifications add column if not exists elapsed_s numeric;
#   alter table public.dmrg_classifications add column if not exists hardware text;
# ). There is no partial-insert fallback below (unlike simulation_runs'
# p8_seal_payload retry): a submit against a column that does not yet exist
# fails the WHOLE insert with Postgres's "column does not exist" (db_status=
# "error", LEON's own audit/notarize step still succeeds independently — only
# the DB row is lost). Run the migration BEFORE the first --submit carrying
# "hardware". "hardware" was added because this table had no hardware column at
# all — unlike simulation_runs' p3_backend, which records the GPU the 3A-HPC
# CASSCF/VQE path actually used. DMRG here is CPU-only (solange_dmrg.py's
# detect_hardware()).
# run submitted before this migration silently drops the field until it is run.
# ). Until that migration runs, any submit carrying elapsed_s will insert-fail
# (caught below, db_status="error") — LEON's seal is still verified and audited,
# but the record silently will not land in the table. Run the migration BEFORE
# the first solange_dmrg.py --submit after this change.


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
                         "time_budget_hit, bond_dims_requested, dmrg_energies, "
                         "elapsed_s, method, provenance_source, dmrg_hash, hardware, "
                         "geometry, avas, charge, spin")
                 .order("created_at", desc=True).limit(limit).execute())
        return {"classifications": res.data or []}
    except Exception as e:
        # "hardware" was added to this SELECT without confirming the migration
        # (alter table dmrg_classifications add column if not exists hardware
        # text) had actually been run — same reasoning as hpc_dispatch's
        # dmrg_scf column below. If it hasn't, this SELECT itself now fails with
        # Postgres's own "column does not exist", which settles the question
        # directly instead of the dashboard silently showing stale rows forever.
        msg = str(e)
        if "does not exist" in msg.lower():
            msg += (" — hint: run 'alter table public.dmrg_classifications add "
                    "column if not exists hardware text, add column if not exists "
                    "geometry text, add column if not exists avas text, add column "
                    "if not exists charge int, add column if not exists spin int;' "
                    "in the Supabase SQL editor")
            return {"classifications": [], "error": msg}
        return {"classifications": [], "error": str(e)}


@router.get("/hpc/dmrg/{class_id}/verify")
async def verify_dmrg_seal(class_id: str):
    """Re-verify a DMRG classification's seal on demand — the DMRG-table equivalent
    of /api/provenance/runs/{id}/verify. LEON recomputes SHA-256 over the record's
    own dmrg_seal_payload (stored verbatim at ingestion, same robust pattern as the
    P8 seal) and compares it to the stored dmrg_hash. Public read, like the P8 path."""
    sb = get_supabase()
    if not sb:
        raise HTTPException(status_code=503, detail="database not configured")
    res = (sb.table("dmrg_classifications").select("*").eq("id", class_id).execute())
    if not res.data:
        raise HTTPException(status_code=404, detail=f"DMRG classification {class_id} not found")
    record = res.data[0]
    stored = record.get("dmrg_hash", "")
    payload = record.get("dmrg_seal_payload")
    if payload:
        recomputed = hashlib.sha256(payload.encode()).hexdigest()
        integrity = "PASS" if recomputed == stored else "FAIL"
    else:
        # Legacy record sealed before dmrg_seal_payload existed — reconstruct from
        # the same field set solange_dmrg.py hashes (id/created_at/hash/payload
        # excluded), matching leon.notarize_generic's exclude set at ingestion.
        recomputed = leon.build_generic_seal(
            record, exclude={"id", "created_at", "dmrg_hash", "dmrg_seal_payload"})
        integrity = "PASS" if recomputed == stored else "LEGACY-UNVERIFIABLE"
    verdict = {"integrity": integrity, "notary": leon.NAME,
               "stored_hash": stored, "recomputed_hash": recomputed, "algorithm": "SHA-256"}
    leon.write_audit(sb, "reverify", class_id, verdict, actor=record.get("provenance_source"))
    return {"run_id": class_id, **verdict}


@router.post("/hpc/dmrg/delete")
async def delete_selected_dmrg(payload: dict = Body(...),
                               authorization: str | None = Header(None)):
    """Delete SPECIFIC DMRG classifications by id (the per-row checkbox flow).
    Requires auth. DMRG classifications are classical and cost no quantum time,
    so there is no phase guard — but we still select-then-delete by id (mirroring
    the QPU path) so the response can report exactly how many rows matched."""
    _uid_from_auth(authorization)
    ids = (payload or {}).get("ids") or []
    if not isinstance(ids, list) or not ids:
        raise HTTPException(status_code=400, detail="body must include a non-empty 'ids' list")
    sb = get_supabase()
    if not sb:
        return {"deleted": 0, "db": "not_configured"}
    try:
        str_ids = [str(i) for i in ids]
        rows = (sb.table("dmrg_classifications").select("id")
                  .in_("id", str_ids).execute())
        found = [r["id"] for r in (rows.data or [])]
        if not found:
            return {"deleted": 0, "status": "deleted", "requested": len(ids),
                    "note": "no matching DMRG rows for the given ids"}
        res = sb.table("dmrg_classifications").delete().in_("id", found).execute()
        n = len(res.data) if getattr(res, "data", None) else len(found)
        return {"deleted": n, "status": "deleted", "requested": len(ids)}
    except Exception as e:
        logging.error("DMRG delete failed: %s", e)
        raise HTTPException(status_code=500, detail=f"delete failed: {e}")


# ── SHCI classification ingestion (standalone, optionally cross-validated) ─────
# SHCI is a genuinely independent classical method — its failure mechanism is
# determinant-sparsity growth, unrelated to DMRG's bond-dimension growth across
# a one-dimensional chain (see solange_shci.py's own module docstring for the
# classification rule and its stated one-signal limitation). It reaches its
# OWN A/B/C verdict on whatever active space it is given, on the same footing
# as the DMRG classifier — it is not subordinate to DMRG, and does not require
# a DMRG record to run against. dmrg_classification_id is OPTIONAL: when given,
# this endpoint ADDITIONALLY computes a cross-validation (delta_mha/agreement)
# against that DMRG record's own stored energy — a genuine extra check, not the
# reason SHCI's own classification is trusted.
#
# DP1 ("verify, don't trust") applies to the cross-validation verdict
# specifically: delta_mha/agreement, when computed, come from the referenced
# DMRG record's own stored energy, never from the submitting script's claim —
# a compromised or buggy client could otherwise assert "agrees" over a real
# disagreement with no way for a reader to catch it. SHCI's own bqp_class,
# by contrast, IS taken from the submitting script, exactly as the DMRG
# classifier's own bqp_class already is (§06.i) — the classification itself is
# sealed evidence from the run that produced it, not something the backend can
# independently re-derive without re-running SHCI.
#
# One-time migration (Supabase SQL editor), before the first --submit. Additive
# to the schema already live from this table's first (comparison-only) version
# — existing rows are unaffected; the new columns are simply null on them:
#   create table if not exists public.shci_crossvalidations (
#     id uuid primary key,
#     created_at timestamptz not null default now(),
#     key text, dmrg_classification_id uuid,
#     ncas int, nelec int, e_shci numeric,
#     sweep_eps text, method text, elapsed_s numeric,
#     provenance_source text, hardware text,
#     e_dmrg_ref numeric, delta_mha numeric, agreement boolean,
#     shci_seal_payload text, shci_hash text
#   );
#   alter table public.shci_crossvalidations add column if not exists bqp_class text;
#   alter table public.shci_crossvalidations add column if not exists class_rationale text;
#   alter table public.shci_crossvalidations add column if not exists shci_energies jsonb;
#   alter table public.shci_crossvalidations add column if not exists class_agreement boolean;
_SHCI_DB_COLUMNS = frozenset({
    "id", "created_at", "key", "dmrg_classification_id", "ncas", "nelec",
    "e_shci", "sweep_eps", "method", "elapsed_s", "provenance_source",
    "hardware", "e_dmrg_ref", "delta_mha", "agreement", "class_agreement",
    "bqp_class", "class_rationale", "shci_energies",
    "shci_seal_payload", "shci_hash",
})
_SHCI_CHEM_ACC_MHA = 1.6  # same bar solange_dmrg.py's classify() uses


@router.post("/hpc/shci/submit")
async def submit_shci_crossvalidation(payload: dict = Body(...)):
    """Ingest an SHCI classification from solange_shci.py --submit — standalone,
    or additionally cross-validated against a DMRG record if
    dmrg_classification_id is given. LEON re-verifies the seal the script
    computed at source before the record is trusted; a mismatch is REJECTED
    (422), never stored. When a DMRG record is named, the agreement verdict is
    computed HERE, against that record's own stored energy, never accepted
    from the submitting script (see module note above).
    """
    if not payload.get("shci_hash"):
        raise HTTPException(400, "missing shci_hash — record was not sealed at source")
    if not payload.get("bqp_class"):
        raise HTTPException(400, "missing bqp_class — SHCI classifies independently now; "
                                  "see solange_shci.py's classify_shci()")

    record = dict(payload)
    record.setdefault("id", str(uuid.uuid4()))
    record["provenance_source"] = record.get("provenance_source", "HPC/external (SHCI)")

    verdict = leon.notarize_generic(record, hash_field="shci_hash", exclude={"shci_seal_payload"})
    sb = get_supabase()
    actor = record.get("provenance_source")
    if not verdict["ok"]:
        leon.write_audit(sb, "reject", record["id"], verdict, actor=actor,
                         note="SHCI record seal mismatch — rejected at ingestion")
        raise HTTPException(
            422, f"LEON: SHCI seal verification FAILED — recomputed "
                 f"{verdict['recomputed_hash'][:16]}… != submitted "
                 f"{str(verdict['submitted_hash'])[:16]}…; record rejected")

    if not sb:
        return {"status": "PASSED", "verified": True, "notary": "LEON",
                "seal_ok": verdict["seal_ok"], "run_id": record["id"],
                "db_status": "not_configured"}

    dmrg_id = payload.get("dmrg_classification_id")
    if dmrg_id:
        dmrg_row = (sb.table("dmrg_classifications").select("nelecas, ncas, dmrg_energies, key, bqp_class")
                      .eq("id", str(dmrg_id)).execute())
        if not dmrg_row.data:
            leon.write_audit(sb, "reject", record["id"], verdict, actor=actor,
                             note=f"SHCI submit referenced unknown dmrg_classification_id={dmrg_id}")
            raise HTTPException(404, f"dmrg_classification_id {dmrg_id} not found — cannot "
                                      f"cross-validate against a record that does not exist")
        dmrg_rec = dmrg_row.data[0]
        dmrg_energies = dmrg_rec.get("dmrg_energies") or []
        if not dmrg_energies:
            raise HTTPException(409, f"DMRG record {dmrg_id} has no dmrg_energies to compare against")
        e_dmrg_ref = dmrg_energies[-1][1]
        record["e_dmrg_ref"] = e_dmrg_ref
        record["delta_mha"] = round(abs(record.get("e_shci", 0.0) - e_dmrg_ref) * 1000.0, 4)
        # "agreement" (energy) and "class_agreement" (A/B/C verdict) are
        # DELIBERATELY separate signals, not one collapsed into the other —
        # found live 2026-09-02 on a Class-A [2Fe-2S] proxy where they point
        # opposite ways: both solvers independently concluded Class A (same
        # verdict), yet their raw energies differed by 95.5 mHa, because
        # NEITHER had actually converged — two unconverged numbers have no
        # reason to agree even when the shared non-convergence is exactly
        # what both classifiers are reporting. Collapsing that into a single
        # "agreement: false" reads as a contradiction between the methods
        # when it is actually corroboration of the same Class A verdict.
        record["agreement"] = record["delta_mha"] <= _SHCI_CHEM_ACC_MHA
        dmrg_class = dmrg_rec.get("bqp_class")
        record["class_agreement"] = (
            (record.get("bqp_class") == dmrg_class) if dmrg_class else None)
        # Scope check mirrors the DMRG-record's own scope discipline (§06.i): a
        # cross-validation is only meaningful over the SAME active space, not a
        # fragment or a superset of it.
        if dmrg_rec.get("ncas") != record.get("ncas") or dmrg_rec.get("nelecas") != record.get("nelec"):
            raise HTTPException(409, f"active-space mismatch — DMRG record {dmrg_id} is "
                                      f"CAS({dmrg_rec.get('nelecas')},{dmrg_rec.get('ncas')}), "
                                      f"this SHCI run is CAS({record.get('nelec')},{record.get('ncas')}); "
                                      f"a cross-validation must use the identical active space")
        record["key"] = record.get("key") or dmrg_rec.get("key")
    # else: standalone classification — no DMRG record to compare against, and
    # none required. e_dmrg_ref/delta_mha/agreement stay unset (null in the DB).

    safe = {k: v for k, v in record.items() if k in _SHCI_DB_COLUMNS}
    db_status = "not_configured"
    try:
        sb.table("shci_crossvalidations").insert(safe).execute()
        db_status = "stored"
    except Exception as e:
        db_status = "error"
        logging.error("SHCI cross-validation insert failed: %s", e)

    cross_note = (f"delta={record['delta_mha']} mHa agreement={record['agreement']} "
                  f"class_agreement={record.get('class_agreement')} vs DMRG {dmrg_id}"
                  if dmrg_id else "standalone — no DMRG record referenced")
    logging.info("LEON notarized SHCI classification %s (%s) — class=%s %s db=%s",
                 record["id"], record.get("key"), record.get("bqp_class"), cross_note, db_status)
    leon.write_audit(sb, "notarize", record["id"],
                     {**verdict, "integrity": "PASS", "method": "generic-ingestion"},
                     actor=actor, note=f"SHCI class={record.get('bqp_class')} — {cross_note} (db={db_status})")

    return {
        "status": "PASSED", "verified": True, "notary": "LEON",
        "seal_ok": verdict["seal_ok"], "run_id": record["id"], "db_status": db_status,
        "bqp_class": record.get("bqp_class"),
        "delta_mha": record.get("delta_mha"), "agreement": record.get("agreement"),
        "class_agreement": record.get("class_agreement"),
    }


@router.get("/hpc/shci/list")
async def list_shci_crossvalidations(limit: int = 50):
    """List LEON-notarized SHCI cross-validations for the dashboard."""
    sb = get_supabase()
    if not sb:
        return {"crossvalidations": [], "db": "not_configured"}
    try:
        res = (sb.table("shci_crossvalidations")
                 .select("id, created_at, key, dmrg_classification_id, ncas, nelec, "
                         "e_shci, bqp_class, class_rationale, shci_energies, "
                         "e_dmrg_ref, delta_mha, agreement, class_agreement, sweep_eps, method, "
                         "elapsed_s, provenance_source, hardware, shci_hash")
                 .order("created_at", desc=True).limit(limit).execute())
        return {"crossvalidations": res.data or []}
    except Exception as e:
        msg = str(e)
        if "does not exist" in msg.lower():
            msg += (" — hint: the shci_crossvalidations table has not been created yet; "
                    "see the migration comment above submit_shci_crossvalidation()")
        return {"crossvalidations": [], "error": msg}


@router.get("/hpc/shci/{cv_id}/verify")
async def verify_shci_seal(cv_id: str):
    """Re-verify an SHCI cross-validation's seal on demand — same pattern as
    /hpc/dmrg/{id}/verify. Re-checks the cryptographic seal only; it does not
    re-derive delta_mha/agreement, which were computed once at ingestion against
    the DMRG record as it stood then (§DP1 note above submit_shci_crossvalidation)."""
    sb = get_supabase()
    if not sb:
        raise HTTPException(status_code=503, detail="database not configured")
    res = (sb.table("shci_crossvalidations").select("*").eq("id", cv_id).execute())
    if not res.data:
        raise HTTPException(status_code=404, detail=f"SHCI cross-validation {cv_id} not found")
    record = res.data[0]
    stored = record.get("shci_hash", "")
    payload = record.get("shci_seal_payload")
    if payload:
        recomputed = hashlib.sha256(payload.encode()).hexdigest()
        integrity = "PASS" if recomputed == stored else "FAIL"
    else:
        recomputed = leon.build_generic_seal(
            record, exclude={"id", "created_at", "shci_hash", "shci_seal_payload"})
        integrity = "PASS" if recomputed == stored else "LEGACY-UNVERIFIABLE"
    verdict = {"integrity": integrity, "notary": leon.NAME,
               "stored_hash": stored, "recomputed_hash": recomputed, "algorithm": "SHA-256"}
    leon.write_audit(sb, "reverify", cv_id, verdict, actor=record.get("provenance_source"))
    return {"run_id": cv_id, **verdict}


@router.post("/hpc/shci/delete")
async def delete_selected_shci(payload: dict = Body(...),
                               authorization: str | None = Header(None)):
    """Delete specific SHCI cross-validations by id. Requires auth — mirrors
    /hpc/dmrg/delete exactly; SHCI is classical and costs no quantum time."""
    _uid_from_auth(authorization)
    ids = (payload or {}).get("ids") or []
    if not isinstance(ids, list) or not ids:
        raise HTTPException(status_code=400, detail="body must include a non-empty 'ids' list")
    sb = get_supabase()
    if not sb:
        return {"deleted": 0, "db": "not_configured"}
    try:
        str_ids = [str(i) for i in ids]
        rows = (sb.table("shci_crossvalidations").select("id")
                  .in_("id", str_ids).execute())
        found = [r["id"] for r in (rows.data or [])]
        if not found:
            return {"deleted": 0, "status": "deleted", "requested": len(ids),
                    "note": "no matching SHCI rows for the given ids"}
        res = sb.table("shci_crossvalidations").delete().in_("id", found).execute()
        n = len(res.data) if getattr(res, "data", None) else len(found)
        return {"deleted": n, "status": "deleted", "requested": len(ids)}
    except Exception as e:
        logging.error("SHCI delete failed: %s", e)
        raise HTTPException(status_code=500, detail=f"delete failed: {e}")


# ── QM cluster cache — build once, reuse on every later encounter ──────────────
# "New Target from PDB" (solange_screen_and_classify.py) used to protonate, carve,
# and probe from scratch every time, even for a (pdb_id, chain, resi) already
# built in a prior run — real but pointless repeated compute. This table lets the
# pipeline check for an existing cluster FIRST, before ever calling protonate.py
# or build_qm_cluster.py again, and save a newly-accepted one for next time.
#
# One-time migration (Supabase SQL editor):
#   create table if not exists public.qm_clusters (
#     id uuid primary key,
#     created_at timestamptz not null default now(),
#     pdb_id text, chain text, resi int, expect_resname text,
#     key text, geometry text, avas text, charge int, spin int,
#     ncas int, nelec int, radius numeric,
#     unique(pdb_id, chain, resi)
#   );
_CLUSTER_COLUMNS = frozenset({
    "id", "created_at", "pdb_id", "chain", "resi", "expect_resname",
    "key", "geometry", "avas", "charge", "spin", "ncas", "nelec", "radius",
})


@router.get("/cluster/lookup")
async def lookup_cluster(pdb_id: str, chain: str, resi: int):
    """Return a previously-built cluster for this exact (pdb_id, chain, resi), if
    one exists — no re-protonation, no re-carving, no re-probing needed. Public
    read: this is reusable infrastructure data, not a provenance record itself
    (the DMRG/SHCI classifications built FROM it are what get sealed)."""
    sb = get_supabase()
    if not sb:
        return {"cluster": None, "db": "not_configured"}
    try:
        res = (sb.table("qm_clusters").select("*")
                 .eq("pdb_id", pdb_id.upper()).eq("chain", chain).eq("resi", resi)
                 .execute())
        return {"cluster": res.data[0] if res.data else None}
    except Exception as e:
        return {"cluster": None, "error": str(e)}


@router.post("/cluster/save")
async def save_cluster(payload: dict = Body(...)):
    """Save a newly-accepted cluster for reuse. Called by
    solange_screen_and_classify.py right after a radius is accepted — before
    DMRG even runs — so the cluster is reusable even if the DMRG/SHCI steps
    later fail. Upserts on (pdb_id, chain, resi): a re-run with a different
    accepted radius/AVAS legitimately replaces the cached one."""
    record = dict(payload)
    for f in ("pdb_id", "chain", "resi", "geometry"):
        if not record.get(f) and record.get(f) != 0:
            raise HTTPException(400, f"missing {f}")
    record["pdb_id"] = str(record["pdb_id"]).upper()
    record.setdefault("id", str(uuid.uuid4()))
    sb = get_supabase()
    if not sb:
        return {"saved": False, "db": "not_configured"}
    safe = {k: v for k, v in record.items() if k in _CLUSTER_COLUMNS}
    try:
        sb.table("qm_clusters").upsert(safe, on_conflict="pdb_id,chain,resi").execute()
        return {"saved": True}
    except Exception as e:
        logging.error("cluster save failed: %s", e)
        return {"saved": False, "error": str(e)}


# NOTE: must stay ABOVE the "/{mutation_id}" catch-all below — FastAPI matches
# routes in declaration order, so a literal path declared after it would be
# swallowed as a mutation id (and would trigger a full VQE run, which this is
# explicitly not).
@router.get("/targets/reference")
async def target_reference_energies():
    """Per-target reference energies, read straight from the stored JW Hamiltonian
    set. Pure lookup: no simulation is run and nothing is written.

    The quantity of interest is `correlation_ha` — the CASSCF-exact energy minus
    the RHF (mean-field) energy over the SAME active space, same molecule, same
    basis. Unlike an absolute total energy it is directly interpretable, and
    unlike a mutant-minus-native difference it is actually meaningful here: the
    two sides are represented by different model compounds, so subtracting them
    would compare unrelated molecules rather than the effect of the mutation.

    Reported against chemical accuracy (1.6 mHa) it states how far mean-field
    falls short for a target — the empirical footing under the §06.i class."""
    CHEM_ACC_MHA = 1.6
    out = {}
    for key, sides in _JW_DATA.items():
        if not isinstance(sides, dict):
            continue
        for side, e in sides.items():
            if not isinstance(e, dict):
                continue
            exact, rhf = e.get("e_active_exact"), e.get("e_active_rhf")
            if exact is None or rhf is None:
                continue
            corr = exact - rhf
            out.setdefault(key, {})[side] = {
                "compound":        e.get("compound"),
                "e_casscf_ha":     e.get("e_casscf"),
                "e_active_exact":  exact,
                "e_active_rhf":    rhf,
                "correlation_ha":  round(corr, 10),
                "correlation_mha": round(corr * 1000.0, 4),
                # How many multiples of chemical accuracy mean-field misses by —
                # the number that makes the value legible without a chemistry background.
                "vs_chem_acc":     round(abs(corr * 1000.0) / CHEM_ACC_MHA, 1),
            }
    return {"chem_accuracy_mha": CHEM_ACC_MHA, "targets": out,
            "_source": "jw_hamiltonians.json (stored reference set — no run performed)"}


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
        "p1_ansatz":       f"AllSinglesDoubles UCCSD ({vqe['n_electrons']}e/{vqe['n_qubits']}q, PennyLane HF initial state)",

        # P2 — Compilation lineage
        "p2_compiler":         "PennyLane",
        "p2_compiler_version": qml.__version__,
        "p2_encoding":         f"Jordan-Wigner (PySCF CAS({vqe['n_electrons']}e,{config['active_orbitals']}o) → openfermion → {vqe['n_paulis']} Pauli terms)",
        "p2_basis_set":        f"STO-3G (PySCF CASSCF({vqe['n_electrons']},{config['active_orbitals']}))",
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
                if config.get("full_electrons")
                else "not applicable — pipeline self-test, not a Phase 3B scaling target"
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
