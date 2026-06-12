"""
generate_expansion_jw.py — Real PySCF CASSCF(2,2) + Jordan-Wigner transforms
for SOLANGE expansion gene LOF targets.

Each expansion gene is assigned the model compound that represents its key
functional residue — the residue whose electronic structure is most disrupted
by the LOF mutation. Native = functional residue. Mutant = disrupted state.

Two new model compounds are added to the existing set:
  - imidazole     (His sidechain) → IDH1 R132H, IDH2 R172H, AXIN2 R815H
  - propionic acid (Glu sidechain) → POLE, ATM catalytic Glu

All geometries from Engh & Huber (1991) idealized amino acid geometry.
PySCF: RHF → CASSCF(2,2)/STO-3G. JW: openfermion InteractionOperator.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# ── New model compound geometries ─────────────────────────────────────────────

# Imidazole (C3H4N2) — His sidechain proxy
# Aromatic N-heterocycle, planar. C–N: 1.382 Å (N3–C4, N1–C2), C–C: 1.370 Å.
# Key for: IDH1/IDH2 neomorphic R132H/R172H, AXIN2 R815H
IMIDAZOLE = [
    ('N', ( 0.000,  0.000,  0.000)),  # N1  (pyridine-like, no H)
    ('C', ( 1.316,  0.000,  0.000)),  # C2  (between two N)
    ('N', ( 1.897,  1.163,  0.000)),  # N3  (pyrrole-like, has H)
    ('C', ( 1.030,  2.122,  0.000)),  # C4
    ('C', (-0.279,  1.832,  0.000)),  # C5
    ('H', ( 1.766, -0.952,  0.000)),  # H2  (on C2)
    ('H', ( 2.965,  1.272,  0.000)),  # H3  (N3–H, pyrrole N)
    ('H', ( 1.367,  3.152,  0.000)),  # H4  (on C4)
    ('H', (-1.238,  2.330,  0.000)),  # H5  (on C5)
]

# Propionic acid (C3H6O2, CH3CH2COOH) — Glu sidechain proxy
# One methylene longer than acetic acid. C–C: 1.520 Å, C=O: 1.214 Å, C–O: 1.364 Å.
# Key for: POLE exonuclease Glu272, ATM Glu catalytic motif
PROPIONIC_ACID = [
    ('C', ( 0.000,  0.000,  0.000)),  # methyl C
    ('C', ( 1.520,  0.000,  0.000)),  # methylene C
    ('C', ( 2.040,  1.424,  0.000)),  # carbonyl C
    ('O', ( 3.256,  1.610,  0.000)),  # C=O
    ('O', ( 1.240,  2.420,  0.000)),  # C–OH
    ('H', ( 0.690,  3.290,  0.000)),  # O–H
    ('H', (-0.390,  1.026,  0.000)),  # methyl H1
    ('H', (-0.390, -0.513,  0.889)),  # methyl H2
    ('H', (-0.390, -0.513, -0.889)),  # methyl H3
    ('H', ( 1.910, -0.513,  0.889)),  # methylene H1
    ('H', ( 1.910, -0.513, -0.889)),  # methylene H2
]

# Methanol (CH3OH) — Ser/Thr sidechain proxy
# C–O: 1.430 Å. O–H: 0.963 Å. C–O–H: 108.5°. C–H: 1.090 Å.
# Key for: VHL Ser65/Thr, FGFR3 Ser kinase-domain, POLE Ser/Thr substrate contact
METHANOL = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('O', ( 1.430,  0.000,  0.000)),
    ('H', ( 1.743,  0.927,  0.000)),  # O–H
    ('H', (-0.390,  1.026,  0.000)),  # C–H1
    ('H', (-0.390, -0.513,  0.889)),  # C–H2
    ('H', (-0.390, -0.513, -0.889)),  # C–H3
]

# Existing geometries (from pyscf_mutations.py — same idealized Engh & Huber)
METHANETHIOL = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('S', ( 0.000,  0.000,  1.819)),
    ('H', ( 1.330,  0.000,  1.971)),
    ('H', ( 1.026,  0.000, -0.363)),
    ('H', (-0.513,  0.889, -0.363)),
    ('H', (-0.513, -0.889, -0.363)),
]
TOLUENE = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('C', ( 0.000,  0.000,  1.512)),
    ('C', ( 1.209,  0.000,  2.210)),
    ('C', ( 1.209,  0.000,  3.602)),
    ('C', ( 0.000,  0.000,  4.300)),
    ('C', (-1.209,  0.000,  3.602)),
    ('C', (-1.209,  0.000,  2.210)),
    ('H', ( 1.026,  0.000, -0.363)),
    ('H', (-0.513,  0.889, -0.363)),
    ('H', (-0.513, -0.889, -0.363)),
    ('H', ( 2.156,  0.000,  1.672)),
    ('H', ( 2.156,  0.000,  4.140)),
    ('H', ( 0.000,  0.000,  5.383)),
    ('H', (-2.156,  0.000,  4.140)),
    ('H', (-2.156,  0.000,  1.672)),
]
P_CRESOL = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('C', ( 0.000,  0.000,  1.512)),
    ('C', ( 1.209,  0.000,  2.210)),
    ('C', ( 1.209,  0.000,  3.602)),
    ('C', ( 0.000,  0.000,  4.300)),
    ('C', (-1.209,  0.000,  3.602)),
    ('C', (-1.209,  0.000,  2.210)),
    ('H', ( 1.026,  0.000, -0.363)),
    ('H', (-0.513,  0.889, -0.363)),
    ('H', (-0.513, -0.889, -0.363)),
    ('H', ( 2.156,  0.000,  1.672)),
    ('H', ( 2.156,  0.000,  4.140)),
    ('H', (-2.156,  0.000,  4.140)),
    ('H', (-2.156,  0.000,  1.672)),
    ('O', ( 0.000,  0.000,  5.664)),
    ('H', ( 0.908,  0.000,  5.985)),
]
FORMAMIDE = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('O', ( 0.000,  1.193,  0.000)),
    ('N', ( 1.121, -0.756,  0.000)),
    ('H', (-0.925, -0.600,  0.000)),
    ('H', ( 2.014, -0.300,  0.000)),
    ('H', ( 1.034, -1.754,  0.000)),
]
ACETIC_ACID = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('C', ( 1.520,  0.000,  0.000)),
    ('O', ( 2.127,  1.051,  0.000)),
    ('O', ( 2.139, -1.215,  0.000)),
    ('H', ( 1.440, -1.890,  0.000)),
    ('H', (-0.390,  1.026,  0.000)),
    ('H', (-0.390, -0.513,  0.889)),
    ('H', (-0.390, -0.513, -0.889)),
]
ACETAMIDE = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('C', ( 1.522,  0.000,  0.000)),
    ('O', ( 2.151,  1.046,  0.000)),
    ('N', ( 2.105, -1.250,  0.000)),
    ('H', ( 3.100, -1.355,  0.000)),
    ('H', ( 1.546, -2.079,  0.000)),
    ('H', (-0.390,  1.026,  0.000)),
    ('H', (-0.390, -0.513,  0.889)),
    ('H', (-0.390, -0.513, -0.889)),
]
GUANIDINE = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('N', ( 1.352,  0.000,  0.000)),
    ('N', (-0.676,  1.171,  0.000)),
    ('N', (-0.676, -1.171,  0.000)),
    ('H', ( 1.857,  0.875,  0.000)),
    ('H', (-0.200,  2.050,  0.000)),
    ('H', (-1.660,  1.110,  0.000)),
    ('H', (-0.200, -2.050,  0.000)),
    ('H', (-1.660, -1.110,  0.000)),
]
ISOBUTANE = [
    ('C', ( 0.000,  0.000,  0.000)),
    ('C', ( 0.000,  0.000,  1.522)),
    ('C', ( 1.436,  0.000, -0.507)),
    ('C', (-0.718,  1.243, -0.507)),
    ('H', ( 0.000,  0.000, -1.090)),
    ('H', ( 1.028,  0.000,  1.885)),
    ('H', (-0.514,  0.890,  1.885)),
    ('H', (-0.514, -0.890,  1.885)),
    ('H', ( 2.121,  0.000,  0.342)),
    ('H', ( 1.607,  0.891, -1.113)),
    ('H', ( 1.607, -0.891, -1.113)),
    ('H', (-1.050,  1.819,  0.285)),
    ('H', ( 0.308,  0.916, -0.870)),
    ('H', (-1.745,  1.526, -0.870)),
]

GEOM = {
    'methanethiol':  METHANETHIOL,
    'toluene':       TOLUENE,
    'p_cresol':      P_CRESOL,
    'formamide':     FORMAMIDE,
    'acetic_acid':   ACETIC_ACID,
    'acetamide':     ACETAMIDE,
    'guanidine':     GUANIDINE,
    'isobutane':     ISOBUTANE,
    'imidazole':     IMIDAZOLE,
    'propionic_acid': PROPIONIC_ACID,
    'methanol':      METHANOL,
}

# ── Expansion gene LOF targets → native/mutant compound assignments ─────────
# Basis: the key functional residue whose electronic structure is disrupted
# by the LOF mutation. References from UniProt, PDB, and cancer mutation databases.
#
# native  = model compound for the functional residue in the wild-type gene
# mutant  = model compound for the residue in the LOF state
#           (specific substitution if known, otherwise loss of functional group)
# note    = biochemical justification
EXPANSION_MODELS = {
    # ── Neuroendocrine / glioma ──────────────────────────────────────────────

    "IDH1_LOF": {
        "gene": "IDH1", "full_electrons": 22, "full_qubits": 44,
        # Neomorphic R132H: Arg132 → His132 at isocitrate-binding site
        # Native: guanidinium (Arg132, electrostatic contact with isocitrate C1-carboxylate)
        # Mutant: imidazole (His132, neomorphic pocket produces 2-HG)
        # Ref: Ward PS et al. (2010) Nature 462:739; Dang L et al. (2009) Nature 462:739
        "native": {"compound": "guanidine",  "residue": "Arg132", "note": "Arg guanidinium contacts isocitrate C1-carboxylate"},
        "mutant": {"compound": "imidazole",  "residue": "His132",  "note": "His imidazole creates neomorphic 2-HG-producing pocket"},
        "pdb": "1T0L",
    },
    "IDH2_LOF": {
        "gene": "IDH2", "full_electrons": 22, "full_qubits": 44,
        # Neomorphic R172H/K: Arg172 → His172 at mitochondrial isocitrate-binding site
        # Same mechanism as IDH1 R132H (active site Arg → His)
        # Ref: Yan H et al. (2009) NEJM 360:765; Losman JA et al. (2013) Science 339:1621
        "native": {"compound": "guanidine",  "residue": "Arg172", "note": "Arg172 contacts isocitrate; same mechanism as IDH1 R132"},
        "mutant": {"compound": "imidazole",  "residue": "His172",  "note": "His172 neomorphic pocket; same JW as IDH1 R132H mutant"},
        "pdb": "1LWD",
    },
    "ATRX_LOF": {
        "gene": "ATRX", "full_electrons": 30, "full_qubits": 60,
        # SNF2 helicase domain ATPase Walker B Asp (Asp2104 in ATRX): catalytic Asp
        # anchors Mg²⁺-ATP. LOF mutations frequently disrupt the Walker B motif.
        # Native: acetic acid (Asp2104 carboxylate, DFG-like Walker B motif)
        # Mutant: acetamide (Asn substitution — common Walker B LOF pattern)
        # Ref: Clynes D et al. (2013) Nat Struct Mol Biol 20:814
        "native": {"compound": "acetic_acid", "residue": "Asp2104", "note": "Walker B Asp2104 carboxylate coordinates Mg²⁺-ATP"},
        "mutant": {"compound": "acetamide",   "residue": "Asn2104", "note": "Asn substitution abolishes Mg²⁺ coordination; ATPase LOF"},
        "pdb": "AlphaFold-Q92793",
    },

    # ── Pan-cancer ────────────────────────────────────────────────────────────

    "POLE_LOF": {
        "gene": "POLE", "full_electrons": 24, "full_qubits": 48,
        # Exonuclease domain Glu272 (ExoIII motif Asp-Glu-Asp triad).
        # POLE exonuclease LOF ultramutators typically hit Asp-Glu contacts.
        # Native: propionic acid (Glu272 γ-carboxylate, one CH2 longer than Asp)
        # Mutant: acetic_acid (shorter carboxylate; represents loss of Glu reach to metal ion)
        # Ref: Palles C et al. (2013) Nat Genet 45:136; Shinbrot E et al. (2014) Genome Res 24:1730
        "native": {"compound": "propionic_acid", "residue": "Glu272", "note": "ExoIII Glu272 γ-carboxylate coordinates 3'-exonuclease Mg²⁺"},
        "mutant": {"compound": "acetic_acid",    "residue": "Asp272", "note": "Asp substitution (shorter reach) impairs metal coordination; proofreading LOF"},
        "pdb": "4M8O",
    },
    "SMARCA4_LOF": {
        "gene": "SMARCA4", "full_electrons": 40, "full_qubits": 80,
        # SWI/SNF ATPase domain Walker B Glu479 (Glu in DExx box, equivalent to DFG Glu).
        # SNF2 helicase Glu479 (human BRG1/SMARCA4) is the Walker B catalytic Glu.
        # Native: propionic acid (Glu479 carboxylate)
        # Mutant: acetic_acid (shorter carboxylate — loss of metal coordination reach)
        # Ref: Kadoch C & Crabtree GR (2015) Sci Adv 1:e1500447
        "native": {"compound": "propionic_acid", "residue": "Glu479", "note": "Walker B Glu479 in DExx box coordinates Mg²⁺-ATP"},
        "mutant": {"compound": "acetic_acid",    "residue": "Asp479", "note": "Asp substitution reduces metal-coordination reach; ATPase LOF"},
        "pdb": "6LTJ",
    },
    "ARID1A_LOF": {
        "gene": "ARID1A", "full_electrons": 28, "full_qubits": 56,
        # ARID domain Trp1815 (ARID domain of BAF250a) is a conserved Trp
        # critical for DNA-minor-groove recognition and chromatin remodeling.
        # Native: indole? We don't have indole — use toluene (closest available aromatic proxy)
        # Actually: ARID domain uses Arg-Thr-Pro wedge with Trp contacts. Key disrupted residue
        # in cancer is often a splice acceptor loss → truncation. For point mutations,
        # Trp→stop or Glu→Gly (truncation of ARID domain).
        # Use toluene (Phe-like aromatic) as best available proxy for aromatic stacking LOF.
        # Ref: Guan B et al. (2011) PNAS 108:17aerr; Mathur R et al. (2017) Nat Genet 49:407
        "native": {"compound": "toluene",    "residue": "Trp1815", "note": "ARID-domain Trp1815 aromatic stacking in DNA minor groove (toluene proxy)"},
        "mutant": {"compound": "isobutane",  "residue": "LOF",     "note": "Loss of aromatic stacking → aliphatic proxy (disrupted ARID-DNA contact)"},
        "pdb": "2L9X",
    },
    "BRCA1_LOF": {
        "gene": "BRCA1", "full_electrons": 32, "full_qubits": 64,
        # RING domain Cys24/Cys27/His41/Cys44 — canonical C3HC4 zinc finger
        # Native Cys: methanethiol (Cys24/Cys44 Zn-coordinating sulfhydryl)
        # Mutant LOF: formamide (backbone-only — loss of Cys sidechain → Gly equivalent)
        # Ref: Meza JE et al. (1999) J Biol Chem 274:5659; Morris JR & Solomon E (2004) Hum Mol Genet 13:807
        "native": {"compound": "methanethiol", "residue": "Cys44",  "note": "RING-domain Cys44 thiolate coordinates BARD1-binding Zn²⁺"},
        "mutant": {"compound": "formamide",    "residue": "LOF",    "note": "Loss of Cys sidechain → backbone amide only; Zn-finger disrupted"},
        "pdb": "1JM7",
    },
    "BRCA2_LOF": {
        "gene": "BRCA2", "full_electrons": 32, "full_qubits": 64,
        # BRC repeat Phe3175 is the canonical Rad51-binding Phe within the DBD.
        # C-terminal OB fold Phe3175 π-stacks with Rad51 Phe82 at the interface.
        # LOF frequently disrupts this aromatic contact.
        # Native: toluene (Phe3175 aromatic)
        # Mutant: isobutane (aliphatic substitution — loss of π-stacking)
        # Ref: Lo T et al. (2003) J Biol Chem 278:14; Pellegrini L et al. (2002) Nature 420:287
        "native": {"compound": "toluene",   "residue": "Phe3175", "note": "OB-fold Phe3175 π-stacks with Rad51 Phe82 at HR interface"},
        "mutant": {"compound": "isobutane", "residue": "LOF",     "note": "Aliphatic substitution abolishes Rad51 π-stacking; HR LOF"},
        "pdb": "1MJE",
    },
    "ATM_LOF": {
        "gene": "ATM", "full_electrons": 28, "full_qubits": 56,
        # PIKK-family PI3K-like kinase domain Asp2870 (DFG motif equivalent).
        # ATM catalytic Asp2870 is required for autophosphorylation and substrate phosphorylation.
        # Native: acetic_acid (Asp2870 carboxylate)
        # Mutant: acetamide (Asn substitution — kinase-dead equivalent)
        # Ref: Lavin MF & Shiloh Y (1997) Annu Rev Biochem 66:169; Barlow C et al. (1996) Cell 86:159
        "native": {"compound": "acetic_acid", "residue": "Asp2870", "note": "PI3K-like DFG Asp2870 catalytic carboxylate"},
        "mutant": {"compound": "acetamide",   "residue": "Asn2870", "note": "Asn substitution → kinase-dead; loss of DNA damage response"},
        "pdb": "AlphaFold-Q13315",
    },
    "TERT_LOF": {
        "gene": "TERT", "full_electrons": 24, "full_qubits": 48,
        # RT active site Asp712 (YMDD → YMDA kinase-dead motif in telomerase RT).
        # Asp712 coordinates the two Mg²⁺ ions for nucleotide transfer.
        # Native: acetic_acid (Asp712 carboxylate)
        # Mutant: acetamide (Asn substitution — telomerase-dead)
        # Ref: Lingner J et al. (1997) Science 276:561; Harrington L et al. (1997) Science 275:973
        "native": {"compound": "acetic_acid", "residue": "Asp712", "note": "RT YMDD motif Asp712 coordinates dual-Mg²⁺ nucleotide transfer"},
        "mutant": {"compound": "acetamide",   "residue": "Asn712", "note": "Asn substitution abolishes Mg²⁺ coordination; telomerase-dead"},
        "pdb": "7LYT",
    },
    "RB1_LOF": {
        "gene": "RB1", "full_electrons": 32, "full_qubits": 64,
        # RB1 large pocket Glu2 and acidic patch: key acidic residues in the A/B domains.
        # The E2F-binding interface is dominated by acidic residues (Asp/Glu).
        # Native: propionic_acid (Glu in E2F-binding acidic patch)
        # Mutant: methanol (Ser substitution — partial loss of E2F contact)
        # Ref: Weinberg RA (1995) Cell 81:323; Rubin SM et al. (2005) PNAS 102:7660
        "native": {"compound": "propionic_acid", "residue": "Glu2",  "note": "E2F-binding acidic patch Glu; key electrostatic E2F contact"},
        "mutant": {"compound": "methanol",       "residue": "Ser",   "note": "Ser substitution reduces E2F-binding affinity; cell cycle checkpoint LOF"},
        "pdb": "2AZE",
    },
    "NF1_LOF": {
        "gene": "NF1", "full_electrons": 28, "full_qubits": 56,
        # RasGAP catalytic Arg1276 (arginine finger) — inserts into Ras active site
        # and stabilizes the transition state for GTP hydrolysis.
        # Native: guanidine (Arg1276 guanidinium finger)
        # Mutant: formamide (backbone/Gly substitution — loss of Arg sidechain)
        # Ref: Scheffzek K et al. (1998) Science 277:333; Cichowski K & Jacks T (2001) Cell 104:593
        "native": {"compound": "guanidine", "residue": "Arg1276", "note": "Catalytic Arg-finger 1276 guanidinium stabilises GTPase transition state"},
        "mutant": {"compound": "formamide", "residue": "Gly1276", "note": "Loss of Arg sidechain → Gly/backbone only; GAP activity abolished"},
        "pdb": "1NF1",
    },
    "NF2_LOF": {
        "gene": "NF2", "full_electrons": 22, "full_qubits": 44,
        # FERM domain: canonical F1-F2-F3 subdomain structure. Key Arg in β1 of F3.
        # Merlin (NF2) FERM Arg341 contacts the plasma membrane phospholipid headgroups.
        # Native: guanidine (Arg341 guanidinium electrostatic contact with membrane)
        # Mutant: acetamide (Gln substitution — common cancer LOF in FERM domain)
        # Ref: Bretscher A et al. (2002) Nat Rev Mol Cell Biol 3:586
        "native": {"compound": "guanidine", "residue": "Arg341", "note": "FERM F3-domain Arg341 guanidinium contacts PI(4,5)P2 headgroup"},
        "mutant": {"compound": "acetamide", "residue": "Gln341", "note": "Gln substitution reduces membrane affinity; scaffold LOF"},
        "pdb": "1H4R",
    },
    "AXIN1_LOF": {
        "gene": "AXIN1", "full_electrons": 24, "full_qubits": 48,
        # DIX domain: Phe631 (human AXIN1) is required for head-to-tail polymerization
        # of AXIN1 dimers within the β-catenin destruction complex.
        # Native: toluene (Phe631 aromatic sidechain in DIX hydrophobic interface)
        # Mutant: isobutane (aliphatic substitution — loss of aromatic packing)
        # Ref: Schwarz-Romond T et al. (2007) Nat Struct Mol Biol 14:484
        "native": {"compound": "toluene",   "residue": "Phe631", "note": "DIX-domain Phe631 aromatic core of head-to-tail polymerization interface"},
        "mutant": {"compound": "isobutane", "residue": "LOF",    "note": "Aliphatic substitution disrupts DIX polymerization; β-catenin destruction complex LOF"},
        "pdb": "1WSP",
    },
    "AXIN2_LOF": {
        "gene": "AXIN2", "full_electrons": 24, "full_qubits": 48,
        # c.2444G>A → p.Arg815His (AXIN2 codon 815, position 2 of codon CGT/CGC → CAT/CAC).
        # Arg815 is in the AXIN2 C-terminal regulatory domain; His substitution disrupts
        # electrostatic contacts with APC and destabilizes the destruction complex.
        # Native: guanidine (Arg815 electrostatic contact with APC SAMP repeats)
        # Mutant: imidazole (His815 — reduced positive charge; weaker APC binding)
        # Ref: Salahshor S & Woodgett JR (2005) Neoplasia 7:867
        "native": {"compound": "guanidine",  "residue": "Arg815", "note": "C-terminal Arg815 electrostatic contact with APC SAMP repeats"},
        "mutant": {"compound": "imidazole",  "residue": "His815", "note": "His815 (c.2444G>A) — reduced charge; weakened APC contact; WNT LOF"},
        "pdb": "AlphaFold-O15169",
    },
    "CDKN2A_LOF": {
        "gene": "CDKN2A", "full_electrons": 20, "full_qubits": 40,
        # INK4 ankyrin repeat Arg58 (p16-INK4A): key electrostatic contact with CDK4.
        # Arg58 guanidinium forms a salt bridge with Asp97 of CDK4 at the binding interface.
        # Native: guanidine (Arg58 guanidinium CDK4 contact)
        # Mutant: methanol (Ser substitution — one of the most common INK4A point mutations)
        # Ref: Kamb A et al. (1994) Science 264:436; Ruas M & Peters G (1998) Biochim Biophys Acta 1378:F115
        "native": {"compound": "guanidine", "residue": "Arg58", "note": "p16-INK4A Arg58 salt bridge with CDK4 Asp97"},
        "mutant": {"compound": "methanol",  "residue": "Ser58", "note": "Ser58 substitution breaks CDK4 salt bridge; CDK4/6 cycle checkpoint LOF"},
        "pdb": "1BI7",
    },

    # ── Urological / RCC ─────────────────────────────────────────────────────

    "VHL_LOF": {
        "gene": "VHL", "full_electrons": 25, "full_qubits": 50,
        # VHL β-domain Pro95 is the key residue in the HIF-α proline hydroxylation
        # recognition site. The VHL HIF-1α-binding domain uses Tyr112 and Ser111
        # to contact the hydroxyl on HIF-1α Pro402/Pro564.
        # Native: methanol (Ser111 hydroxyl contacts hydroxyproline on HIF-1α)
        # Mutant: formamide (Gly backbone — common cancer missense abolishes Ser contact)
        # Ref: Minervini R et al. (2020) PNAS 117:29001; Ivan M et al. (2001) Science 292:464
        "native": {"compound": "methanol",  "residue": "Ser111", "note": "β-domain Ser111 OH contacts HIF-1α hydroxyproline Pro564"},
        "mutant": {"compound": "formamide", "residue": "Gly111", "note": "Loss of Ser sidechain → Gly; disrupts HIF-1α recognition; HIF pathway deregulated"},
        "pdb": "1LM8",
    },
    "BAP1_LOF": {
        "gene": "BAP1", "full_electrons": 35, "full_qubits": 70,
        # UCH (ubiquitin carboxyl-terminal hydrolase) domain: Cys91 is the catalytic Cys
        # in the Cys-His-Asp/Asn catalytic triad.
        # Native: methanethiol (Cys91 catalytic thiolate)
        # Mutant: formamide (backbone — loss of Cys sidechain → Gly equivalent)
        # Ref: Ventii KH et al. (2008) Cancer Res 68:6953; Murali R et al. (2013) Lancet Oncol 14:e254
        "native": {"compound": "methanethiol", "residue": "Cys91",  "note": "UCH-domain Cys91 catalytic thiolate; deubiquitinase active site"},
        "mutant": {"compound": "formamide",    "residue": "Gly91",  "note": "Loss of Cys91 sidechain → backbone only; deubiquitinase abolished"},
        "pdb": "3KVF",
    },
    "PBRM1_LOF": {
        "gene": "PBRM1", "full_electrons": 28, "full_qubits": 56,
        # BD2 bromodomain: Tyr1242 (PBRM1 BD2) is the conserved Tyr in the WPF shelf
        # that contacts acetylated Lys on histones via π–cation stacking.
        # Native: p_cresol (Tyr1242 phenolic hydroxyl + aromatic ring for Kac recognition)
        # Mutant: toluene (loss of OH → Phe equivalent; weakened Kac binding)
        # Ref: Xu Y et al. (2012) Structure 20:513; Kadoch C et al. (2013) Nat Genet 45:592
        "native": {"compound": "p_cresol", "residue": "Tyr1242", "note": "BD2 WPF-shelf Tyr1242 phenol contacts acetyl-Lys on H3/H4"},
        "mutant": {"compound": "toluene",  "residue": "Phe1242", "note": "Loss of OH → Phe; reduced Kac recognition; chromatin targeting LOF"},
        "pdb": "3G0L",
    },
    "SETD2_LOF": {
        "gene": "SETD2", "full_electrons": 30, "full_qubits": 60,
        # SET domain: Tyr1666 (human SETD2) is the conserved SET-domain Tyr (i-SET helix)
        # that positions H3K36 for methyl transfer via a π-cation interaction with the
        # substrate Lys ε-amino group.
        # Native: p_cresol (Tyr1666 phenol for π-cation / H3K36 substrate positioning)
        # Mutant: toluene (Phe substitution — loss of OH disrupts substrate positioning)
        # Ref: Bhatt DL (2022); Wagner EJ & Bhatt DL (2012) Mol Cell 46:736
        "native": {"compound": "p_cresol", "residue": "Tyr1666", "note": "i-SET Tyr1666 phenol π-cation with H3K36 substrate amino group"},
        "mutant": {"compound": "toluene",  "residue": "Phe1666", "note": "Phe substitution loses H3K36 orientation contact; H3K36me3 abolished"},
        "pdb": "5JLB",
    },
    "FGFR3_LOF": {
        "gene": "FGFR3", "full_electrons": 24, "full_qubits": 48,
        # Kinase domain DFG motif Asp641 (human FGFR3): catalytic Asp coordinating Mg²⁺-ATP.
        # In urothelial carcinoma LOF context, inactivating mutations hit the kinase domain.
        # Native: acetic_acid (Asp641 carboxylate in DFG motif)
        # Mutant: methanol (Ser641 substitution — common GOF but LOF possible; or Gly)
        # Ref: Helsten T et al. (2016) Clin Cancer Res 22:259; Chesi M et al. (1997) Nat Genet 16:260
        "native": {"compound": "acetic_acid", "residue": "Asp641", "note": "DFG-motif Asp641 carboxylate coordinates Mg²⁺-ATP γ-phosphate transfer"},
        "mutant": {"compound": "methanol",    "residue": "Ser641", "note": "Ser substitution reduces Mg²⁺ coordination; kinase domain LOF"},
        "pdb": "4K33",
    },
    "TSC1_LOF": {
        "gene": "TSC1", "full_electrons": 26, "full_qubits": 52,
        # Hamartin (TSC1) Arg692 in the coiled-coil domain mediates TSC1–TSC2 heterodimerisation.
        # This electrostatic contact is essential for the TSC1/TSC2 GAP complex stability.
        # Native: guanidine (Arg692 electrostatic contact with TSC2 acidic patch)
        # Mutant: formamide (Gly/backbone equivalent — loss of Arg sidechain)
        # Ref: Chong-Kopera H et al. (2006) J Biol Chem 281:29542
        "native": {"compound": "guanidine", "residue": "Arg692", "note": "Coiled-coil Arg692 salt bridge with TSC2 acidic patch; TSC complex stability"},
        "mutant": {"compound": "formamide", "residue": "Gly692", "note": "Loss of Arg sidechain destabilises TSC1-TSC2 complex; mTORC1 deregulated"},
        "pdb": "AlphaFold-Q92574",
    },
    "TSC2_LOF": {
        "gene": "TSC2", "full_electrons": 30, "full_qubits": 60,
        # GAP domain Arg1743 (equivalent to NF1 catalytic Arg-finger) inserts into Rheb
        # active site to stabilise the transition state for GTP hydrolysis.
        # Native: guanidine (Arg1743 Arg-finger guanidinium)
        # Mutant: formamide (Gly substitution — loss of Arg finger = GAP-dead)
        # Ref: Inoki K et al. (2003) Cell 115:577; Garami A et al. (2003) Mol Cell 11:1457
        "native": {"compound": "guanidine", "residue": "Arg1743", "note": "GAP-domain Arg1743 catalytic Arg-finger contacts Rheb GTP γ-phosphate"},
        "mutant": {"compound": "formamide", "residue": "Gly1743", "note": "Arg-finger deletion → Rheb-GAP activity abolished; mTORC1 constitutively active"},
        "pdb": "5EJO",
    },
}


# ── PySCF CASSCF(2,2) + JW transform ─────────────────────────────────────────

def geom_to_string(atom_list):
    return '\n'.join(f'{sym}  {x:.4f}  {y:.4f}  {z:.4f}' for sym, (x, y, z) in atom_list)


def run_casscf(compound_name, verbose=0):
    """Run RHF → CASSCF(2,2)/STO-3G for the given compound."""
    from pyscf import gto, scf, mcscf, ao2mo

    atom_list = GEOM[compound_name]
    mol = gto.Mole()
    mol.atom   = geom_to_string(atom_list)
    mol.basis  = 'sto-3g'
    mol.charge = 0
    mol.spin   = 0
    mol.verbose = verbose
    mol.build()

    mf = scf.RHF(mol)
    mf.max_cycle = 300
    mf.conv_tol  = 1e-10
    e_rhf = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"RHF did not converge for {compound_name}")

    mc = mcscf.CASSCF(mf, ncas=2, nelecas=2)
    mc.conv_tol       = 1e-9
    mc.conv_tol_grad  = 1e-5
    mc.max_cycle_macro = 150
    e_casscf = mc.kernel()[0]
    if not mc.converged:
        print(f"  WARNING: CASSCF did not fully converge for {compound_name}", file=sys.stderr)

    h1e, ecore = mc.get_h1eff()
    h2e_compressed = mc.get_h2eff()
    h2e = ao2mo.restore(1, h2e_compressed, mc.ncas)

    return {
        'e_rhf':    e_rhf,
        'e_casscf': e_casscf,
        'ecore':    float(ecore),
        'h1e':      h1e,
        'h2e':      h2e,
        'converged': bool(mc.converged),
    }


def build_jw_terms(h1e, h2e, ecore):
    """
    Build Jordan-Wigner Pauli terms from CAS(2e,2o) active-space integrals.

    Uses openfermion InteractionOperator → jordan_wigner transform.
    Returns (terms, e_active_exact, e_active_rhf) where:
      terms           = list of {pauli, coeff} for PennyLane
      e_active_exact  = exact diagonalisation energy of active-space Hamiltonian
      e_active_rhf    = RHF energy contribution from active space
    """
    from openfermion import InteractionOperator, jordan_wigner
    from openfermion.linalg import get_sparse_operator
    import scipy.sparse.linalg as spla

    # Construct InteractionOperator from 1e/2e integrals (2 spatial orbitals, 4 spin-orbitals)
    # Convention: h2e_of[i,j,k,l] = <ij|kl> = h2e[i,k,j,l] (openfermion uses chemist notation)
    n_orb = 2
    one_body = np.zeros((2 * n_orb, 2 * n_orb))
    two_body = np.zeros((2 * n_orb, 2 * n_orb, 2 * n_orb, 2 * n_orb))

    # Spin-orbital mapping: 0=α0, 1=α1, 2=β0, 3=β1
    for p in range(n_orb):
        for q in range(n_orb):
            # Alpha spin
            one_body[2*p,   2*q]   = h1e[p, q]
            # Beta spin
            one_body[2*p+1, 2*q+1] = h1e[p, q]

    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    val = 0.5 * h2e[p, q, r, s]
                    # αα
                    two_body[2*p,   2*q,   2*r,   2*s]   = val
                    # ββ
                    two_body[2*p+1, 2*q+1, 2*r+1, 2*s+1] = val
                    # αβ
                    two_body[2*p,   2*q+1, 2*r,   2*s+1] = val
                    two_body[2*p+1, 2*q,   2*r+1, 2*s]   = val

    ham_op = InteractionOperator(constant=0.0, one_body_tensor=one_body, two_body_tensor=two_body)
    jw_op  = jordan_wigner(ham_op)

    # Exact diagonalisation of 4-qubit Hamiltonian
    sparse = get_sparse_operator(jw_op)
    e_active_exact = float(spla.eigsh(sparse, k=1, which='SA')[0][0])

    # RHF active-space energy: h1e[0,0] + h1e[1,1] + 0.5*(J_aa + 2*J_ab - K_ab)
    # For 2 electrons in 2 orbitals (doubly occupied α=0, β=0 in RHF):
    e_active_rhf = (h1e[0, 0] + h1e[1, 1]
                    + 0.5 * (h2e[0,0,0,0] + 2*h2e[0,0,1,1] - h2e[0,1,1,0]))

    # Convert JW terms to {pauli, coeff} format for PennyLane
    terms = []
    for term, coeff in jw_op.terms.items():
        if abs(coeff) < 1e-12:
            continue
        if len(term) == 0:
            pauli_str = 'I'
        else:
            pauli_str = ' '.join(f'{op}{qubit}' for qubit, op in sorted(term))
        terms.append({'coeff': float(coeff.real), 'pauli': pauli_str})

    return terms, e_active_exact, float(e_active_rhf)


def compute_jw_entry(compound_name, residue_note, cached_casscf=None, verbose=0):
    """Run PySCF + JW for one compound. Returns the full jw_hamiltonians.json entry dict."""
    if cached_casscf is not None:
        r = cached_casscf
    else:
        r = run_casscf(compound_name, verbose=verbose)

    terms, e_active_exact, e_active_rhf = build_jw_terms(r['h1e'], r['h2e'], r['ecore'])

    return {
        'compound':       compound_name,
        'residue_note':   residue_note,
        'ecore':          r['ecore'],
        'e_casscf':       r['e_casscf'],
        'e_active_exact': e_active_exact,
        'e_active_rhf':   e_active_rhf,
        'n_paulis':       len(terms),
        'terms':          terms,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    jw_path = Path(__file__).parent / 'backend' / 'jw_hamiltonians.json'
    with open(jw_path) as f:
        jw_data = json.load(f)

    # Cache PySCF results by compound name to avoid recomputing the same molecule
    casscf_cache = {}

    print(f"\n{'='*70}")
    print(f"  SOLANGE — Expansion gene JW Hamiltonian generation")
    print(f"  PySCF CASSCF(2,2)/STO-3G + openfermion Jordan-Wigner transform")
    print(f"  Total targets: {len(EXPANSION_MODELS)}")
    print(f"{'='*70}\n")

    for gene_id, cfg in EXPANSION_MODELS.items():
        if gene_id in jw_data:
            print(f"  [{gene_id}] already in jw_hamiltonians.json — skipping")
            continue

        print(f"\n  ─── {gene_id} ({cfg['gene']}, {cfg['full_electrons']}e / {cfg['full_qubits']}q) ───")
        print(f"  PDB: {cfg['pdb']}")

        jw_entry = {}
        for side in ('native', 'mutant'):
            compound = cfg[side]['compound']
            residue  = cfg[side]['residue']
            note     = cfg[side]['note']
            print(f"  [{side.upper()}] {residue} → model: {compound}")
            print(f"         {note}")

            if compound not in casscf_cache:
                print(f"         Running PySCF CASSCF(2,2)/STO-3G for {compound}...", end=' ', flush=True)
                casscf_cache[compound] = run_casscf(compound, verbose=0)
                print(f"done. e_casscf = {casscf_cache[compound]['e_casscf']:.8f} Ha")
            else:
                print(f"         Using cached PySCF result for {compound}: e_casscf = {casscf_cache[compound]['e_casscf']:.8f} Ha")

            jw_entry[side] = compute_jw_entry(
                compound_name  = compound,
                residue_note   = note,
                cached_casscf  = casscf_cache[compound],
            )
            print(f"         JW transform: {jw_entry[side]['n_paulis']} Pauli terms")
            print(f"         e_active_exact = {jw_entry[side]['e_active_exact']:.8f} Ha")

        jw_data[gene_id] = jw_entry

        # Save after every gene so partial results survive interruption
        with open(jw_path, 'w') as f:
            json.dump(jw_data, f, indent=2)
        print(f"  ✓ {gene_id} saved to jw_hamiltonians.json")

    print(f"\n{'='*70}")
    print(f"  Expansion gene JW Hamiltonians — complete")
    print(f"  Total entries in jw_hamiltonians.json: {len(jw_data)}")
    print(f"{'='*70}\n")


if __name__ == '__main__':
    main()
