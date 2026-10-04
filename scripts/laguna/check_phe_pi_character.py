"""
Check whether the 16 virtual MOs selected for the TP53 C275F mutant's CAS(36,34)
active space carry real Phe275 ring pi/pi* character, or are generic virtuals
unrelated to the introduced mutation.

Run this on Laguna, in the same directory as the mutant files
(tp53_c275f_mutant_cluster.xyz, tp53_c275f_mutant_cas34_FIXED_mo.npy).

Method:
  1. Auto-detect the phenyl ring from geometry alone (6 carbons, mutually
     bonded in a planar hexagon) -- no need to know atom indices in advance.
  2. For each of the 16 virtual MO columns, compute the Mulliken-style
     population on (a) all AOs centered on the ring, (b) just the ring
     carbons' 2p AOs (the pi-system's building block).
  3. High ring/2p weight on the virtuals -> they are measuring the mutation's
     own chemistry. Low/flat weight -> they are generic, and the active-space
     selection needs to be revisited (does NOT probe what we think it probes).

Caveat, stated explicitly: summing all three 2p Cartesian components (not just
the one perpendicular to the ring plane) does not perfectly isolate pi from
in-plane sigma-p character. It is a legitimate first-pass screen; if the
result is ambiguous (e.g. 30-50%), the follow-up is to project onto the
ring's own normal vector (already computed here as the smallest SVD singular
vector) to isolate true out-of-plane pi weight.
"""
import itertools
import numpy as np
from pyscf import gto
from scipy.spatial.distance import pdist, squareform

# --- must match the solange_dmrg.py run exactly ---
XYZ = "tp53_c275f_mutant_cluster.xyz"
BASIS = "6-31g"
CHARGE = 0
SPIN = 0
MO_PATH = "tp53_c275f_mutant_cas34_FIXED_mo.npy"
NCORE = 341
N_OCC = 18
N_VIR = 16

BOHR_PER_ANGSTROM = 1.8897259886
CC_BOND_MIN = 1.30 * BOHR_PER_ANGSTROM
CC_BOND_MAX = 1.50 * BOHR_PER_ANGSTROM
PLANARITY_TOL = 0.15  # s[2]/s[0] from SVD; small = flat ring


def find_phenyl_rings(mol):
    coords = mol.atom_coords()  # Bohr
    elements = [mol.atom_symbol(i) for i in range(mol.natm)]
    carbon_idx = [i for i, e in enumerate(elements) if e == "C"]
    sub = np.array([coords[i] for i in carbon_idx])
    D = squareform(pdist(sub))
    adj = (D > CC_BOND_MIN) & (D < CC_BOND_MAX)

    candidates = []
    n = len(carbon_idx)
    for combo in itertools.combinations(range(n), 6):
        block = adj[np.ix_(combo, combo)]
        if np.all(block.sum(axis=1) == 2):  # each atom bonded to exactly 2 others -> a cycle
            pts = sub[list(combo)]
            centroid = pts.mean(axis=0)
            _, s, vt = np.linalg.svd(pts - centroid)
            planarity = s[2] / s[0]
            if planarity < PLANARITY_TOL:
                ring_atoms = [carbon_idx[c] for c in combo]
                normal = vt[2]
                candidates.append((ring_atoms, planarity, normal))
    candidates.sort(key=lambda x: x[1])
    return candidates


def main():
    mol = gto.M(atom=XYZ, basis=BASIS, charge=CHARGE, spin=SPIN)
    S = mol.intor("int1e_ovlp")
    mo = np.load(MO_PATH)

    rings = find_phenyl_rings(mol)
    print(f"Found {len(rings)} planar 6-carbon ring candidate(s) in {mol.natm} atoms")
    for ring_atoms, planarity, _ in rings[:5]:
        print(f"  atoms {ring_atoms} (0-based)  planarity={planarity:.4f}")
    if not rings:
        raise SystemExit(
            "No phenyl ring auto-detected -- check CC_BOND_MIN/MAX or PLANARITY_TOL, "
            "or confirm the mutant cluster actually contains the Phe275 ring atoms."
        )

    ring_atoms, planarity, normal = rings[0]
    print(f"\nUsing ring atoms (0-based): {ring_atoms}  (planarity={planarity:.4f})")

    ao_labels = mol.ao_labels(fmt=False)  # (atom_idx, elem, nl, m)
    ring_all_ao = [i for i, (a, *_rest) in enumerate(ao_labels) if a in ring_atoms]
    ring_2p_ao = [i for i, (a, elem, nl, m) in enumerate(ao_labels)
                  if a in ring_atoms and nl.strip() == "2p"]
    print(f"{len(ring_2p_ao)} ring-carbon 2p AOs, {len(ring_all_ao)} total ring AOs "
          f"out of {mol.nao} AOs total")

    vir_start = NCORE + N_OCC
    vir_end = vir_start + N_VIR
    print(f"\nVirtual MO columns in mo_coeff: [{vir_start}:{vir_end}]\n")
    print(f"{'col':>5}  {'ring_2p_%':>10}  {'ring_all_%':>11}")
    SC = S @ mo
    totals_2p, totals_all = [], []
    for col in range(vir_start, vir_end):
        c = mo[:, col]
        pop = c * SC[:, col]  # Mulliken population per AO, sums to 1.0 for a normalized MO
        p_2p = pop[ring_2p_ao].sum()
        p_all = pop[ring_all_ao].sum()
        totals_2p.append(p_2p)
        totals_all.append(p_all)
        print(f"{col:5d}  {100*p_2p:9.1f}%  {100*p_all:10.1f}%")

    print(f"\nMean across the 16 virtuals: ring_2p={100*np.mean(totals_2p):.1f}%  "
          f"ring_all={100*np.mean(totals_all):.1f}%")
    print("\nInterpretation: high and consistent ring weight across most of the 16 "
          "virtuals means they genuinely probe the Phe275 ring's electronic structure. "
          "Low or scattered weight means the selection is picking up generic "
          "correlation elsewhere in the 187-atom cluster, not the mutation site -- "
          "in which case the active-space construction needs revisiting before any "
          "R(mutant) number from this run can be trusted.")


if __name__ == "__main__":
    main()
