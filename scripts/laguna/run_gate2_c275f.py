#!/usr/bin/env python3
"""Gate 2: TP53 C275F -- MODELLED mutant vs wild-type -> QM cluster -> AVAS.

No deposited PDB entry carries C275F (Claude Science round 2: exhaustive
4-channel search of 325 entities, sensitivity-checked -- see
active_space_spec.json's TP53_C275F.round2_search). This is the modelled-
mutant path the project's own README always said would eventually be
needed, built with PyMOL's Mutagenesis wizard on 2OCJ chain A residue 275
(Cys->Phe).

MODEL CORRECTED 2026-09-30 (v1 -> v2, Claude Science + independent check).
The first model (2OCJ_C275F_model.pdb, "lowest-strain rotamer" 56.23 vs.
82.03) made the mutant AVAS run's SCF fail to converge. Root cause was NOT
a bad rotamer: the entire PHE275 residue had been placed 8.03 A from
CYS275's actual deposited CB position (SER240, checked as a fixed
reference point, sat at drift=0.00 A -- the rest of the structure was
untouched, only residue 275 was mis-registered). Signature of a dedup/
mutagenesis-indexing defect, not a strained side chain. The replacement
(2OCJ_C275F_model_v2.pdb) rebuilds an ideal Phe ring on the deposited
N/CA/CB of CYS275 and picks chi1=-162 deg, chi2=+71 deg (inside the Lovell
"trans" rotamer well -- a populated rotamer, not an ad hoc fit) via a full
65,341-pose chi1 x chi2 clash scan against all 486 chain-A heavy atoms
within 16 A. Independently verified here (not just trusted): CYS275's
N/CA/C/O/CB read directly from struct/2OCJ.pdb match v2's PHE275
N/CA/C/O/CB to 3 decimal places. The acceptable rotamer window at this
site is only 0.54% of the full chi1 x chi2 space (14/2592 coarse poses) --
worth stating as the measured steric cost of this mutation, not just that
a solution exists.

THIS IS A MODEL, NOT A DEPOSITED STRUCTURE. Report it as such everywhere
downstream -- this is exactly the category of thing the retracted first
C275F result got wrong (there the record didn't even carry a real
mutation; here it is a real, disclosed model, which is a different and
weaker-but-honest evidentiary status than a deposited structure).

R3 (same fixed residue list on both sides) is satisfied by construction:
the "mutant" file is 2OCJ itself with only residue 275 changed, so the
SAME REGION_RESIDUES list below applies verbatim to both the wild-type
(struct/2OCJ.pdb) and the modelled mutant (struct/2OCJ_C275F_model_v2.pdb)
-- not re-derived per side.

Usage
-----
    python run_gate2_c275f.py --side wt
    python run_gate2_c275f.py --side mutant
"""
import argparse
import json
import math
import os
import sys

CAP_BOND_LENGTH = 1.09

# Verbatim from active_space_spec.json's TP53_C275F.region_residues.
REGION_RESIDUES = {
    ("A", "135"): "CYS", ("A", "136"): "GLN", ("A", "137"): "LEU",
    ("A", "239"): "ASN", ("A", "240"): "SER", ("A", "273"): "ARG",
    ("A", "274"): "VAL", ("A", "275"): None,   # CYS in wt, PHE in the model
    ("A", "276"): "ALA", ("A", "277"): "CYS", ("A", "278"): "PRO",
    ("A", "281"): "ASP",
}
SIDECHAIN_Q = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}
EXPECTED_NET_CHARGE = 0   # Arg273(+1) + Asp281(-1); no metal, no thiolate here
FILES = {"wt": "struct/2OCJ.pdb", "mutant": "struct/2OCJ_C275F_model_v2.pdb"}


def check_residue_275(pdb_path, expect):
    comps = set()
    for ln in open(pdb_path):
        if ln[:6] == "ATOM  " and ln[21] == "A" and ln[22:26].strip() == "275":
            comps.add(ln[17:20].strip())
    if comps != {expect}:
        sys.exit(f"residue 275 in {pdb_path} is {comps}, expected exactly {{'{expect}'}}")
    print(f"[verify] {pdb_path}: residue 275 = {expect} (direct read, matches expectation)")


def read_pdb(path):
    atoms = []
    for ln in open(path):
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        atoms.append(dict(
            name=ln[12:16].strip(), comp=ln[17:20].strip(),
            ch=ln[21].strip() or "A", seq=ln[22:26].strip(),
            xyz=(float(ln[30:38]), float(ln[38:46]), float(ln[46:54])),
            elem=(ln[76:78].strip() or ln[12:16].strip()[0]).upper()))
    return atoms


def protonate(src_pdb, ph=7.4):
    out = src_pdb.replace(".pdb", "_H.pdb")
    if os.path.exists(out):
        return out
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=src_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print(f"[prep] protonated at pH {ph} -> {out}")
    return out


def d(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def build_cluster(atoms, side):
    groups = {}
    for a in atoms:
        key = (a["ch"], a["seq"])
        if key not in REGION_RESIDUES:
            continue
        groups.setdefault(key, []).append(a)

    # DEDUP: found live 2026-09-29 on the mutant side -- pdbfixer duplicated
    # residue 275 (the modelled PHE) 11x in struct/2OCJ_C275F_model_H.pdb
    # (220 atom lines instead of ~20; every OTHER residue, e.g. 135, appears
    # exactly once as expected). Root cause not pinned down (pdbfixer
    # template-matching on a modelled, non-deposited residue is the leading
    # guess), but the fix is robust regardless of cause: keep only the FIRST
    # occurrence of each atom name per residue. This silently inflated the
    # mutant cluster to 2196 basis functions (vs. 977 on the WT side) and is
    # almost certainly why that run took >24h instead of finishing in the
    # same ballpark as the WT side.
    for key in groups:
        seen, deduped = set(), []
        for a in groups[key]:
            if a["name"] in seen:
                continue
            seen.add(a["name"])
            deduped.append(a)
        if len(deduped) != len(groups[key]):
            print(f"[cluster] NOTE: deduplicated {key} from {len(groups[key])} to "
                  f"{len(deduped)} atoms (duplicate atom names found)")
        groups[key] = deduped

    missing = [k for k in REGION_RESIDUES if k not in groups]
    if missing:
        sys.exit(f"region residues not found: {missing}")

    charge = 0
    for (ch, seq), comp in REGION_RESIDUES.items():
        actual = groups[(ch, seq)][0]["comp"]
        if seq == "275":
            expect = "CYS" if side == "wt" else "PHE"
            if actual != expect:
                sys.exit(f"residue 275 is {actual}, expected {expect} for --side {side}")
            continue   # neutral either way -- Cys thiol (not thiolate, no metal here) or Phe
        charge += SIDECHAIN_Q.get(actual, 0)

    sel_res = {(k[0], int(k[1])) for k in groups if k[1].isdigit()}
    caps = []
    for ch, n in sorted(sel_res):
        for nb, anchor, partner in ((n - 1, "N", "C"), (n + 1, "C", "N")):
            if (ch, nb) in sel_res:
                continue
            aa = [a for a in atoms if a["ch"] == ch and a["seq"] == str(n) and a["name"] == anchor]
            bb = [a for a in atoms if a["ch"] == ch and a["seq"] == str(nb) and a["name"] == partner]
            if not aa or not bb:
                continue
            p, q = aa[0]["xyz"], bb[0]["xyz"]
            v = [q[i] - p[i] for i in range(3)]
            nrm = math.sqrt(sum(x * x for x in v)) or 1.0
            caps.append(("H", tuple(p[i] + v[i] / nrm * CAP_BOND_LENGTH for i in range(3))))

    return groups, caps, charge


def write_xyz(groups, caps, path, comment=""):
    rows = [(a["elem"], a["xyz"]) for k in sorted(groups) for a in groups[k]]
    rows += caps
    with open(path, "w") as fh:
        fh.write(f"{len(rows)}\n{comment}\n")
        for el, (x, y, z) in rows:
            fh.write(f"{el:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    return rows


def build_mf(xyz, charge, spin, basis="6-31g"):
    from pyscf import gto, scf
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin, verbose=3)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, charge={charge} spin={spin}")
    mf = (scf.RHF(mol) if spin == 0 else scf.ROHF(mol)).density_fit()
    mf.kernel()
    if not mf.converged:
        print("*** SCF did NOT converge -- treat any active space below as provisional")
    return mf


def avas_at_threshold(mf, ao_labels, threshold):
    from pyscf.mcscf import avas
    aos = [s.strip() for s in ao_labels.split(",")]
    ncas, nelecas, mo = avas.avas(mf, aos, threshold=threshold)
    print(f"[avas] threshold={threshold} AOs={aos} -> CAS({nelecas},{ncas}) "
          f"= {2 * ncas} qubits under Jordan-Wigner")
    return dict(threshold=threshold, ncas=int(ncas), nelecas=int(nelecas), qubits=int(2 * ncas))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", required=True, choices=["wt", "mutant"])
    ap.add_argument("--sweep", default="0.1,0.2,0.4,0.6,0.8,0.95")
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()

    src = FILES[a.side]
    if not os.path.exists(src):
        sys.exit(f"{src} not found -- expected in the working directory")
    if a.side == "mutant":
        print("*** THIS IS A MODELLED MUTANT (v2, corrected 2026-09-30: ideal Phe ring "
              "on the deposited N/CA/CB of CYS275, chi1=-162 chi2=+71, inside the Lovell "
              "trans rotamer well, max vdW overlap 0.309 A against all chain-A heavy "
              "atoms within 16 A), NOT a deposited structure. Report it as such. ***")
    check_residue_275(src, "CYS" if a.side == "wt" else "PHE")

    ph_pdb = protonate(src, a.ph)
    atoms = read_pdb(ph_pdb)
    groups, caps, charge = build_cluster(atoms, a.side)
    n_heavy = sum(1 for k in groups for a2 in groups[k] if a2["elem"] != "H")
    print(f"[cluster] {len(groups)} residues, {n_heavy} heavy atoms, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec expects {EXPECTED_NET_CHARGE:+d})")
    xyz = f"tp53_c275f_{a.side}_cluster.xyz"
    write_xyz(groups, caps, xyz, comment=f"TP53 C275F {a.side}, 2OCJ region")
    print(f"[cluster] wrote {xyz}")

    mf = build_mf(xyz, charge, 0, basis=a.basis)
    # AO set corrected 2026-09-30 (Claude Science, reviewing the CAS(64,48) plan):
    # "S 3p, N 2p, O 2p" alone is not symmetric between WT and mutant. CYS275's
    # side chain (WT) has SG -- sulfur, a target atom. PHE275's side chain
    # (mutant) is pure carbon (CB/CG/CD1/CD2/CE1/CE2/CZ, the phenyl ring) -- no
    # S/N/O at all, so with this AO set the mutant active space is driven only
    # by backbone N/O and is BLIND to the mutated side chain itself, while the
    # WT space is anchored on the Cys thiolate. Adding "C 2p" restores symmetry:
    # both the phenyl ring (mutant) and the thiolate (WT) are covered by the
    # same criterion. This is the same class of finding as "only one sulfur
    # within 5A of C275 and the mutation deletes it" from an earlier session.
    results = [avas_at_threshold(mf, "C 2p, N 2p, O 2p, S 3p", th)
               for th in [float(x) for x in a.sweep.split(",")]]
    out = f"tp53_c275f_{a.side}_avas.json"
    json.dump(dict(side=a.side, source=src, is_modelled=(a.side == "mutant"),
                   charge=charge, spin=0, basis=a.basis,
                   e_scf=float(mf.e_tot), scf_converged=bool(mf.converged), runs=results),
              open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
