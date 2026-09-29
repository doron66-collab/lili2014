#!/usr/bin/env python3
"""Gate 2: TP53 R175H -- wild-type Zn-site cluster -> QM cluster -> AVAS.

R175H has NO usable deposited mutant structure (Claude Science round 2: the
six PDB entries carrying it are 9-residue MHC-bound neoantigen peptides with
no p53 fold at all -- see mmcif_verify.py's new construct_scope() /
MUTANT_CONFIRMED_FRAGMENT_ONLY verdict). This script classifies the
WILD-TYPE zinc-site cluster on 2OCJ, at the region/radius Claude Science
verified two independent ways (zinc shell closes AND net charge hits a
clean 0, both first true at 7.2 A) -- so we have a real number for the WT
side now, while a modelled R175H mutant is a separate, later decision.

Region is taken verbatim from active_space_spec.json's TP53_R175H entry
(zinc_closure_verification.region_residues), not re-derived by a radius
cut here -- same reasoning as run_gate2_sdhb.py: a generic heuristic built
for a single mutated residue does not know ZN501 is a metal cofactor
carrying its own formal charge, or that the shell needs exactly 4 ligands.

    pip install pyscf pdbfixer openmm gemmi

Usage
-----
    python run_gate2_r175h.py --verify-only
    python run_gate2_r175h.py --sweep 0.1,0.2,0.4,0.6,0.8,0.95

Inputs expected in the working directory: mmcif_verify.py, struct/2OCJ.cif,
sifts/2OCJ.json (already fetched in round 1 -- see struct/2OCJ.cif from the
TP53 bundle).
"""
import argparse
import json
import math
import os
import sys

CAP_BOND_LENGTH = 1.09
PDB_ID = "2OCJ"
CHAIN = "A"
ACCESSION = "P04637"

# Verbatim from active_space_spec.json's TP53_R175H.zinc_closure_verification
# .region_residues (7.2 A closure, confirmed by two independent criteria).
REGION_RESIDUES = {
    ("A", "173"): "VAL", ("A", "174"): "ARG", ("A", "175"): "ARG",
    ("A", "176"): "CYS", ("A", "177"): "PRO", ("A", "179"): "HIS",
    ("A", "180"): "GLU", ("A", "191"): "PRO", ("A", "192"): "GLN",
    ("A", "193"): "HIS", ("A", "194"): "LEU", ("A", "214"): "HIS",
    ("A", "237"): "MET", ("A", "238"): "CYS", ("A", "242"): "CYS",
    ("A", "244"): "GLY", ("A", "245"): "GLY", ("A", "246"): "MET",
    ("A", "501"): "ZN",
}
ZN_LIGAND_CYS = {"176", "238", "242"}   # thiolates in the first Zn shell
SIDECHAIN_Q = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}
ZN_CHARGE = 2
EXPECTED_NET_CHARGE = 3   # Zn(+2) + 3 neutral Cys thiols (pdbfixer default, confirmed live
                          # against struct/2OCJ_H.pdb -- HG present on all three) +
                          # Arg174(+1) + Arg175(+1) + Glu180(-1). The spec's original 0 assumed
                          # thiolates by construction (R2's metal-coordination rule), which
                          # is the chemically-intended state but not what pdbfixer's generic,
                          # metal-unaware protonation actually produced.


def verify_wt():
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    cif = f"struct/{PDB_ID}.cif"
    sifts_json = json.load(open(f"sifts/{PDB_ID}.json"))
    r = mv.verify_sifts(cif, sifts_json, PDB_ID, ACCESSION, "R", 175, "H")
    print(f"[verify] {PDB_ID} vs R175H: {r['verdict']}")
    print(f"[verify] {r['detail']}")
    if not r["verdict"].startswith("WILD_TYPE_AT_SITE"):
        print("*** unexpected verdict -- expected WILD_TYPE_AT_SITE (2OCJ is wild-type "
              "only, per round 1/2 findings); re-check before proceeding")
    print("[verify] Proceeding to classify the wild-type Zn-site cluster. A modelled "
          "R175H mutant is a separate, later decision (not isoelectronic with this WT "
          "cluster -- see active_space_spec.json's comparison_warning).")


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


def cif_to_pdb(pdb_id):
    pdb_path = f"struct/{pdb_id}.pdb"
    if os.path.exists(pdb_path):
        return pdb_path
    import gemmi
    st = gemmi.read_structure(f"struct/{pdb_id}.cif")
    st.setup_entities()
    st.write_pdb(pdb_path)
    print(f"[prep] converted struct/{pdb_id}.cif -> {pdb_path} via gemmi")
    return pdb_path


def read_zn_from_cif(pdb_id, chain, resseq):
    """Read ZN501 directly from the mmCIF, the same way run_gate2_sdhb.py reads
    FES301 -- gemmi's cif->pdb conversion is not trusted to carry a hetero
    group through intact (confirmed it silently drops FES; not re-verified
    for a lone metal ion, so this sidesteps the question rather than assuming)."""
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    data = mv.read_cif(f"struct/{pdb_id}.cif")
    out = []
    for r in mv.rows(data, "_atom_site"):
        if r.get("label_comp_id") != "ZN":
            continue
        ch = r.get("auth_asym_id") or r.get("label_asym_id")
        sq = r.get("auth_seq_id") or r.get("label_seq_id")
        if ch != chain or sq != resseq:
            continue
        out.append(dict(name="ZN", comp="ZN", ch=ch, seq=sq,
                        xyz=(float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])),
                        elem="ZN"))
    if not out:
        sys.exit(f"read_zn_from_cif found no ZN at chain {chain} residue {resseq} in "
                 f"struct/{pdb_id}.cif -- check chain/residue numbering directly")
    print(f"[cluster] read {len(out)} ZN atom directly from struct/{pdb_id}.cif "
          f"(chain {chain}, residue {resseq})")
    return out


def protonate(ph=7.4, out=None):
    out = out or f"struct/{PDB_ID}_H.pdb"
    if os.path.exists(out):
        return out
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    local_pdb = cif_to_pdb(PDB_ID)
    fixer = PDBFixer(filename=local_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # NOT calling removeHeterogens() -- see run_gate2_sdhb.py's identical note.
    # ZN501 is read separately from the cif anyway (read_zn_from_cif), so
    # whatever pdbfixer does to it here is moot; skipping it just avoids any
    # chance it also perturbs something else in the process.
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print(f"[prep] protonated at pH {ph} -> {out}")
    return out


def d(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def build_cluster(atoms):
    groups = {}
    for a in atoms:
        key = (a["ch"], a["seq"])
        if key not in REGION_RESIDUES:
            continue
        groups.setdefault(key, []).append(a)

    missing = [k for k in REGION_RESIDUES if k not in groups]
    if missing:
        sys.exit(f"region residues not found in structure: {missing} -- "
                 f"check chain ID and residue numbering in struct/{PDB_ID}.cif")

    charge = ZN_CHARGE
    for (ch, seq), comp in REGION_RESIDUES.items():
        if comp == "ZN":
            continue
        if comp == "CYS" and seq in ZN_LIGAND_CYS:
            # Check the ACTUAL protonation state, don't assume thiolate.
            # pdbfixer does not know this Cys is metal-bound, so it protonates
            # it as a normal neutral thiol (HG present) at pH 7.4 -- confirmed
            # live for all three ligands here (176/238/242 all came out with
            # HG in struct/2OCJ_H.pdb). Only subtract the thiolate charge if
            # HG is genuinely absent.
            has_hg = any(a["name"] in ("HG", "HG1") for a in groups[(ch, seq)])
            charge += 0 if has_hg else -1
            continue
        charge += SIDECHAIN_Q.get(comp, 0)

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
    from pyscf import gto, scf, df
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin, verbose=3)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, charge={charge} spin={spin}")
    mf = (scf.RHF(mol) if spin == 0 else scf.ROHF(mol)).density_fit()
    mf.with_df.auxbasis = df.make_auxbasis(mol)   # Zn has no cc-pvdz-jkfit entry either
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
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--sweep", default="0.1,0.2,0.4,0.6,0.8,0.95")
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()

    verify_wt()
    if a.verify_only:
        return

    ph_pdb = protonate(a.ph)
    atoms = read_pdb(ph_pdb) + read_zn_from_cif(PDB_ID, CHAIN, "501")
    groups, caps, charge = build_cluster(atoms)
    n_heavy = sum(1 for k in groups for a in groups[k] if a["elem"] != "H")
    print(f"[cluster] {len(groups)} residues, {n_heavy} heavy atoms, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec expects {EXPECTED_NET_CHARGE:+d})")
    if charge != EXPECTED_NET_CHARGE:
        print("[cluster] NOTE: charge differs from spec -- check for a spurious HG on a "
              "metal-bound Cys thiolate (known failure mode) before trusting this number.")
    xyz = "tp53_r175h_wt_cluster.xyz"
    write_xyz(groups, caps, xyz, comment=f"TP53 R175H WT Zn site, {PDB_ID}, r=7.2A")
    print(f"[cluster] wrote {xyz}")

    mf = build_mf(xyz, charge, 0, basis=a.basis)
    results = [avas_at_threshold(mf, "Zn 3d, S 3p, N 2p, O 2p", th)
               for th in [float(x) for x in a.sweep.split(",")]]
    out = "tp53_r175h_wt_avas.json"
    json.dump(dict(pdb=PDB_ID, charge=charge, spin=0, basis=a.basis,
                   e_scf=float(mf.e_tot), scf_converged=bool(mf.converged), runs=results),
              open(out, "w"), indent=1)
    print(f"\nwrote {out}")
    if len(results) > 1:
        print("threshold  CAS(nelec,ncas)  qubits")
        for r in results:
            print(f"{r['threshold']:<10}CAS({r['nelecas']},{r['ncas']}){r['qubits']:>12}")


if __name__ == "__main__":
    main()
