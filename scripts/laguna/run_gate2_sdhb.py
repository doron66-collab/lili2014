#!/usr/bin/env python3
"""Gate 2: SDHB [2Fe-2S] (S1, Cys101) wild-type cluster -> QM cluster -> AVAS.

Companion to run_gate2_avas.py, for a genuinely different case: here we
classify the WILD-TYPE metal-cluster site (no mutant structure exists or is
wanted -- see active_space_spec_sdhb.json's fe_s_cluster / mutant_hypothesis
fields). The region is NOT rediscovered by a radius cut here -- it is taken
verbatim from Claude Science's verified spec (region_residues, R1/R2/R4
reasoning already done and checked against the structure), because for a
4-ligand metal cluster a generic distance-based charge/closure heuristic
(built for a single mutated residue) does not know FES is a metal cofactor
carrying its own formal charge. Getting that wrong here would silently
reintroduce exactly the kind of error DP1 exists to catch.

    pip install pyscf pdbfixer openmm   # same as run_gate2_avas.py

Usage
-----
    # verify WT status only (does not build/run anything)
    python run_gate2_sdhb.py --verify-only

    # build the cluster and run one spin state
    python run_gate2_sdhb.py --spin 0     # S=0, antiferromagnetic ground state
    python run_gate2_sdhb.py --spin 10    # S=5, ferromagnetic (isoelectronic)

Inputs expected in the working directory: mmcif_verify.py, struct/9KC4.cif,
sifts/9KC4.json (fetched separately -- see chat).
"""
import argparse
import json
import math
import os
import sys

CAP_BOND_LENGTH = 1.09
PDB_ID = "9KC4"
CHAIN = "B"
ACCESSION = "P21912"

# Verbatim from active_space_spec_sdhb.json's SDHB_C101Y.region_residues /
# region_rule -- not re-derived. FES301 is the [2Fe-2S] cofactor itself.
REGION_RESIDUES = {
    ("B", "92"): "SER", ("B", "93"): "CYS", ("B", "94"): "ARG",
    ("B", "95"): "GLU", ("B", "96"): "GLY", ("B", "97"): "ILE",
    ("B", "98"): "CYS", ("B", "99"): "GLY", ("B", "100"): "SER",
    ("B", "101"): "CYS", ("B", "113"): "CYS", ("B", "301"): "FES",
}
CYS_LIGANDS = {"93", "98", "101", "113"}   # thiolates, part of the cluster core
SIDECHAIN_Q = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}
FES_CORE_CHARGE = 2   # [Fe2S2]2+ oxidised core, per spec's R4 assignment
EXPECTED_NET_CHARGE = -2   # FES core +2, four Cys thiolates -4, ARG94 +1, GLU95 -1


def verify_wt():
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    cif = f"struct/{PDB_ID}.cif"
    sifts_json = json.load(open(f"sifts/{PDB_ID}.json"))
    for mut1, label in (("Y", "C101Y"), ("S", "C101S")):
        r = mv.verify_sifts(cif, sifts_json, PDB_ID, ACCESSION, "C", 101, mut1)
        print(f"[verify] {PDB_ID} vs {label}: {r['verdict']}")
        print(f"[verify] {r['detail']}")
        if not r["verdict"].startswith("WILD_TYPE_AT_SITE"):
            print(f"*** unexpected verdict for {label} -- re-check before proceeding "
                  f"(expected WILD_TYPE_AT_SITE, per Claude Science's report)")
    print("[verify] Proceeding to classify the wild-type cluster, per gate2.py's "
          "own worked example and the request's explicit scope (§ 'What NOT to do').")


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


def d(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def read_fes_from_cif(pdb_id, chain, resseq):
    """Read FES301's own atoms directly from the mmCIF's _atom_site loop.

    Needed because gemmi's write_pdb() (cif_to_pdb, above) silently drops the
    FES HETATM records during cif->pdb conversion -- confirmed by grepping
    struct/9KC4.pdb for "FES" and finding nothing, though the cif's own HET/
    LINK records show it at chain B, residue 301. No protonation step is
    needed for these atoms (Fe/S do not carry added hydrogens), so reading
    them straight from the original, unmodified cif is not just a workaround
    -- it is the more correct source for this specific group."""
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    data = mv.read_cif(f"struct/{pdb_id}.cif")
    out = []
    for r in mv.rows(data, "_atom_site"):
        if r.get("label_comp_id") != "FES":
            continue
        ch = r.get("auth_asym_id") or r.get("label_asym_id")
        sq = r.get("auth_seq_id") or r.get("label_seq_id")
        if ch != chain or sq != resseq:
            continue
        out.append(dict(
            name=r.get("label_atom_id", "?"), comp="FES", ch=ch, seq=sq,
            xyz=(float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])),
            elem=(r.get("type_symbol") or r.get("label_atom_id", "?")[0]).upper()))
    if not out:
        sys.exit(f"read_fes_from_cif found no FES atoms at chain {chain} residue {resseq} "
                 f"in struct/{pdb_id}.cif -- check chain/residue numbering directly")
    print(f"[cluster] read {len(out)} FES atoms directly from struct/{pdb_id}.cif "
          f"(chain {chain}, residue {resseq}) -- not from the protonated PDB")
    return out


def build_cluster(atoms):
    """Extract exactly REGION_RESIDUES (fixed set, see module docstring),
    cap cut peptide bonds, and compute net charge from explicit assignment
    (FES core + thiolates + standard sidechain charges) -- not a generic
    metal-distance heuristic, since FES is a named cofactor, not an element."""
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

    charge = FES_CORE_CHARGE
    for (ch, seq), comp in REGION_RESIDUES.items():
        if comp == "FES":
            continue
        if comp == "CYS" and seq in CYS_LIGANDS:
            charge += -1   # metal-bound thiolate, by construction (spec R4)
            continue
        charge += SIDECHAIN_Q.get(comp, 0)

    # cap cut peptide bonds (same scheme as run_gate2_avas.py)
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


def cif_to_pdb(pdb_id):
    """9KC4 is a large cryo-EM entry -- RCSB has no legacy .pdb for it, so
    PDBFixer(pdbid=...)'s automatic fetch 404s. Convert the mmCIF we already
    have locally instead (same fix as run_gate2_avas.py needed)."""
    pdb_path = f"struct/{pdb_id}.pdb"
    if os.path.exists(pdb_path):
        return pdb_path
    import gemmi
    st = gemmi.read_structure(f"struct/{pdb_id}.cif")
    st.setup_entities()
    st.write_pdb(pdb_path)
    print(f"[prep] converted struct/{pdb_id}.cif -> {pdb_path} via gemmi")
    return pdb_path


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
    # NOTE: deliberately NOT calling fixer.removeHeterogens() here. It strips
    # multi-atom hetero groups (it keeps lone metal ions like Zn, but treats a
    # compound cofactor like FES as a removable "ligand") -- which silently
    # deleted the [2Fe-2S] cluster itself in an earlier run. Any leftover
    # waters are harmless: build_cluster() below only extracts atoms whose
    # (chain, resseq) is in REGION_RESIDUES, so stray waters are just ignored.
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print(f"[prep] protonated at pH {ph} -> {out}")
    return out


def build_mf(xyz, charge, spin, basis="ccpvdz"):
    from pyscf import gto, scf, df
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin, verbose=3)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, charge={charge} spin={spin}")
    mf = (scf.RHF(mol) if spin == 0 else scf.ROHF(mol)).density_fit()
    # The default JK-fit auxiliary basis (cc-pvdz-jkfit) has no entry for Fe --
    # transition metals are outside standard Dunning aux-basis coverage. Build
    # an even-tempered (ETB) auxiliary basis instead, which covers every
    # element actually present rather than requiring a hand-picked one.
    mf.with_df.auxbasis = df.make_auxbasis(mol)
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
    ap.add_argument("--spin", type=int, choices=[0, 10], help="pyscf spin_2S: 0 or 10 (see spec)")
    ap.add_argument("--sweep", default="0.1,0.2,0.4,0.6,0.8,0.95")
    ap.add_argument("--basis", default="ccpvdz",
                     help="Fe-S clusters need polarization on S/Fe; 6-31g is too thin for "
                          "the cluster core (spec did not pin a basis -- flag this choice)")
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()

    verify_wt()
    if a.verify_only:
        return
    if a.spin is None:
        sys.exit("--spin 0 or --spin 10 required (see active_space_spec_sdhb.json "
                 "candidate_spins) unless --verify-only")

    ph_pdb = protonate(a.ph)
    atoms = read_pdb(ph_pdb) + read_fes_from_cif(PDB_ID, CHAIN, "301")
    groups, caps, charge = build_cluster(atoms)
    n_heavy = sum(1 for k in groups for a in groups[k] if a["elem"] != "H")
    print(f"[cluster] {len(groups)} residues, {n_heavy} heavy atoms, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec expects {EXPECTED_NET_CHARGE:+d})")
    if charge != EXPECTED_NET_CHARGE:
        print(f"[cluster] NOTE: charge differs from spec -- check whether pdbfixer assigned "
              f"a spurious HG to a metal-bound Cys thiolate (known failure mode; see "
              f"run_gate2_avas.py's build_cluster comments) before trusting this number.")
    xyz = "sdhb_c101_cluster.xyz"
    write_xyz(groups, caps, xyz, comment=f"SDHB [2Fe-2S] S1 site, {PDB_ID}, spin={a.spin}")
    print(f"[cluster] wrote {xyz}")

    mf = build_mf(xyz, charge, a.spin, basis=a.basis)
    results = [avas_at_threshold(mf, "Fe 3d, S 3p", th) for th in [float(x) for x in a.sweep.split(",")]]
    out = f"sdhb_c101_spin{a.spin}_avas.json"
    json.dump(dict(pdb=PDB_ID, charge=charge, spin=a.spin, basis=a.basis,
                   e_scf=float(mf.e_tot), scf_converged=bool(mf.converged), runs=results),
              open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
