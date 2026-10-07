#!/usr/bin/env python3
"""Gate 2: SDHB [3Fe-4S] (S3, Cys196/243/249) wild-type cluster -> QM cluster -> AVAS.

Sibling to run_gate2_sdhb.py (S1, [2Fe-2S]) -- NOT a replacement for it. Claude
Science's 2026-10-07 consult retracted the S1 [2Fe-2S] recommendation from the
round before: every active-space construction of S1 is either small enough that
DMRG solves it numerically exactly (dE(M)=0, R undefined) or, at its largest
(Fe 3d + all S 3p, CAS(46,28)), 82% filled -- the same orbital-starvation defect
that invalidated the original eight Class B verdicts. S1 is kept in the sibling
script ONLY as a solver-validation case now (CAS(10,10) should be FCI-exact;
DMRG and SHCI must agree with it to ~1e-9), not as a Class A candidate.

S3 [3Fe-4S] is the actual candidate this round, for reasons verified from
PDB 8GS8's own coordinates (not estimated): the target active space (Fe 3d +
bridging S 3p = 15+12 = 27 AOs) is pool==cut against the cluster's own 27
orbitals (no selection-pressure confound); the exact Schmidt rank (~2.59e6) is
four orders of magnitude past what M=4000 captures (0.155%), so dE(M) is
actually measurable and R is well-defined -- the only SDHB construction where
that's true; and the deposited geometry matches 24 matched-ligation X-ray
references (<=1.2 A) to 0.1 sigma (Fe-Fe +0.005 A), usable as deposited with NO
relaxation needed (S1 and S2 both failed this check: 13.5 sigma and 4.6 sigma
respectively, and are NOT usable as deposited -- see sdhb_cluster_geometry_8gs8.csv).

KNOWN LIGATION BUG TO AVOID: ligands are NOT assignable by sequential cysteine
motif number. Cys253 sits in the protein's third sequential Cys motif but
coordinates the [4Fe-4S] (S2) cluster; Cys196 sits in the second motif but
coordinates THIS [3Fe-4S] (S3) cluster. REGION_RESIDUES below is the geometric
assignment (measured Fe-S distance < 2.8 A), not a motif-order guess.

Ground state is a HALF-INTEGER spin (S=1/2, 15 d-electrons: 3x Fe(III) d5,
high-spin-coupled then antiferromagnetically combined) -- odd electron count,
so there IS no S=0 state here; see build_lowspin_casci() below, generalized
from run_gate2_sdhb.py's build_s0_casci() for an odd na-nb difference of 1
instead of 0. The high-spin control is S=15/2 (spin=15, fully ferromagnetic,
single ROHF determinant).

    pip install pyscf pdbfixer openmm   # same as run_gate2_sdhb.py

Usage
-----
    python run_gate2_sdhb_s3.py --verify-only
    python run_gate2_sdhb_s3.py                # builds both spin states

Inputs expected in the working directory: mmcif_verify.py (same directory as
run_gate2_sdhb.py); struct/8GS8.cif and sifts/8GS8.json are auto-fetched.
"""
import argparse
import json
import math
import os
import sys

CAP_BOND_LENGTH = 1.09
PDB_ID = "8GS8"
CHAIN = "B"
ACCESSION = "P21912"   # same SDHB protein as S1/S2, different PDB entry

# Geometric assignment (measured Fe-S < 2.8 A in 8GS8), NOT motif-order --
# see module docstring's "KNOWN LIGATION BUG" warning. F3S303 is the [3Fe-4S]
# cofactor itself.
REGION_RESIDUES = {
    ("B", "196"): "CYS", ("B", "243"): "CYS", ("B", "249"): "CYS",
    ("B", "303"): "F3S",
}
CYS_LIGANDS = {"196", "243", "249"}   # thiolates, part of the cluster core
SIDECHAIN_Q = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}
# [3Fe-4S]1+ core: 3 Fe(III) + 4 inorganic S(2-) = 3(+3) + 4(-2) = +1
# (contrast run_gate2_sdhb.py's FES_CORE_CHARGE=+2 for [2Fe-2S]2+ -- a
# genuinely different formal core charge, not a typo).
F3S_CORE_CHARGE = 1
EXPECTED_NET_CHARGE = -2   # F3S core +1, three Cys thiolates -3


def fetch(pdb_id, cache="struct"):
    """Same pattern as run_gate2_sdhb.py / run_gate2_avas.py."""
    import urllib.request
    os.makedirs(cache, exist_ok=True)
    cif = os.path.join(cache, f"{pdb_id}.cif")
    if not os.path.exists(cif):
        print(f"[fetch] downloading {cif} from RCSB...")
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id}.cif", cif)
    sif = os.path.join("sifts", f"{pdb_id}.json")
    os.makedirs("sifts", exist_ok=True)
    if not os.path.exists(sif):
        print(f"[fetch] downloading {sif} from EBI SIFTS...")
        with open(sif, "wb") as fh:
            fh.write(urllib.request.urlopen(
                f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"
                f"{pdb_id.lower()}", timeout=60).read())
    return cif, json.load(open(sif))


def verify_wt():
    """Checks the clinically-documented C196Y (PPGL4, rs876658367) variant
    position is wild-type (CYS) in 8GS8 -- same discipline as
    run_gate2_sdhb.py's C101Y/C101S checks for S1, adapted to S3's own
    clinically-anchored ligand position (see module docstring)."""
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    cif, sifts_json = fetch(PDB_ID)
    for mut1, label in (("Y", "C196Y"),):
        r = mv.verify_sifts(cif, sifts_json, PDB_ID, ACCESSION, "C", 196, mut1)
        print(f"[verify] {PDB_ID} vs {label}: {r['verdict']}")
        print(f"[verify] {r['detail']}")
        if not r["verdict"].startswith("WILD_TYPE_AT_SITE"):
            print(f"*** unexpected verdict for {label} -- re-check before proceeding "
                  f"(expected WILD_TYPE_AT_SITE, per Claude Science's report)")
    print("[verify] Proceeding to classify the wild-type [3Fe-4S] (S3) cluster.")


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


def read_f3s_from_cif(pdb_id, chain, resseq):
    """Read F3S303's own atoms directly from the mmCIF's _atom_site loop --
    same reasoning as run_gate2_sdhb.py's read_fes_from_cif(): do not rely on
    gemmi's cif->pdb conversion to preserve (or drop) a compound cofactor
    HETATM group either way. Read it directly, and run_gate2_sdhb_s3.py's own
    main() below explicitly excludes any F3S read from the protonated PDB to
    avoid the exact double-counting bug found live in the S1 script
    2026-10-07 ("RuntimeError: Ill geometry")."""
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    data = mv.read_cif(f"struct/{pdb_id}.cif")
    out = []
    for r in mv.rows(data, "_atom_site"):
        if r.get("label_comp_id") != "F3S":
            continue
        ch = r.get("auth_asym_id") or r.get("label_asym_id")
        sq = r.get("auth_seq_id") or r.get("label_seq_id")
        if ch != chain or sq != resseq:
            continue
        out.append(dict(
            name=r.get("label_atom_id", "?"), comp="F3S", ch=ch, seq=sq,
            xyz=(float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])),
            elem=(r.get("type_symbol") or r.get("label_atom_id", "?")[0]).upper()))
    if not out:
        sys.exit(f"read_f3s_from_cif found no F3S atoms at chain {chain} residue {resseq} "
                 f"in struct/{pdb_id}.cif -- check chain/residue numbering directly")
    print(f"[cluster] read {len(out)} F3S atoms directly from struct/{pdb_id}.cif "
          f"(chain {chain}, residue {resseq}) -- not from the protonated PDB")
    return out


def build_cluster(atoms):
    """Same approach as run_gate2_sdhb.py's build_cluster(): explicit charge
    assignment from the cofactor's own formal core charge + thiolates, not a
    generic metal-distance heuristic."""
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

    charge = F3S_CORE_CHARGE
    for (ch, seq), comp in REGION_RESIDUES.items():
        if comp == "F3S":
            continue
        if comp == "CYS" and seq in CYS_LIGANDS:
            charge += -1
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
    # Deliberately NOT calling fixer.removeHeterogens() -- same reasoning as
    # run_gate2_sdhb.py: it can silently strip compound cofactors like F3S.
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print(f"[prep] protonated at pH {ph} -> {out}")
    return out


def build_mf(xyz, charge, spin, basis="def2-tzvp"):
    """High-spin (S=15/2, spin=15) ROHF reference -- the single-determinant,
    genuinely well-behaved state. Same reasoning as run_gate2_sdhb.py's
    build_mf(): never build the low-spin (S=1/2) state's own independent
    mean-field, since it's open-shell/multi-reference despite having a
    well-defined spin quantum number."""
    from pyscf import gto, scf
    if spin != 15:
        raise ValueError("build_mf() only builds the high-spin (S=15/2, spin=15) reference -- "
                          "see this function's docstring. The S=1/2 state is evaluated by "
                          "build_lowspin_casci() on these same orbitals.")
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin, verbose=3)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, charge={charge} spin={spin}")
    mf = scf.ROHF(mol).density_fit()
    # Same disk-space and auxbasis fixes as run_gate2_sdhb.py, applied here too --
    # this cluster is smaller (27 target orbitals vs 28) but the full QM cluster
    # (all heavy+capping atoms, not just the active space) is comparable in size,
    # so the same ERI-cache blowup risk applies.
    from pyscf import lib
    dftmp_dir = os.path.join(os.getcwd(), "pyscf_tmp")
    os.makedirs(dftmp_dir, exist_ok=True)
    lib.param.TMPDIR = dftmp_dir
    print(f"[scratch] PySCF density-fitting scratch redirected to {dftmp_dir}")
    mf.with_df.auxbasis = "def2-universal-jkfit"
    mf.kernel()
    if not mf.converged:
        print("*** SCF did NOT converge -- treat any active space below as provisional")
    ss, mult = mf.spin_square()
    print(f"[spin] high-spin reference: <S^2>={ss:.4f} (expect 63.75 for S=15/2), "
          f"2S+1={mult:.4f} (expect 16.0)")
    if abs(ss - 63.75) > 0.5:
        print("*** <S^2> does not match the expected S=15/2 value -- reference may have "
              "converged to a different spin state than intended; check before trusting "
              "any orbitals derived from it")
    return mf


def avas_at_threshold(mf, ao_labels, threshold):
    from pyscf.mcscf import avas
    aos = [s.strip() for s in ao_labels.split(",")]
    ncas, nelecas, mo = avas.avas(mf, aos, threshold=threshold)
    print(f"[avas] threshold={threshold} AOs={aos} -> CAS({nelecas},{ncas}) "
          f"= {2 * ncas} qubits under Jordan-Wigner")
    spec = dict(threshold=threshold, ncas=int(ncas), nelecas=int(nelecas), qubits=int(2 * ncas))
    return spec, mo


def build_lowspin_casci(mf_highspin, mo_coeff, charge, ncas, nelecas):
    """S=1/2 (ground, antiferromagnetically coupled) state, CASCI on the
    HIGH-SPIN reference's own orbitals. Generalized from run_gate2_sdhb.py's
    build_s0_casci() for an ODD na-nb difference (=1, not =0): S3 has 15
    d-electrons (odd), so there is no S=0 state here at all -- the true
    ground state carries a half-integer spin quantum number by construction,
    not by choice.
    """
    from pyscf import mcscf
    mol_ls = mf_highspin.mol.copy()
    mol_ls.spin = 1
    mol_ls.charge = charge
    mol_ls.build()
    assert nelecas % 2 == 1, (
        f"nelecas={nelecas} is even -- S=1/2 (odd-electron-count) CASCI is not "
        f"well-defined on this active space; re-check the AVAS selection before proceeding "
        f"(expected 15 d-electrons for [3Fe-4S], an odd number)")
    na = (nelecas + 1) // 2
    nb = (nelecas - 1) // 2
    mc = mcscf.CASCI(mf_highspin, ncas, (na, nb))
    mc.mol = mol_ls
    mc.mo_coeff = mo_coeff
    mc.kernel()
    # Same DFCASCI.spin_square() AttributeError fix as run_gate2_sdhb.py's
    # build_s0_casci() -- go through fcisolver directly.
    ss, mult = mc.fcisolver.spin_square(mc.ci, mc.ncas, mc.nelecas)
    print(f"[spin] S=1/2 CASCI: <S^2>={ss:.4f} (expect 0.75), 2S+1={mult:.4f} (expect 2.0)")
    if abs(ss - 0.75) > 0.2:
        print("*** <S^2> does not match the expected S=1/2 value for this CASCI solution -- "
              "the solver may have found a different spin state within the same orbital "
              "space; check the CI vector / active-space symmetry before trusting the energy")
    return mc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--threshold", type=float, default=0.2,
                     help="AVAS threshold on the high-spin reference. Fe 3d (15 AOs, 3 Fe) + "
                          "bridging S 3p (12 AOs) = 27 target AOs vs 27 cluster orbitals -- "
                          "pool==cut exactly, per Claude Science's 2026-10-07 consult; nothing "
                          "for a threshold sweep to select between.")
    ap.add_argument("--basis", default="def2-tzvp",
                     help="def2-tzvp (Claude Science's preferred choice, not just the def2-svp "
                          "minimum) -- S3 is the actual candidate this round, worth the extra "
                          "cost. R is computed within SDHB against its own matched control, so "
                          "this does not need to match any other target's basis.")
    ap.add_argument("--fe-semicore", choices=["core", "active"], default="core",
                     help="Same explicit decision point as run_gate2_sdhb.py -- see that "
                          "script's own flag for the full rationale.")
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()

    verify_wt()
    if a.verify_only:
        return

    ph_pdb = protonate(a.ph)
    # Same double-counting guard as run_gate2_sdhb.py's fix for the FES case --
    # exclude any F3S read from the protonated PDB; read_f3s_from_cif()'s direct
    # cif read is the only source of those atoms.
    all_pdb_atoms = read_pdb(ph_pdb)
    pdb_atoms = [a for a in all_pdb_atoms if a["comp"] != "F3S"]
    dropped = len(all_pdb_atoms) - len(pdb_atoms)
    if dropped:
        print(f"[cluster] NOTE: {dropped} F3S atom(s) were present in the protonated PDB -- "
              f"excluded; read_f3s_from_cif()'s direct cif read is the single source of truth.")
    atoms = pdb_atoms + read_f3s_from_cif(PDB_ID, CHAIN, "303")

    # Root cause of the 2026-10-07 "Electron number 312 and spin 15 are not
    # consistent" crash, per Claude Science's consult reply: pdbfixer's
    # addMissingHydrogens() protonates every CYS as a neutral thiol (-SH) by
    # default -- it has no notion of a metal-coordinating thiolate (-S-)
    # unless the residue is explicitly typed CYM, which nothing here does.
    # The three ligating cysteines therefore each carried one spurious HG on
    # their SG. EXPECTED_NET_CHARGE above could never catch this: it's the
    # formal recipe (core charge + per-residue contributions), compared
    # against a `charge` variable computed the same way -- a constant
    # checked against itself, not against the real atom list.
    spurious_hg = [a for a in atoms
                   if a["comp"] == "CYS" and a["seq"] in CYS_LIGANDS and a["name"] == "HG"]
    if spurious_hg:
        print(f"[cluster] NOTE: stripping {len(spurious_hg)} spurious thiol H (HG) from "
              f"ligating Cys residue(s) {sorted(a['seq'] for a in spurious_hg)} -- these "
              f"coordinate the [3Fe-4S] cluster as thiolates and must not carry that "
              f"hydrogen (see cluster_spin_guard.py).")
        atoms = [a for a in atoms
                 if not (a["comp"] == "CYS" and a["seq"] in CYS_LIGANDS and a["name"] == "HG")]
    groups, caps, charge = build_cluster(atoms)
    n_heavy = sum(1 for k in groups for a in groups[k] if a["elem"] != "H")
    print(f"[cluster] {len(groups)} residues, {n_heavy} heavy atoms, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec expects {EXPECTED_NET_CHARGE:+d})")
    if charge != EXPECTED_NET_CHARGE:
        print(f"[cluster] NOTE: charge differs from spec -- check whether pdbfixer assigned "
              f"a spurious HG to a metal-bound Cys thiolate before trusting this number.")
    xyz = "sdhb_s3_cluster.xyz"
    write_xyz(groups, caps, xyz, comment=f"SDHB [3Fe-4S] S3 site, {PDB_ID}, Fe3d+bridging-S3p AVAS")
    print(f"[cluster] wrote {xyz}")
    print(f"[core] Fe 3s/3p treated as {a.fe_semicore} (--fe-semicore)")

    # Preflight -- run before gto.M(), not after the SCF. Same atom-list-
    # ordering `write_xyz` itself uses (sorted(groups) then each group's own
    # list order, caps appended last), so sg_indices line up with the XYZ
    # file's own atom order.
    import cluster_spin_guard as csg
    ordered = [atm for k in sorted(groups) for atm in groups[k]]
    sg_indices = [i for i, atm in enumerate(ordered) if atm["name"] == "SG"]
    symbols = [atm["elem"] for atm in ordered] + ["H"] * len(caps)
    coords = [atm["xyz"] for atm in ordered] + [c[1] for c in caps]
    nelectron = csg.sum_nuclear_charge(symbols) - charge
    csg.preflight(symbols, coords, nelectron, charge, spin=15,
                  sg_indices=sg_indices, iron_oxidation_states=[3, 3, 3],
                  core_charge=F3S_CORE_CHARGE, label="SDHB S3")
    print(f"[guard] preflight passed: N={nelectron}, charge={charge:+d}, spin=15, "
          f"{len(sg_indices)} ligating SG checked for residual H")

    mf_hs = build_mf(xyz, charge, spin=15, basis=a.basis)
    spec, mo = avas_at_threshold(mf_hs, "Fe 3d, S 3p", a.threshold)
    ncas, nelecas = spec["ncas"], spec["nelecas"]

    mo_path = f"{a.out_dir}/sdhb_s3_mo_coeff_highspin.npy"
    import numpy as np
    np.save(mo_path, mo)
    print(f"[avas] saved shared orbitals -> {mo_path}")

    mc_ls = build_lowspin_casci(mf_hs, mo, charge, ncas, nelecas)

    out = f"{a.out_dir}/sdhb_s3_both_spins.json"
    json.dump(dict(
        pdb=PDB_ID, charge=charge, basis=a.basis, fe_semicore=a.fe_semicore,
        threshold=a.threshold, ncas=ncas, nelecas=nelecas, qubits=spec["qubits"],
        high_spin=dict(spin=15, e_scf=float(mf_hs.e_tot), scf_converged=bool(mf_hs.converged),
                        spin_square_expected=63.75),
        casci_lowspin=dict(spin=1, e_tot=float(mc_ls.e_tot), spin_square_expected=0.75),
        mo_coeff_path=mo_path,
    ), open(out, "w"), indent=1)
    print(f"\nwrote {out}")
    print(f"\n[next] P3 FIRST (not this target): validate DMRG/SHCI against exact FCI on "
          f"S1's CAS(10,10) -- cheap, must agree to ~1e-9 before anything here is readable. "
          f"Then feed ncas={ncas}, nelecas={nelecas}, --load-orbitals {mo_path} into "
          f"solange_dmrg.py's --geometry path for the bond-dimension ladder (both spin "
          f"states), plus a size/filling-matched chemistry-free control at CAS({nelecas},{ncas}), "
          f"same basis, before reading any R_spin value. Pre-registered predictions P1/P2 "
          f"(R_spin > 1 for S3 and S2; R ordering S2 > S3 > S1 by iron count) are falsifiable -- "
          f"report the result either way.")


if __name__ == "__main__":
    main()
