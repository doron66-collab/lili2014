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


def fetch(pdb_id, cache="struct"):
    """Same pattern as run_gate2_avas.py's fetch() -- added here because this
    script previously assumed struct/9KC4.cif and sifts/9KC4.json were already
    present from an earlier manual fetch, which is exactly what failed live
    2026-10-07 (FileNotFoundError on a fresh checkout)."""
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
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    cif, sifts_json = fetch(PDB_ID)
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
    """Mean-field reference, HIGH-SPIN ONLY (spin=10, S=5).

    Claude Science's 2026-10-07 consult flagged the prior version of this
    function as the pipeline's actual blocking defect: it built an
    independent RHF reference for the spin=0 (antiferromagnetic) state.
    [2Fe-2S]2+ at S=0 is Fe(III)/Fe(III) antiferromagnetically coupled --
    MULTI-REFERENCE even though formally closed-shell -- so a single-
    determinant RHF reference for that state is not a valid starting point
    (not merely low-quality: the wrong kind of wavefunction for what's being
    asked of it). Standard practice for AF-coupled metal clusters (see e.g.
    Noodleman broken-symmetry DFT) is instead to build ONE well-behaved
    single-determinant reference on the HIGH-SPIN state (S=5 here, ferro-
    magnetic, spin=10, genuinely single-reference), run AVAS on THAT
    reference, and reuse the resulting orbitals for CASCI/DMRG on BOTH spin
    states -- never build S=0's own independent mean-field at all. This
    function therefore only ever runs at spin=10 now; see build_s0_casci()
    below for how the S=0 state is actually evaluated.
    """
    from pyscf import gto, scf, df
    if spin != 10:
        raise ValueError("build_mf() only builds the high-spin (S=5, spin=10) reference now -- "
                          "see this function's docstring. The spin=0 state is evaluated by "
                          "build_s0_casci() on these same orbitals, not by its own mean-field.")
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin, verbose=3)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, charge={charge} spin={spin}")
    mf = scf.ROHF(mol).density_fit()
    # Found live 2026-10-07: at 150 atoms / 1570 basis functions, the density-fitting
    # Cholesky ERI cache hit ~24GB and filled the SLURM job's node-local /tmp (typically
    # small), crashing with OSError "No space left on device" partway through writing.
    # Redirect PySCF's scratch to the CURRENT working directory instead -- wherever this
    # script is actually being run from (beegfs or similar networked storage on Laguna),
    # which has far more headroom than node-local /tmp. Must be set before density_fit()
    # actually builds its scratch file, i.e. before mf.kernel() below.
    from pyscf import lib
    dftmp_dir = os.path.join(os.getcwd(), "pyscf_tmp")
    os.makedirs(dftmp_dir, exist_ok=True)
    lib.param.TMPDIR = dftmp_dir
    print(f"[scratch] PySCF density-fitting scratch redirected to {dftmp_dir} "
          f"(was defaulting to node-local /tmp, too small for this system's ~24GB+ "
          f"ERI cache) -- clean this directory up after the run, it is not auto-deleted")
    # Also switch the auxiliary basis from an auto-generated even-tempered (ETB) one to
    # the standard def2-universal-jkfit, which properly covers Fe/S (unlike cc-pvdz-jkfit)
    # without PySCF having to improvise one -- same convention solange_dmrg.py already
    # uses elsewhere in this project (--df-auxbasis default). Comparable size to the ETB
    # basis in local testing, so this is about using a vetted basis, not about the disk
    # issue above (which is the scratch-location fix, not this).
    mf.with_df.auxbasis = "def2-universal-jkfit"
    mf.kernel()
    if not mf.converged:
        print("*** SCF did NOT converge -- treat any active space below as provisional")
    ss, mult = mf.spin_square()
    print(f"[spin] high-spin reference: <S^2>={ss:.4f} (expect 30.0 for S=5), "
          f"2S+1={mult:.4f} (expect 11.0)")
    if abs(ss - 30.0) > 0.1:
        print("*** <S^2> does not match the expected S=5 value -- reference may have "
              "converged to a different spin state than intended; check before trusting "
              "any orbitals derived from it")
    return mf


def avas_at_threshold(mf, ao_labels, threshold):
    """AVAS on the high-spin reference. Returns (spec_dict, mo_coeff) --
    the orbitals themselves are what gets reused for the S=0 CASCI below,
    not just the (ncas, nelecas) counts."""
    from pyscf.mcscf import avas
    aos = [s.strip() for s in ao_labels.split(",")]
    ncas, nelecas, mo = avas.avas(mf, aos, threshold=threshold)
    print(f"[avas] threshold={threshold} AOs={aos} -> CAS({nelecas},{ncas}) "
          f"= {2 * ncas} qubits under Jordan-Wigner")
    spec = dict(threshold=threshold, ncas=int(ncas), nelecas=int(nelecas), qubits=int(2 * ncas))
    return spec, mo


def build_s0_casci(mf_highspin, mo_coeff, charge, ncas, nelecas):
    """S=0 (antiferromagnetic) state, CASCI on the HIGH-SPIN reference's own
    orbitals -- not a fresh RHF mean-field. This is the fix for the defect
    above: both spin states are read off the identical orbital set, so any
    difference between them is attributable to spin coupling alone, not to
    starting from two different (and for S=0, invalid) references.

    nelecas here must be expressed as (n_alpha, n_beta) for the S=0 CASCI
    call, not the single nelecas count AVAS reports for the high-spin case --
    PySCF's CASCI takes a tuple when spin != (nelec_total convention), and at
    S=0 alpha=beta=nelecas//2 exactly (even electron count expected for this
    cluster; fails loudly via the assertion below if that's not the case,
    rather than silently building the wrong state).
    """
    from pyscf import mcscf
    mol_s0 = mf_highspin.mol.copy()
    mol_s0.spin = 0
    mol_s0.charge = charge
    mol_s0.build()
    assert nelecas % 2 == 0, (
        f"nelecas={nelecas} is odd -- S=0 (closed-shell-electron-count) CASCI is not "
        f"well-defined on this active space; re-check the AVAS selection before proceeding")
    na = nb = nelecas // 2
    mc = mcscf.CASCI(mf_highspin, ncas, (na, nb))
    mc.mol = mol_s0
    mc.mo_coeff = mo_coeff
    mc.kernel()
    # mc.spin_square() itself raises AttributeError here: mcscf.CASCI() built on a
    # density-fitted mean-field (as build_mf() above always produces) returns a
    # DFCASCI instance, which does not inherit spin_square() the way plain CASCI does.
    # Found live 2026-10-07 testing this exact pattern on a stretched-H2 analog before
    # trusting it on the real cluster -- go through fcisolver.spin_square() directly,
    # which works on both variants.
    ss, mult = mc.fcisolver.spin_square(mc.ci, mc.ncas, mc.nelecas)
    print(f"[spin] S=0 CASCI: <S^2>={ss:.4f} (expect 0.0), 2S+1={mult:.4f} (expect 1.0)")
    if abs(ss) > 0.1:
        print("*** <S^2> does not match the expected S=0 value for this CASCI solution -- "
              "the solver may have found a different spin state within the same orbital "
              "space; check the CI vector / active-space symmetry before trusting the energy")
    return mc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--threshold", type=float, default=0.2,
                     help="AVAS threshold on the high-spin reference. Claude Science's point 2: "
                          "the Fe 3d + S 3p target AO pool (10+18=28 AOs) equals the cluster's own "
                          "28 orbitals exactly -- pool==cut, so there is nothing for a threshold "
                          "sweep to select between. Kept as a single value, not swept, and recorded "
                          "as such; do not read a sweep's worth of meaning into this one number.")
    ap.add_argument("--basis", default="def2-svp",
                     help="Changed from ccpvdz (this script's own prior default) per Claude "
                          "Science's 2026-10-07 point 4: Fe 3d needs better polarization than "
                          "Dunning-family bases provide at this size; def2-svp is the stated "
                          "minimum, def2-tzvp the preferred choice if affordable. This does NOT "
                          "need to match any other target's basis (e.g. TP53 C275F's 6-31g) -- R "
                          "is computed within SDHB against its own size-matched control, so the "
                          "control just needs to share THIS basis, not any other system's.")
    ap.add_argument("--fe-semicore", choices=["core", "active"], default="core",
                     help="Explicit decision point per Claude Science's point 3: whether Fe 3s/3p "
                          "semi-core orbitals are frozen (ncore) or included in the active space. "
                          "Default 'core' (standard practice; 3s/3p correlation effects are "
                          "normally small relative to 3d) -- recorded explicitly here rather than "
                          "left to whatever AVAS/PySCF defaults happen to do, per the point's own "
                          "instruction that this needs one line in the record, not silence.")
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()

    verify_wt()
    if a.verify_only:
        return

    ph_pdb = protonate(a.ph)
    # Found live 2026-10-07: "RuntimeError: Ill geometry" from duplicate atom
    # coordinates (4 pairs, exactly the FES cofactor's own atom count). The old
    # comment on read_fes_from_cif() assumed gemmi's write_pdb() always drops
    # FES HETATM records during cif->pdb conversion -- that assumption just
    # proved wrong (gemmi/pdbfixer version-dependent, apparently, since this
    # exact code path ran clean in an earlier round). Stop relying on the
    # conversion dropping them; explicitly exclude any FES read from the PDB
    # instead, so read_fes_from_cif()'s own direct cif read is the ONLY source
    # of those 4 atoms regardless of what the PDB conversion does or doesn't do.
    all_pdb_atoms = read_pdb(ph_pdb)
    pdb_atoms = [a for a in all_pdb_atoms if a["comp"] != "FES"]
    dropped = len(all_pdb_atoms) - len(pdb_atoms)
    if dropped:
        print(f"[cluster] NOTE: {dropped} FES atom(s) were present in the protonated PDB "
              f"this run (contradicting this script's prior assumption that gemmi drops "
              f"them) -- excluded here; the single source of truth for FES is "
              f"read_fes_from_cif()'s direct cif read, below.")
    atoms = pdb_atoms + read_fes_from_cif(PDB_ID, CHAIN, "301")
    groups, caps, charge = build_cluster(atoms)
    n_heavy = sum(1 for k in groups for a in groups[k] if a["elem"] != "H")
    print(f"[cluster] {len(groups)} residues, {n_heavy} heavy atoms, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec expects {EXPECTED_NET_CHARGE:+d})")
    if charge != EXPECTED_NET_CHARGE:
        print(f"[cluster] NOTE: charge differs from spec -- check whether pdbfixer assigned "
              f"a spurious HG to a metal-bound Cys thiolate (known failure mode; see "
              f"run_gate2_avas.py's build_cluster comments) before trusting this number.")
    xyz = "sdhb_c101_cluster.xyz"
    # Single cluster geometry for BOTH spin states -- same core, same charge, same basis,
    # same active space, per Claude Science's own framing of this as the built-in control
    # (R_spin = dE(S=0)/dE(S=5) at matched M -- the only difference between the two legs
    # is spin coupling, nothing else).
    write_xyz(groups, caps, xyz, comment=f"SDHB [2Fe-2S] S1 site, {PDB_ID}, Fe3d+S3p AVAS")
    print(f"[cluster] wrote {xyz}")
    print(f"[core] Fe 3s/3p treated as {a.fe_semicore} (--fe-semicore)")

    # Step 1: high-spin (S=5) reference + AVAS -- the only mean-field build in this script.
    mf_hs = build_mf(xyz, charge, spin=10, basis=a.basis)
    spec, mo = avas_at_threshold(mf_hs, "Fe 3d, S 3p", a.threshold)
    ncas, nelecas = spec["ncas"], spec["nelecas"]

    mo_path = f"{a.out_dir}/sdhb_s1_mo_coeff_highspin.npy"
    import numpy as np
    np.save(mo_path, mo)
    print(f"[avas] saved shared orbitals -> {mo_path} (reused for S=0 below, "
          f"and for any downstream solange_dmrg.py --load-orbitals run)")

    # Step 2: S=0 (antiferromagnetic) CASCI on those SAME orbitals -- no independent
    # mean-field, per the fix documented in build_mf()'s and build_s0_casci()'s docstrings.
    mc_s0 = build_s0_casci(mf_hs, mo, charge, ncas, nelecas)

    out = f"{a.out_dir}/sdhb_s1_both_spins.json"
    json.dump(dict(
        pdb=PDB_ID, charge=charge, basis=a.basis, fe_semicore=a.fe_semicore,
        threshold=a.threshold, ncas=ncas, nelecas=nelecas, qubits=spec["qubits"],
        high_spin=dict(spin=10, e_scf=float(mf_hs.e_tot), scf_converged=bool(mf_hs.converged),
                        spin_square_expected=30.0),
        casci_s0=dict(spin=0, e_tot=float(mc_s0.e_tot), spin_square_expected=0.0),
        mo_coeff_path=mo_path,
    ), open(out, "w"), indent=1)
    print(f"\nwrote {out}")
    print(f"\n[next] feed ncas={ncas}, nelecas={nelecas}, --load-orbitals {mo_path} into "
          f"solange_dmrg.py's --geometry path for the actual bond-dimension ladder (both spin "
          f"states), plus the matched size/filling control (adamantane C10H16 or n-C12H26 at "
          f"CAS({nelecas},{ncas}), same basis) before reading any R value.")


if __name__ == "__main__":
    main()
