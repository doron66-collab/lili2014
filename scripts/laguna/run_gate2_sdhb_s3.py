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
# Optional backbone-only extension (--extend-backbone), per Claude Science's
# 2026-10-10 reply: the 43-atom model drops 10 of the 12 N-H...S hydrogen bonds
# the protein donates to the cluster sulfurs (measured in 8GS8, N...S < 3.8 A;
# see s3_model_extension.csv), which makes the sulfurs far more electron-rich
# than in the protein and ligand oxidation competitive with three Fe(III).
# Backbone N/H/CA/HA/C/O plus CB (side chain capped with H at CB) restores
# all 12 donors without bringing in any titratable side chain (HIS244's ring
# stays out, so no protonation decision is needed). Neutral: charge unchanged.
BACKBONE_EXTENSION = {("B", s) for s in ("198", "199", "244", "245", "246", "247", "248", "250")}
BACKBONE_NAMES = {"N", "H", "CA", "HA", "C", "O", "CB"}
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


def build_cluster(atoms, extend_backbone=False):
    """Same approach as run_gate2_sdhb.py's build_cluster(): explicit charge
    assignment from the cofactor's own formal core charge + thiolates, not a
    generic metal-distance heuristic. extend_backbone adds BACKBONE_EXTENSION
    (backbone + CB only, neutral) -- see that constant's comment."""
    groups = {}
    for a in atoms:
        key = (a["ch"], a["seq"])
        if key in REGION_RESIDUES:
            groups.setdefault(key, []).append(a)
        elif extend_backbone and key in BACKBONE_EXTENSION and (
                a["name"] in BACKBONE_NAMES or a["name"].startswith("HB")):
            groups.setdefault(key, []).append(a)
    if extend_backbone:
        missing_bb = [k for k in BACKBONE_EXTENSION if k not in groups]
        if missing_bb:
            sys.exit(f"backbone-extension residues not found: {missing_bb}")

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

    # Side-chain caps for backbone-only residues: every heavy atom bonded to CB
    # (other than CA) that was left out becomes an H on the CB->X vector.
    n_sc_caps = 0
    if extend_backbone:
        for key in sorted(BACKBONE_EXTENSION):
            res_atoms = [a for a in atoms if (a["ch"], a["seq"]) == key]
            cb = [a for a in res_atoms if a["name"] == "CB"]
            if not cb:
                continue
            p = cb[0]["xyz"]
            for x in res_atoms:
                if x["elem"] == "H" or x["name"] in BACKBONE_NAMES:
                    continue
                v = [x["xyz"][i] - p[i] for i in range(3)]
                d = math.sqrt(sum(t * t for t in v))
                if d < 1.9:
                    caps.append(("H", tuple(p[i] + v[i] / d * CAP_BOND_LENGTH for i in range(3))))
                    n_sc_caps += 1
        print(f"[cluster] backbone extension: {len(BACKBONE_EXTENSION)} residues, "
              f"{n_sc_caps} side-chain H caps at CB")

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


def build_mf(xyz, charge, spin, basis="def2-tzvp", level_shift=0.3, max_cycle=80,
             chkfile="sdhb_s3_scf.chk", guess_basis="def2-svp", newton_max_cycle=60,
             xc=None, solvent_eps=None, damp=0.0):
    """High-spin (S=15/2, spin=15) ROHF reference -- the single-determinant,
    genuinely well-behaved state. Same reasoning as run_gate2_sdhb.py's
    build_mf(): never build the low-spin (S=1/2) state's own independent
    mean-field, since it's open-shell/multi-reference despite having a
    well-defined spin quantum number.

    Convergence (added 2026-10-08): the plain ROHF here failed to converge
    on two consecutive real runs, and the two runs did NOT agree with each
    other -- E_SCF differed by 78 mHa (-7308.872 vs -7308.794) and AVAS then
    picked different active spaces (CAS(53,38) vs CAS(53,39)) from the same
    input. An unconverged, non-reproducible reference makes every number
    downstream of it unreportable, regardless of how good the DMRG is.
    Three standard measures, in order of preference:
      1. chkfile -- if a previous run converged and saved its orbitals, start
         from those: reproducible by construction, and much faster.
      2. Otherwise, converge first in a smaller basis (guess_basis), project
         that density into the target basis, and use it as the starting
         guess -- a much better start for a 3-Fe open-shell cluster than
         PySCF's default atomic-superposition guess.
      3. level_shift on the virtual orbitals throughout -- damps the
         occupied/virtual oscillation that is the usual failure mode for
         transition-metal ROHF. A level shift changes the PATH, not the
         converged fixed point: verified on an O2-triplet toy before using it
         here (shifted and unshifted converged energies agree to 1e-13 Ha;
         reloading from the chkfile reproduces to 5e-13 Ha).
    """
    from pyscf import gto, scf, lib
    from pyscf.scf import addons
    import numpy as np
    if spin != 15:
        raise ValueError("build_mf() only builds the high-spin (S=15/2, spin=15) reference -- "
                          "see this function's docstring. The S=1/2 state is evaluated by "
                          "build_lowspin_casci() on these same orbitals.")
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    # Same disk-space fix as run_gate2_sdhb.py -- set BEFORE any DF object is
    # built (the small-basis guess below builds one too).
    dftmp_dir = os.path.join(os.getcwd(), "pyscf_tmp")
    os.makedirs(dftmp_dir, exist_ok=True)
    lib.param.TMPDIR = dftmp_dir
    print(f"[scratch] PySCF density-fitting scratch redirected to {dftmp_dir}")

    def _rohf(m):
        # xc: orbitals from restricted-open-shell DFT (ROKS) instead of ROHF.
        # solvent_eps: ddCOSMO dielectric screening. Both per Claude Science's
        # 2026-10-10 reply: the ROHF minimum is a ligand-radical state of a
        # bare dianion in vacuum; HF's missing dynamic correlation
        # systematically destabilizes compact high-spin d5 configurations, and
        # the functional only has to produce ORBITALS -- DMRG recomputes the
        # energy inside the active space (BP86 recommended: stable on Fe-S).
        if xc:
            from pyscf import dft
            f = dft.ROKS(m).density_fit()
            f.xc = xc
        else:
            f = scf.ROHF(m).density_fit()
        f.with_df.auxbasis = "def2-universal-jkfit"
        f.level_shift = level_shift
        f.max_cycle = max_cycle
        f.damp = damp
        if solvent_eps:
            f = f.ddCOSMO()
            f.with_solvent.eps = solvent_eps
        return f

    def _as_reference(f):
        """The object every downstream step uses (AVAS, CAS integrals, both gates).

        For a plain ROHF run that is the SCF itself. For ROKS and/or ddCOSMO it is
        a gas-phase, density-fitted ROHF object carrying the SAME orbitals, whose
        e_tot is the Hartree-Fock energy of that determinant -- the quantity the
        CAS Hamiltonian actually reproduces, and therefore the right reference for
        the embedding gate and the [embedding-hs] check (E_KS is not: it is a
        different functional of the same orbitals). The active-space Hamiltonian
        is gas phase; R is a within-system ratio on one shared orbital set, so a
        constant environment term cancels in it.
        """
        if not xc and not solvent_eps:
            return f
        ref = scf.ROHF(f.mol).density_fit()
        ref.with_df.auxbasis = "def2-universal-jkfit"
        ref.mo_coeff, ref.mo_occ, ref.mo_energy = f.mo_coeff, f.mo_occ, f.mo_energy
        ref.e_tot = ref.energy_tot(ref.make_rdm1())
        ref.converged = f.converged
        ref.e_scf_source = float(f.e_tot)
        print(f"[scf] orbitals from {'ROKS/' + xc if xc else 'ROHF'}"
              f"{f' + ddCOSMO(eps={solvent_eps})' if solvent_eps else ''}: E_SCF={f.e_tot:.8f}; "
              f"gas-phase HF energy of that determinant (gate reference) = {ref.e_tot:.8f}")
        return ref

    def _converge(f, dm0, label):
        """DIIS + level shift first (cheap per iteration); if that runs out of
        cycles, CONTINUE from exactly where it stopped with second-order SCF
        (PySCF's co-iterative augmented-Hessian `newton`) -- far more robust
        for multi-iron open-shell references, at a higher cost per step.
        Added 2026-10-08 after the def2-svp stage itself failed to converge
        in 300 level-shifted DIIS cycles on the real cluster: the problem is
        the SCF landscape, not the basis size, and more DIIS cycles alone
        would not have fixed it."""
        f.kernel(dm0)
        if f.converged:
            return f
        print(f"[scf] {label}: DIIS+level-shift did not converge in {f.max_cycle} cycles "
              f"(E={f.e_tot:.8f}) -- continuing from that point with second-order SCF (newton)")
        g = f.newton()
        g.max_cycle = newton_max_cycle
        g.kernel(f.mo_coeff, f.mo_occ)
        print(f"[scf] {label}: second-order SCF converged={g.converged} E={g.e_tot:.8f}")
        return g

    def _report_s2(f):
        ss, mult = f.spin_square()
        print(f"[spin] high-spin reference: <S^2>={ss:.4f} (expect 63.75 for S=15/2), "
              f"2S+1={mult:.4f} (expect 16.0)")
        if abs(ss - 63.75) > 0.5:
            print("*** <S^2> does not match the expected S=15/2 value -- reference may have "
                  "converged to a different spin state than intended; check before trusting "
                  "any orbitals derived from it")

    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin, verbose=4)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, charge={charge} spin={spin}")
    mf = _rohf(mol)

    dm0 = None
    if chkfile and os.path.exists(chkfile):
        # Restart goes STRAIGHT to second-order SCF from the saved orbitals --
        # no DIIS, no level shift. Found 2026-10-10 (job 2350814): the old path
        # rebuilt only a density matrix from the chkfile and re-ran level-
        # shifted DIIS; cycle 1 sat at the converged point (|g|=3.6e-6,
        # dE=-6.6e-11), then DIIS walked AWAY from it over 80 cycles to an
        # energy 0.71 Ha higher, and newton spent the rest of a 12-hour
        # allocation climbing back. Newton from converged orbitals finishes in
        # one or two macro-iterations, so this costs nothing when the saved
        # state is good and still converges when it isn't.
        try:
            from pyscf.scf import chkfile as _chk
            saved = _chk.load(chkfile, "scf")
            print(f"[scf] restarting from saved orbitals in {chkfile} with second-order SCF "
                  f"(saved E={float(saved['e_tot']):.8f})")
            g = _rohf(mol)
            g.level_shift = 0.0
            g.chkfile = chkfile
            g = g.newton()
            g.max_cycle = newton_max_cycle
            g.kernel(saved["mo_coeff"], saved["mo_occ"])
            print(f"[scf] restart: second-order SCF converged={g.converged} E={g.e_tot:.8f}")
            print(f"[scf] converged={g.converged} E_SCF={g.e_tot:.8f}")
            if g.converged:
                g = _as_reference(g)
                _report_s2(g)
                return g
            print("[scf] restart did not converge -- falling back to the full two-stage path")
        except Exception as e:
            print(f"[scf] could not restart from {chkfile} ({type(e).__name__}: {e}) -- ignoring it")
    if dm0 is None and guess_basis:
        print(f"[scf] stage 1: converging {'ROKS/' + xc if xc else 'ROHF'} in {guess_basis} as a starting guess "
              f"(level_shift={level_shift}, max_cycle={max_cycle})")
        mol_s = gto.M(atom="\n".join(lines), basis=guess_basis, charge=charge, spin=spin, verbose=4)
        mf_s = _converge(_rohf(mol_s), None, f"stage 1 ({guess_basis})")
        print(f"[scf] stage 1 ({guess_basis}): converged={mf_s.converged} E={mf_s.e_tot:.8f}")
        dm0 = np.array([addons.project_dm_nr2nr(mol_s, d, mol) for d in mf_s.make_rdm1()])

    if chkfile:
        mf.chkfile = chkfile
    print(f"[scf] stage 2: {'ROKS/' + xc if xc else 'ROHF'} in {basis} (level_shift={level_shift}, "
          f"damp={damp}, max_cycle={max_cycle})")
    mf = _converge(mf, dm0, f"stage 2 ({basis})")
    print(f"[scf] converged={mf.converged} E_SCF={mf.e_tot:.8f}")
    if not mf.converged:
        print("*** SCF did NOT converge -- treat any active space below as provisional")
    mf = _as_reference(mf)
    _report_s2(mf)
    return mf


def avas_at_threshold(mf, ao_labels, threshold, openshell_option=3):
    # openshell_option=3 keeps every singly-occupied ROHF orbital in the
    # active space. PySCF's default (2) projects SOMOs together with the
    # doubly-occupied orbitals, so "core" orbitals can carry SOMO character
    # and the ROHF determinant is no longer inside the CAS. Found 2026-10-09:
    # job 2350799's S=15/2 DMRG leg converged 3.6 Ha ABOVE E_ROHF on option-2
    # orbitals (diag_sdhb_s3_embedding.py reproduces the leak on a small
    # Fe-S model: beta core occupation 0.985 and +31.6 mHa with option 2,
    # exactly 1.000 and 0.000 mHa with option 3).
    from pyscf.mcscf import avas
    import copy
    import numpy as np
    # Option 3 slices the MO list by POSITION (doubly occupied = the first
    # nocc-spin columns, singly occupied = the next spin columns), not by
    # mo_occ. Second-order SCF / ROKS can return a SOMO energy-ordered among
    # the virtuals, and AVAS then silently mixes it into the virtual space.
    # Found 2026-10-10 on a small Fe-S model with ROKS orbitals: 0.008 e leaked
    # to virtuals and the embedding check missed by +19.9 mHa; sorting by
    # occupation first gave exactly 0.000 mHa.
    order = np.argsort(-np.asarray(mf.mo_occ), kind="stable")
    mf = copy.copy(mf)
    mf.mo_coeff, mf.mo_occ = mf.mo_coeff[:, order], mf.mo_occ[order]
    mf.mo_energy = mf.mo_energy[order]
    aos = [s.strip() for s in ao_labels.split(",")]
    ncas, nelecas, mo = avas.avas(mf, aos, threshold=threshold,
                                  openshell_option=openshell_option)
    print(f"[avas] openshell_option={openshell_option} threshold={threshold} AOs={aos} -> CAS({nelecas},{ncas}) "
          f"= {2 * ncas} qubits under Jordan-Wigner")
    spec = dict(threshold=threshold, ncas=int(ncas), nelecas=int(nelecas), qubits=int(2 * ncas))
    return spec, mo


def build_lowspin_casci(mf_highspin, mo_coeff, charge, ncas, nelecas,
                         bond_dims=(250, 500, 1000, 2000), scratch="./tmp_dmrg_s3",
                         n_threads=4, max_minutes=None):
    """S=1/2 (ground, antiferromagnetically coupled) state, on the HIGH-SPIN
    reference's own orbitals. Generalized from run_gate2_sdhb.py's
    build_s0_casci() for an ODD na-nb difference (=1, not =0): S3 has 15
    d-electrons (odd), so there is no S=0 state here at all -- the true
    ground state carries a half-integer spin quantum number by construction,
    not by choice.

    CORRECTED 2026-10-07: this used to call mc.kernel() -- PySCF's default
    EXACT FCI solver -- on CAS(53,38). That is not a resource/memory tuning
    problem, it is combinatorially impossible (confirmed live: block2's own
    link-string index alone needed 5.67 TiB). Exact FCI is only tractable to
    roughly CAS(16,16); every active space in this project past that size
    (including this script's own module docstring, which already said this
    CAS would land around 27 orbitals -- itself already past that limit)
    goes through DMRG (block2), same as C275F and everything else. This was
    an unflagged planning gap, not a deliberate placeholder -- the module
    docstring's own numbers should have caught it before the script was
    first run.

    Does NOT call mc.kernel(). Pulls the active-space integrals (h1e, h2e,
    ecore) directly off the un-run CASCI object via get_h1eff()/get_h2eff()
    -- the exact same accessor pattern solange_dmrg.py's own CASSCF path
    already uses to hand integrals to run_dmrg() (see that function's own
    comment on why mo_coeff=None there is safe) -- then calls
    solange_dmrg.run_dmrg() directly with spin=1.
    """
    from pyscf import mcscf, ao2mo
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
    target_spin = na - nb
    # The guard Doron asked for explicitly: this function is ONLY valid for
    # the odd-na-nb-difference-of-1 case it was built for. A future caller
    # (copy-paste onto a different cluster) passing an nelecas/spin pair this
    # function didn't derive itself must fail loudly here, not run DMRG in
    # the wrong symmetry sector and report a confident, wrong energy.
    assert target_spin == 1, (
        f"derived spin (na-nb)={target_spin} from nelecas={nelecas}, but this function "
        f"only implements the S=1/2 (spin=1) case it was built for -- do not reuse it "
        f"for a different target without re-deriving na/nb and this assertion together")
    mc = mcscf.CASCI(mf_highspin, ncas, (na, nb))
    mc.mol = mol_ls
    mc.mo_coeff = mo_coeff
    h1e, ecore = mc.get_h1eff(mo_coeff=mo_coeff)
    h2e = ao2mo.restore(1, mc.get_h2eff(mo_coeff), ncas)
    print(f"[dmrg] CAS({nelecas},{ncas}) is far beyond exact FCI's ~CAS(16,16) practical "
          f"ceiling -- solving via DMRG (block2), spin={target_spin}, bond_dims={list(bond_dims)}")

    import sys as _sys
    _sys.path.insert(0, os.getcwd())
    import solange_dmrg
    energies, s_max, stop_reason, discarded_weights, sweep_history = solange_dmrg.run_dmrg(
        h1e, h2e, ecore, ncas, nelecas, list(bond_dims), scratch=scratch,
        n_threads=n_threads, max_minutes=max_minutes, spin=target_spin)

    # <S^2> is NOT independently measured here -- same reasoning as
    # dmrgscf_block2.py's own Block2FCISolver.spin_square(): SU2 symmetry
    # mode is spin-adapted by construction, so a converged state in the
    # spin=1 sector IS an S=1/2 eigenstate by construction, not by
    # verification. That existing analytical-reporting pattern is reused
    # here rather than guessing at an unverified pyblock2 API call for a
    # real spin-expectation measurement (none was found in this codebase to
    # copy; DP1 -- do not invent one for an HPC-only library this sandbox
    # cannot run live).
    s_expected = target_spin / 2.0
    print(f"[spin] S=1/2 DMRG: <S^2>={s_expected * (s_expected + 1):.4f} (expect 0.75, "
          f"analytically from the SU2 spin=1 sector requested -- not independently measured)")

    return dict(energies=energies, s_max=s_max, stop_reason=stop_reason,
                discarded_weights=discarded_weights, sweep_history=sweep_history,
                e_tot=(energies[-1][1] if energies else None), ncas=ncas, nelecas=nelecas,
                spin=target_spin)


def build_highspin_dmrg(mf_highspin, mo_coeff, charge, ncas, nelecas,
                        bond_dims=(250, 500, 1000, 2000), scratch="./tmp_dmrg_s3_hs",
                        n_threads=4, max_minutes=None):
    """S=15/2 (fully ferromagnetic) leg, DMRG on the SAME shared orbitals and
    the SAME active space as build_lowspin_casci() -- the second spin leg the
    protocol needs for R_spin, and (added 2026-10-08) the cleanest available
    check on the integral bookkeeping itself.

    Why it is a check: the high-spin ROHF determinant lies INSIDE this
    spin=15 active-space sector (AVAS's core/active split keeps the 15
    singly-occupied orbitals active -- confirmed by the even core electron
    count, 309-53=256), so the exact active-space energy here can only be at
    or BELOW E_ROHF. A converged DMRG energy that lands well ABOVE E_ROHF
    means h1e/h2e/ecore are wrong, not that the chemistry is interesting.
    That is exactly the open question raised by the first real S=1/2 run
    (M=250 landed 1.6 Ha above E_ROHF): if this leg also lands ~1.6 Ha high,
    the bug is in the integrals; if it lands at/below E_ROHF, the S=1/2 gap
    was the cold-start artifact already documented for R175H.

    Separate scratch directory on purpose: run_dmrg()'s resume check keys on
    (ncas, nelecas[, spin]) and both legs share ncas/nelecas.
    """
    from pyscf import mcscf, ao2mo
    target_spin = 15
    assert nelecas % 2 == 1 and nelecas >= target_spin, (
        f"nelecas={nelecas} cannot host a spin={target_spin} state")
    na = (nelecas + target_spin) // 2
    nb = (nelecas - target_spin) // 2
    assert na - nb == target_spin
    mc = mcscf.CASCI(mf_highspin, ncas, (na, nb))
    mc.mo_coeff = mo_coeff
    h1e, ecore = mc.get_h1eff(mo_coeff=mo_coeff)
    h2e = ao2mo.restore(1, mc.get_h2eff(mo_coeff), ncas)
    print(f"[dmrg-hs] CAS({nelecas},{ncas}) high-spin leg, spin={target_spin}, "
          f"bond_dims={list(bond_dims)}, scratch={scratch}")
    import sys as _sys
    _sys.path.insert(0, os.getcwd())
    import solange_dmrg
    energies, s_max, stop_reason, discarded_weights, sweep_history = solange_dmrg.run_dmrg(
        h1e, h2e, ecore, ncas, nelecas, list(bond_dims), scratch=scratch,
        n_threads=n_threads, max_minutes=max_minutes, spin=target_spin)
    e_last = energies[-1][1] if energies else None
    e_scf = float(mf_highspin.e_tot)
    if e_last is not None:
        gap_mha = (e_last - e_scf) * 1000.0
        print(f"[embedding-hs] E_DMRG(S=15/2, M={energies[-1][0]})={e_last:.8f}  "
              f"E_ROHF={e_scf:.8f}  difference={gap_mha:+.3f} mHa "
              f"(must be <= ~0 if the integrals are right)")
        if gap_mha > 1.6:
            print("*** E_DMRG(high-spin) is ABOVE the ROHF energy it contains -- the active-"
                  "space integrals (h1e/h2e/ecore) or the core/active split are suspect. Do NOT "
                  "interpret the S=1/2 leg until this is explained. ***")
    return dict(energies=energies, s_max=s_max, stop_reason=stop_reason,
                discarded_weights=discarded_weights, sweep_history=sweep_history,
                e_tot=e_last, ncas=ncas, nelecas=nelecas, spin=target_spin)


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
    ap.add_argument("--bond-dims", default="250,500,1000,2000",
                     help="DMRG bond-dimension ladder for the S=1/2 low-spin solve (CAS this "
                          "large has no exact-FCI answer to fall back on -- see "
                          "build_lowspin_casci()'s docstring). Same default ladder as "
                          "solange_dmrg.py's own --bond-dims.")
    ap.add_argument("--dmrg-scratch", default="./tmp_dmrg_s3",
                     help="Reused across resubmissions to resume the MPS instead of "
                          "restarting the bond-dimension ladder from scratch -- same "
                          "mechanism as solange_dmrg.py's --scratch.")
    ap.add_argument("--scf-level-shift", type=float, default=0.3,
                     help="Virtual-orbital level shift (Ha) for the high-spin ROHF. Changes the "
                          "convergence path, not the converged answer.")
    ap.add_argument("--scf-max-cycle", type=int, default=80,
                     help="DIIS cycles before handing over to second-order SCF (newton). Was 300: on "
                          "the real cluster DIIS+level shift failed to converge in 300 cycles even in "
                          "def2-svp, and at ~2.3 min/cycle in def2-tzvp that is ~11.5 h spent before "
                          "the far more robust newton step even starts.")
    ap.add_argument("--xc", default=None,
                     help="Take orbitals from restricted-open-shell DFT with this functional "
                          "(e.g. BP86) instead of ROHF. Claude Science 2026-10-10.")
    ap.add_argument("--solvent-eps", type=float, default=None,
                     help="ddCOSMO dielectric screening during the SCF (protein interior ~4).")
    ap.add_argument("--extend-backbone", action="store_true",
                     help="Add backbone + CB of the 8 residues donating the missing N-H...S "
                          "hydrogen bonds (see BACKBONE_EXTENSION).")
    ap.add_argument("--mixed-basis", action="store_true",
                     help="def2-tzvp on Fe and S only, def2-svp on every other atom.")
    ap.add_argument("--scf-damp", type=float, default=0.0,
                     help="Density damping factor for the DIIS stages (0-1). For a large model "
                          "whose early DIIS oscillates by hundreds of Ha (job 2350824).")
    ap.add_argument("--scf-only", action="store_true",
                     help="Stop after the SCF, the state gate, AVAS and the embedding gate -- "
                          "no DMRG. For checking which electronic state a model setup lands in.")
    ap.add_argument("--scf-chkfile", default=None,
                     help="Saved SCF orbitals. Reused as the starting point on the next run if "
                          "present -- makes the reference reproducible across runs, and lets a "
                          "run that ran out of cycles continue instead of restarting. Default "
                          "is derived from the model/method flags, so a different setup never "
                          "restarts from another setup's orbitals.")
    ap.add_argument("--allow-unconverged-scf", action="store_true",
                     help="Proceed to AVAS/DMRG even if the high-spin SCF did not converge. Off by "
                          "default: two real runs with an unconverged reference gave different "
                          "active spaces from the same input, so nothing downstream is reportable. "
                          "Use only for a deliberate pipeline smoke test.")
    ap.add_argument("--skip-highspin-dmrg", action="store_true",
                     help="Skip the S=15/2 DMRG leg. Not recommended: it is both the second spin "
                          "leg R_spin needs and the integral-bookkeeping check.")
    ap.add_argument("--dmrg-threads", type=int, default=4)
    ap.add_argument("--dmrg-max-minutes", type=float, default=None,
                     help="Wall-clock budget for the DMRG ladder; stops requesting larger "
                          "bond dims once exceeded rather than being killed mid-sweep with "
                          "nothing recorded. Omit for no limit.")
    a = ap.parse_args()

    # One tag per model/method setup, used in every file name this run reads or
    # writes -- the old single 'sdhb_s3_scf.chk' now holds the ligand-radical
    # ROHF state, and no other setup may restart from it (or resume another
    # setup's DMRG MPS).
    tag = ("_bb" if a.extend_backbone else "") + ("_mix" if a.mixed_basis else "") + \
          (f"_{a.xc.lower()}" if a.xc else "") + (f"_eps{a.solvent_eps:g}" if a.solvent_eps else "")
    if a.scf_chkfile is None:
        a.scf_chkfile = f"sdhb_s3_scf{tag}.chk"
    if tag and a.dmrg_scratch == "./tmp_dmrg_s3":
        a.dmrg_scratch = f"./tmp_dmrg_s3{tag}"
    basis = ({"Fe": a.basis, "S": a.basis, "default": "def2-svp"} if a.mixed_basis else a.basis)
    print(f"[setup] tag='{tag or '(plain)'}' chkfile={a.scf_chkfile} dmrg_scratch={a.dmrg_scratch} "
          f"basis={basis} xc={a.xc or 'none (ROHF)'} solvent_eps={a.solvent_eps}")

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
    groups, caps, charge = build_cluster(atoms, extend_backbone=a.extend_backbone)
    n_heavy = sum(1 for k in groups for a in groups[k] if a["elem"] != "H")
    print(f"[cluster] {len(groups)} residues, {n_heavy} heavy atoms, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec expects {EXPECTED_NET_CHARGE:+d})")
    if charge != EXPECTED_NET_CHARGE:
        print(f"[cluster] NOTE: charge differs from spec -- check whether pdbfixer assigned "
              f"a spurious HG to a metal-bound Cys thiolate before trusting this number.")
    xyz = f"sdhb_s3_cluster{'_bb' if a.extend_backbone else ''}.xyz"
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

    mf_hs = build_mf(xyz, charge, spin=15, basis=basis, level_shift=a.scf_level_shift,
                     max_cycle=a.scf_max_cycle, chkfile=a.scf_chkfile,
                     xc=a.xc, solvent_eps=a.solvent_eps, damp=a.scf_damp)
    if not mf_hs.converged and not a.allow_unconverged_scf:
        sys.exit(f"\n*** REFUSING to continue: high-spin SCF did not converge "
                 f"(E={mf_hs.e_tot:.8f}). Its orbitals were saved to {a.scf_chkfile} -- "
                 f"re-running this same command resumes from them. Pass "
                 f"--allow-unconverged-scf only for a deliberate pipeline smoke test; "
                 f"nothing from such a run is reportable. ***")
    # State gate: the SCF has several S=15/2 solutions, and the lowest one
    # found so far is NOT the three-high-spin-Fe(III) state the active space is
    # designed for (see diag_sdhb_s3_spinpop.state_check). <S^2> can't tell
    # them apart, so read where the spin sits and refuse before AVAS/DMRG.
    from diag_sdhb_s3_spinpop import spin_populations, state_check
    _pop = spin_populations(mf_hs.mol, mf_hs.mo_coeff, mf_hs.mo_occ)
    for _ia in range(mf_hs.mol.natm):
        if mf_hs.mol.atom_symbol(_ia) == "Fe":
            print(f"[state] Fe atom {_ia} spin population {_pop[_ia]:+.3f}")
    _ok, _why = state_check(mf_hs.mol, _pop)
    if not _ok:
        sys.exit("\n*** REFUSING: the converged SCF is not the intended three-high-spin-Fe(III) "
                 "state (" + "; ".join(_why) + "). An active space built on it would describe a "
                 "different electronic state. ***")
    spec, mo = avas_at_threshold(mf_hs, "Fe 3d, S 3p", a.threshold)
    ncas, nelecas = spec["ncas"], spec["nelecas"]

    # Embedding gate BEFORE any DMRG: the ROHF determinant, evaluated with the
    # same h1e/h2e/ecore the DMRG legs will use, must reproduce E_ROHF. If it
    # doesn't, the core/active split excludes the reference and every DMRG
    # energy downstream is meaningless -- refuse in minutes, not after 7 hours.
    from diag_sdhb_s3_embedding import check as _embedding_check
    e_det = _embedding_check(mf_hs, mo, ncas, nelecas, "pre-DMRG embedding gate")
    if abs(e_det - mf_hs.e_tot) * 1000.0 > 1.0:
        sys.exit(f"\n*** REFUSING: ROHF determinant in the CAS integrals gives {e_det:.8f}, "
                 f"E_ROHF={mf_hs.e_tot:.8f} ({(e_det - mf_hs.e_tot) * 1000:+.3f} mHa). The "
                 f"active space does not contain the reference; no DMRG run on it is "
                 f"interpretable. ***")

    mo_path = f"{a.out_dir}/sdhb_s3_mo_coeff_highspin{tag}.npy"
    import numpy as np
    np.save(mo_path, mo)
    print(f"[avas] saved shared orbitals -> {mo_path}")
    if a.scf_only:
        print(f"\n[scf-only] state gate and embedding gate PASSED for setup '{tag or '(plain)'}': "
              f"CAS({nelecas},{ncas}) = {2 * ncas} qubits, ncore={(mf_hs.mol.nelectron - nelecas) // 2}. "
              f"Stopping before DMRG (--scf-only).")
        return

    bond_dims = [int(x) for x in a.bond_dims.split(",")]

    def _leg_json(r):
        return dict(spin=r["spin"], bond_dims=[m for m, _ in r["energies"]],
                    energies=[e for _, e in r["energies"]], e_tot=r["e_tot"],
                    s_max=r["s_max"], stop_reason=r["stop_reason"],
                    discarded_weights=r["discarded_weights"])

    # High-spin leg FIRST: it's the integral-bookkeeping check (see
    # build_highspin_dmrg's docstring) -- if it fails, the S=1/2 leg's hours
    # would be spent on integrals already shown to be wrong.
    hs_result = None
    if not a.skip_highspin_dmrg:
        hs_result = build_highspin_dmrg(
            mf_hs, mo, charge, ncas, nelecas, bond_dims=bond_dims,
            scratch=a.dmrg_scratch + "_hs", n_threads=a.dmrg_threads,
            max_minutes=a.dmrg_max_minutes)

    dmrg_result = build_lowspin_casci(
        mf_hs, mo, charge, ncas, nelecas, bond_dims=bond_dims,
        scratch=a.dmrg_scratch, n_threads=a.dmrg_threads, max_minutes=a.dmrg_max_minutes)

    out = f"{a.out_dir}/sdhb_s3_both_spins{tag}.json"
    json.dump(dict(
        pdb=PDB_ID, charge=charge, basis=str(basis), fe_semicore=a.fe_semicore,
        setup=dict(tag=tag, xc=a.xc, solvent_eps=a.solvent_eps,
                   extend_backbone=a.extend_backbone, mixed_basis=a.mixed_basis,
                   e_scf_source=getattr(mf_hs, "e_scf_source", None)),
        threshold=a.threshold, ncas=ncas, nelecas=nelecas, qubits=spec["qubits"],
        high_spin=dict(spin=15, e_scf=float(mf_hs.e_tot), scf_converged=bool(mf_hs.converged),
                        spin_square_expected=63.75),
        dmrg_highspin=(_leg_json(hs_result) if hs_result else None),
        dmrg_lowspin=dict(
            spin=1, spin_square_expected=0.75,
            bond_dims=[m for m, _ in dmrg_result["energies"]],
            energies=[e for _, e in dmrg_result["energies"]],
            e_tot=dmrg_result["e_tot"], s_max=dmrg_result["s_max"],
            stop_reason=dmrg_result["stop_reason"],
            discarded_weights=dmrg_result["discarded_weights"],
        ),
        mo_coeff_path=mo_path,
    ), open(out, "w"), indent=1)
    print(f"\nwrote {out}")
    if not mf_hs.converged:
        print("\n*** REMINDER: the high-spin SCF reference did NOT converge -- every number "
              "above (orbitals, DMRG energies, S_max) is provisional until that's fixed or "
              "re-checked. Do not report this as a final result as-is.")
    print(f"\n[next] P3 FIRST (not this target): validate DMRG/SHCI against exact FCI on "
          f"S1's CAS(10,10) -- cheap, must agree to ~1e-9 before anything here is readable. "
          f"The S=1/2 DMRG ladder for THIS cluster (CAS({nelecas},{ncas})) now runs inline, "
          f"above -- it no longer needs a separate solange_dmrg.py --geometry invocation. "
          f"Still needed before reading any R_spin value: (1) the high-spin reference "
          f"converging cleanly (see the REMINDER above if it didn't this run), and (2) a "
          f"size/filling-matched chemistry-free negative control at CAS({nelecas},{ncas}), "
          f"same basis -- not yet built by this script. Pre-registered predictions P1/P2 "
          f"(R_spin > 1 for S3 and S2; R ordering S2 > S3 > S1 by iron count) are falsifiable -- "
          f"report the result either way.")


if __name__ == "__main__":
    main()
