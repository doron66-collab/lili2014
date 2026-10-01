#!/usr/bin/env python3
"""Gate 2: verified-structure -> protonated QM cluster -> AVAS active space.

Run this on a machine with normal network access (e.g. under Claude Code).
It performs the step that could not be run in the sandbox: the actual AVAS
selection that yields (n_elec, n_orb).

The verification gate is enforced in code: AVAS will NOT run unless the stated
mutation is confirmed present in the deposited coordinates. That is the whole
point -- the retracted C275F result passed every record-consistency check and
failed only this one.

    pip install pyscf pdbfixer openmm requests
    # pdbfixer via conda is often easier:
    #   conda install -c conda-forge pdbfixer openmm pyscf

Usage
-----
    # list what is runnable
    python run_gate2_avas.py --list

    # the one target that is ready as cited
    python run_gate2_avas.py --target TP53_Y220C

    # a target needing coordination closure
    python run_gate2_avas.py --target TP53_G245S --radius 8.4

    # skip AVAS, just build and inspect the cluster
    python run_gate2_avas.py --target TP53_Y220C --cluster-only

Inputs expected in the working directory:
    active_space_spec.json   (charge / spin / AVAS AO set per target)
    mmcif_verify.py          (the verification gate)
Both are produced by the Gate 2 chemistry pass and are saved as artifacts.
"""
import argparse
import json
import math
import os
import sys
import urllib.request

CAP_BOND_LENGTH = 1.09          # matches build_qm_cluster.py
METALS = {"ZN", "MG", "MN", "FE", "NI", "CU", "CA", "CO"}
SIDECHAIN_Q = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}

# UniProt accession + wild-type/position/mutant, keyed as in targets.json
TARGET_ID = {
    "TP53_C275F":     ("P04637", "C", 275, "F"),
    "TP53_Y220C":     ("P04637", "Y", 220, "C"),
    "TP53_R175H":     ("P04637", "R", 175, "H"),
    "TP53_G245S":     ("P04637", "G", 245, "S"),
    "TP53_R249S":     ("P04637", "R", 249, "S"),
    "TP53_R282W":     ("P04637", "R", 282, "W"),
    "KEAP1_G333C":    ("Q14145", "G", 333, "C"),
    "SETD2_R1625C":   ("Q9BYW2", "R", 1625, "C"),
    "SMARCA4_R1192C": ("P51532", "R", 1192, "C"),
    "ARID1A_R1020S":  ("O14497", "R", 1020, "S"),
    "STK11_D194N":    ("Q15831", "D", 194, "N"),
}


def fetch(pdb_id, cache="struct"):
    os.makedirs(cache, exist_ok=True)
    cif = os.path.join(cache, f"{pdb_id}.cif")
    if not os.path.exists(cif):
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id}.cif", cif)
    sif = os.path.join("sifts", f"{pdb_id}.json")
    os.makedirs("sifts", exist_ok=True)
    if not os.path.exists(sif):
        with open(sif, "wb") as fh:
            fh.write(urllib.request.urlopen(
                f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"
                f"{pdb_id.lower()}", timeout=60).read())
    return cif, json.load(open(sif))


# ---------------------------------------------------------------- the gate
def gate(pdb_id, accession, wt, pos, mut):
    """Refuse to proceed unless the mutation is in the coordinates."""
    sys.path.insert(0, os.getcwd())
    import mmcif_verify as mv
    cif, sifts_json = fetch(pdb_id)
    r = mv.verify_sifts(cif, sifts_json, pdb_id, accession, wt, pos, mut)
    print(f"[gate] {pdb_id} {wt}{pos}{mut}: {r['verdict']}")
    print(f"[gate] {r['detail']}")
    if not r["verdict"].startswith("MUTANT_CONFIRMED"):
        print(f"\n*** REFUSING to build an active space. The mutation is not "
              f"present in {pdb_id}. Verdict: {r['verdict']}.\n"
              f"*** Model the mutation first, then re-run against the modelled "
              f"structure -- and record that the geometry is a model.")
        return None
    return r


# ------------------------------------------------------- protonation + cluster
def protonate(pdb_id, ph=7.4, out=None):
    """Protonate the WHOLE structure before extracting, not after.

    Protonating a fragment is unreliable; protonating first and then cutting
    carries correct hydrogens into the cluster.
    """
    out = out or f"struct/{pdb_id}_H.pdb"
    if os.path.exists(out):
        return out
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        sys.exit("pdbfixer/openmm required for protonation:\n"
                 "  conda install -c conda-forge pdbfixer openmm")
    # PATCHED: no network egress in this session (files.rcsb.org blocked by
    # org policy). Use the locally-verified structure instead of fetching --
    # struct/{pdb_id}.pdb was converted from the checksum-verified mmCIF
    # (struct/{pdb_id}.cif, matching manifest.json's sha256_cif) via gemmi,
    # not re-downloaded, so it is the identical verified coordinate set.
    local_pdb = f"struct/{pdb_id}.pdb"
    if not os.path.exists(local_pdb):
        sys.exit(f"no local structure at {local_pdb} -- network fetch is "
                  f"blocked in this environment, place a verified PDB/CIF "
                  f"conversion there first")
    fixer = PDBFixer(filename=local_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not build unresolved loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)   # keeps metals, drops waters
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out, "w") as fh:
        # keepIds=True: PDBFixer/OpenMM renumbers chains+residues sequentially
        # from 1 by default on write, which silently breaks every downstream
        # lookup by original author residue number (region_residues, the
        # mutation site itself). Preserve the original numbering instead.
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print(f"[prep] protonated at pH {ph} -> {out}")
    return out


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


def build_cluster(atoms, chain, resnum, radius):
    """Radius selection + coordination closure. Returns (residues, caps, charge)."""
    centre = [a for a in atoms
              if a["ch"] == chain and a["seq"] == str(resnum) and a["name"] == "CA"]
    if not centre:
        sys.exit(f"no CA for {chain}{resnum}")
    c = centre[0]["xyz"]
    groups = {}
    for a in atoms:
        if a["comp"] in ("HOH", "WAT"):
            continue
        groups.setdefault((a["ch"], a["seq"], a["comp"]), []).append(a)

    inside = {k for k, v in groups.items()
              if any(d(x["xyz"], c) <= radius for x in v if x["elem"] != "H")}

    # rule R2: coordination closure
    thiolates = set()
    for k, v in groups.items():
        if k[2] not in METALS:
            continue
        m = v[0]
        shell = {(a["ch"], a["seq"]) for a in atoms
                 if a["elem"] in ("N", "O", "S") and d(a["xyz"], m["xyz"]) <= 2.8}
        if k in inside or any((x[0], x[1]) in shell for x in inside):
            inside.add(k)
            for kk in list(groups):
                if (kk[0], kk[1]) in shell:
                    inside.add(kk)
                    if kk[2] == "CYS":
                        thiolates.add(kk)
            print(f"[R2] coordination closure pulled in {k[2]}{k[1]} "
                  f"and its shell {sorted(f'{x[1]}' for x in shell)}")

    # net charge from the protonation actually present
    charge = 0
    for k in inside:
        if k[2] in METALS:
            charge += 2
        elif k[2] == "CYS":
            has_hg = any(a["name"] in ("HG", "HG1") for a in groups[k])
            charge += 0 if has_hg else -1
        else:
            charge += SIDECHAIN_Q.get(k[2], 0)

    # rule R4: cap cut peptide bonds with H
    sel_res = {(k[0], int(k[1])) for k in inside if k[1].lstrip("-").isdigit()}
    caps = []
    for ch, n in sorted(sel_res):
        for nb, anchor, partner in ((n - 1, "N", "C"), (n + 1, "C", "N")):
            if (ch, nb) in sel_res:
                continue
            aa = [a for a in atoms if a["ch"] == ch and a["seq"] == str(n)
                  and a["name"] == anchor]
            bb = [a for a in atoms if a["ch"] == ch and a["seq"] == str(nb)
                  and a["name"] == partner]
            if not aa or not bb:
                continue
            p, q = aa[0]["xyz"], bb[0]["xyz"]
            v = [q[i] - p[i] for i in range(3)]
            nrm = math.sqrt(sum(x * x for x in v)) or 1.0
            caps.append(("H", tuple(p[i] + v[i] / nrm * CAP_BOND_LENGTH
                                    for i in range(3))))
    return inside, groups, caps, charge, thiolates


def write_xyz(inside, groups, caps, path, comment=""):
    rows = [(a["elem"], a["xyz"]) for k in sorted(inside) for a in groups[k]]
    rows += caps
    with open(path, "w") as fh:
        fh.write(f"{len(rows)}\n{comment}\n")
        for el, (x, y, z) in rows:
            fh.write(f"{el:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    return rows


# ---------------------------------------------------------------------- AVAS
# PATCHED: the SCF (the expensive part -- ~45 min wall clock at 883 basis
# functions for TP53_Y220C) does not depend on the AVAS threshold at all.
# The original run_avas() reran SCF from scratch for every threshold in a
# --sweep, so an N-point sweep cost N full SCF runs for no reason. Split so
# a sweep runs SCF exactly once and reuses the converged mf object.
#
# chkfile (2026-10-01, Claude Science's sign-off on CAS(36,35)/(36,34)): every
# caller of build_mf() -- avas_mp2_select.py, localize_active_space.py, and
# solange_dmrg.py's --load-orbitals path -- independently reruns the full
# density-fitted SCF on the same ~1000-function cluster from scratch (the
# 30-45 min cost on Laguna), even when nothing about the molecule changed and
# only the ORBITALS differ downstream (a different active-space cut, a
# localized vs canonical rotation). That wall-clock is exactly what Science's
# canonical-vs-localized three-point probe plan assumes is paid ONCE. Passing
# --chkfile lets pyscf's own mechanism (mf.chkfile + init_guess='chkfile')
# seed SCF from the previously converged density on a repeat call: this still
# runs and checks mf.converged for real (never trusts a stored result without
# re-verifying it), it just converges in ~1-2 cycles instead of from scratch
# when the geometry/charge/spin/basis are unchanged -- the only case this is
# used for.
def build_mf(xyz, charge, spin, basis="6-31g", chkfile=None):
    import os
    from pyscf import gto, scf
    lines = [l.strip() for l in open(xyz).read().splitlines()[2:] if l.strip()]
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin,
                verbose=3)
    print(f"[avas] {mol.natm} atoms, {mol.nao} basis functions, "
          f"charge={charge} spin={spin}")
    mf = (scf.RHF(mol) if spin == 0 else scf.ROHF(mol)).density_fit()
    if chkfile:
        mf.chkfile = chkfile
        if os.path.exists(chkfile):
            mf.init_guess = 'chkfile'
            print(f"[avas] seeding SCF from existing chkfile {chkfile} "
                  f"(still converges and verifies mf.converged -- not trusted blindly)")
    mf.kernel()
    if not mf.converged:
        print("*** SCF did NOT converge -- any active space below is provisional")
    return mf


def avas_at_threshold(mf, ao_labels, threshold):
    from pyscf.mcscf import avas
    aos = [s.strip() for s in ao_labels.split(",")]
    ncas, nelecas, mo = avas.avas(mf, aos, threshold=threshold)
    print(f"\n[avas] threshold={threshold}  AOs={aos}")
    print(f"[avas] CAS({nelecas},{ncas})  ->  {2 * ncas} qubits under "
          f"Jordan-Wigner")
    mol = mf.mol
    return dict(threshold=threshold, ncas=int(ncas), nelecas=int(nelecas),
                qubits=int(2 * ncas), nao=int(mol.nao),
                scf_converged=bool(mf.converged), e_scf=float(mf.e_tot))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target")
    ap.add_argument("--spec", default="active_space_spec.json")
    ap.add_argument("--radius", type=float, default=6.0)
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--threshold", type=float, default=0.2)
    ap.add_argument("--sweep", help="comma-separated AVAS thresholds, "
                                    "e.g. 0.1,0.15,0.2,0.3")
    ap.add_argument("--cluster-only", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()

    spec = json.load(open(a.spec))["targets"]
    if a.list or not a.target:
        print(f"{'target':<17}{'status':<30}{'charge':>7} {'spin':>5}  AVAS AOs")
        for k, v in spec.items():
            if "net_charge" not in v:
                continue
            print(f"{k:<17}{v['status']:<30}{v['net_charge']:>+7} "
                  f"{v['spin_2S']:>5}  {v['avas_ao_set']}")
        print("\nOnly READY / READY_AFTER_STRUCTURE_SWAP targets will pass the gate.")
        return

    t = spec[a.target]
    if "net_charge" not in t:
        sys.exit(f"{a.target} has no chemistry spec: {t.get('blocked_reason')}")
    pdb_id = t["structure_used_for_spec"].split()[0]
    chain = t["structure_used_for_spec"].split()[2]
    resnum = int(t["structure_used_for_spec"].split()[-1])
    acc, wt, pos, mut = TARGET_ID[a.target]

    if gate(pdb_id, acc, wt, pos, mut) is None:
        sys.exit(2)

    ph_pdb = protonate(pdb_id, a.ph)
    atoms = read_pdb(ph_pdb)
    inside, groups, caps, charge, thio = build_cluster(atoms, chain, resnum, a.radius)
    print(f"[cluster] {len(inside)} residues, {len(caps)} capping H, "
          f"net charge {charge:+d} (spec said {t['net_charge']:+d})")
    print("[cluster] " + " ".join(sorted(f"{k[2]}{k[1]}" for k in inside)))
    if charge != t["net_charge"]:
        print("[cluster] NOTE: charge differs from the spec -- the spec assumed "
              "standard pH 7.4 states; this run used the protonation pdbfixer "
              "produced. Reconcile before quoting a number.")
    xyz = f"{a.target.lower()}_cluster.xyz"
    write_xyz(inside, groups, caps, xyz, comment=f"{a.target} {pdb_id} r={a.radius}")
    print(f"[cluster] wrote {xyz}")
    if a.cluster_only:
        return

    mf = build_mf(xyz, charge, t["spin_2S"], basis=a.basis)
    results = []
    for th in ([float(x) for x in a.sweep.split(",")] if a.sweep else [a.threshold]):
        results.append(avas_at_threshold(mf, t["avas_ao_set"], th))
    out = f"{a.target.lower()}_avas.json"
    json.dump(dict(target=a.target, pdb=pdb_id, radius=a.radius,
                   basis=a.basis, charge=charge, spin=t["spin_2S"],
                   ao_set=t["avas_ao_set"], runs=results),
              open(out, "w"), indent=1)
    print(f"\nwrote {out}")
    if len(results) > 1:
        print("threshold  CAS(nelec,ncas)  qubits")
        for r in results:
            print(f"{r['threshold']:<10}CAS({r['nelecas']},{r['ncas']})"
                  f"{r['qubits']:>12}")


if __name__ == "__main__":
    main()
