#!/usr/bin/env python3
"""
build_qm_cluster.py — extract a QM-ready active-site cluster from a PDB
structure for solange_dmrg.py's --geometry/--avas path.

WHAT THIS DOES (mechanical, deterministic):
  1. Parse ATOM/HETATM records from a local PDB file.
  2. Locate a center point (a residue's CA, a named atom, or a HETATM like a
     metal ion).
  3. Select every residue with at least one atom within --radius of that
     center (whole-residue selection — no partial sidechains).
  4. Where a selected residue's backbone neighbor (resseq-1 / resseq+1, same
     chain) was NOT selected, cap the cut N-C(prev) or C-N(next) bond with a
     hydrogen placed along the original bond direction at a standard 1.09 A
     C-H length — the common "truncate at the peptide bond, cap with H"
     scheme used for QM-cluster models.
  5. Write an .xyz file, ready for solange_dmrg.py --geometry.

WHAT THIS DOES NOT DO (still needs a scientist's judgment call — see the
project's own honesty note in solange_dmrg.py: "the CLUSTER itself ... is
the chemist's input, NOT auto-generated here"):
  - Decide whether --radius is the right cutoff for a given question.
  - Decide whether the metal ion / cofactor should be in scope.
  - Assign protonation states, tautomers, or spin state.
  - Model a POINT MUTATION (e.g. Cys275Phe) — do that in the structure
    BEFORE running this script (PyMOL's Mutagenesis wizard is the standard,
    well-documented way: load the PDB, `wizard mutagenesis`, pick the
    residue, choose CYS -> PHE, pick the lowest-clash rotamer, save). This
    script then extracts whatever structure you hand it — native or mutant,
    it does not know or care which.
  - Energy-minimize the capped cluster. The capping H position is a
    reasonable starting geometry, not a relaxed one; if solange_dmrg.py's
    CASSCF step converges poorly, a short classical (MM or DFT) relaxation
    of the cap positions first is the standard fix.

USAGE:
  # center on a residue's CA
  python build_qm_cluster.py --pdb 1TSR.pdb --chain A --resi 275 \\
      --radius 6.0 --out tp53_c275_cluster.xyz

  # center on a HETATM (e.g. the structural zinc ion)
  python build_qm_cluster.py --pdb 1TSR.pdb --hetero ZN --radius 6.0 \\
      --out tp53_zn_cluster.xyz

  # after PyMOL mutagenesis, same command on the mutant file:
  python build_qm_cluster.py --pdb 1TSR_C275F.pdb --chain A --resi 275 \\
      --radius 6.0 --out tp53_c275f_cluster.xyz

Then: python solange_dmrg.py --geometry tp53_c275f_cluster.xyz \\
          --avas "<see printed suggestion>" --key TP53_C275F --submit
"""

import argparse
import json
import re
import shlex
import sys
import urllib.request
from pathlib import Path

import numpy as np

RCSB_DOWNLOAD = "https://files.rcsb.org/download"
_HERE = Path(__file__).resolve()


def _targets_path():
    """targets.json is the single source of truth for which structure a target
    uses (CLAUDE.md). Searched rather than assumed, same as simulate.py and
    pdb.py do, so this works from a checkout that ships only one of the two."""
    for p in (_HERE.parents[2] / "targets.json", _HERE.parents[2] / "backend" / "targets.json"):
        if p.exists():
            return p
    return None


def resolve_target(key):
    """Look up a mutation key in targets.json -> (pdb_id, residue_number).

    Returns the residue number ONLY when the target actually names one. A
    generic loss-of-function key (KEAP1_LOF) has no single residue by
    definition — the dissertation records those genes as having "no single-
    residue anchor" — so this returns None for it rather than inventing a
    center. The caller must then supply --resi or --hetero explicitly.

    Chain is deliberately NOT guessed: targets.json stores no chain ID, and
    picking one (say "A") would be a silent structural assumption of exactly
    the kind this pipeline must not make.
    """
    p = _targets_path()
    if p is None:
        raise SystemExit("targets.json not found — pass --pdb/--pdb-id and --resi directly.")
    muts = (json.loads(p.read_text(encoding="utf-8")).get("mutations") or {})
    if key not in muts:
        raise SystemExit(f"'{key}' is not in targets.json. Known: {', '.join(sorted(muts))}")
    entry = muts[key]
    pdb_id = entry.get("pdb")
    if not pdb_id:
        raise SystemExit(f"{key} has no 'pdb' field in targets.json.")
    # Residue number from the canonical name ("TP53 p.Cys275Phe" -> 275); fall
    # back to the key's own suffix ("TP53_C275F" -> 275). Both are absent for a
    # generic LOF target, which is the honest answer, not a failure.
    resi = None
    m = re.search(r"p\.[A-Za-z]{3}(\d+)", entry.get("name") or "")
    if m:
        resi = int(m.group(1))
    else:
        m = re.match(r"^[A-Z0-9]+_[A-Z](\d+)[A-Z]$", key)
        if m:
            resi = int(m.group(1))
    return pdb_id, resi


def fetch_structure(pdb_id, cache_dir="./pdb_cache"):
    """Download {pdb_id}.cif from RCSB into cache_dir, or reuse it if present.

    mmCIF rather than legacy PDB: the .pdb format cannot represent structures
    past 62 chains or 99,999 atoms and RCSB no longer issues it for large
    entries, so requesting it would work for today's small targets and fail
    silently-later on a bigger one. parse_structure() already handles both.
    """
    pdb_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise SystemExit(f"'{pdb_id}' is not a 4-character PDB ID.")
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / f"{pdb_id}.cif"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"using cached {dest}")
        return str(dest)
    url = f"{RCSB_DOWNLOAD}/{pdb_id}.cif"
    print(f"downloading {url} …")
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read()
    except Exception as e:
        raise SystemExit(f"could not download {pdb_id} from RCSB ({e}). "
                         f"Download it manually and pass --pdb <file>.")
    dest.write_bytes(data)
    print(f"wrote {dest} ({len(data)/1024:.0f} KB)")
    return str(dest)

# Standard C-H capping bond length (Angstrom).
CAP_BOND_LENGTH = 1.09

# Element inference fallback when columns 77-78 are blank (older PDB files).
_NAME_ELEMENT_HINTS = {"CA": "C", "CB": "C", "N": "N", "O": "O", "S": "S",
                        "ZN": "ZN", "MG": "MG", "FE": "FE", "MN": "MN"}


def parse_structure(path):
    """Dispatch to the PDB or mmCIF parser by file extension. Both return the
    same atom-dict shape so the rest of the script doesn't care which format
    it got — RCSB serves legacy .pdb only for older depositions, .cif always."""
    suffix = Path(path).suffix.lower()
    if suffix == ".cif":
        return parse_cif(path)
    return parse_pdb(path)


def parse_cif(path):
    """Minimal mmCIF _atom_site loop parser — just enough to get coordinates
    + auth_* identifiers (the numbering scheme that matches the PDB literature
    and RCSB's own web viewer, as opposed to label_* which renumbers from 1
    per entity). Not a general CIF parser: assumes one _atom_site loop, values
    on one line each (true for every RCSB-deposited structure)."""
    lines = Path(path).read_text().splitlines()

    # Phase 1: locate the _atom_site loop's tag block. Find a "loop_" line
    # immediately followed by "_atom_site." tags (there is exactly one such
    # loop in a standard RCSB mmCIF file).
    tag_start = None
    for i, raw in enumerate(lines):
        if raw.strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].strip().startswith("_atom_site."):
            tag_start = i + 1
            break
    if tag_start is None:
        raise SystemExit(f"no _atom_site loop found in {path} — is this a valid mmCIF file?")

    tags = []
    j = tag_start
    while j < len(lines) and lines[j].strip().startswith("_atom_site."):
        tags.append(lines[j].strip().split(".", 1)[1])
        j += 1

    # Phase 2: data rows run from here until a blank line, a comment ("#"),
    # a new tag/loop_, or EOF — whichever comes first.
    rows = []
    while j < len(lines):
        line = lines[j].strip()
        if not line or line.startswith("#") or line.startswith("_") or line == "loop_":
            break
        rows.append(shlex.split(line))
        j += 1

    if not tags or not rows:
        raise SystemExit(f"no _atom_site loop found in {path} — is this a valid mmCIF file?")
    idx = {t: n for n, t in enumerate(tags)}
    required = ["group_PDB", "auth_atom_id", "auth_comp_id", "auth_asym_id",
                "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol"]
    missing = [t for t in required if t not in idx]
    if missing:
        raise SystemExit(f"mmCIF _atom_site loop is missing expected columns: {missing}")

    atoms = []
    for r in rows:
        if len(r) != len(tags):
            continue  # defensive: malformed/truncated row
        record = r[idx["group_PDB"]]
        if record not in ("ATOM", "HETATM"):
            continue
        alt_idx = idx.get("label_alt_id")
        if alt_idx is not None and r[alt_idx] not in (".", "?", "A"):
            continue  # keep only the primary alternate location, same as parse_pdb
        icode_idx = idx.get("pdbx_PDB_ins_code")
        atoms.append({
            "record": record,
            "serial": int(r[idx["id"]]) if "id" in idx else len(atoms) + 1,
            "name": r[idx["auth_atom_id"]].strip("'\""),
            "resname": r[idx["auth_comp_id"]],
            "chain": r[idx["auth_asym_id"]],
            "resseq": int(r[idx["auth_seq_id"]]),
            "icode": "" if icode_idx is None or r[icode_idx] in (".", "?") else r[icode_idx],
            "xyz": np.array([float(r[idx["Cartn_x"]]), float(r[idx["Cartn_y"]]), float(r[idx["Cartn_z"]])]),
            "element": r[idx["type_symbol"]].upper(),
        })
    if not atoms:
        raise SystemExit(f"parsed the _atom_site loop in {path} but found zero ATOM/HETATM rows")
    return atoms


def parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue  # keep only the primary alternate location
            element = line[76:78].strip().upper()
            name = line[12:16].strip()
            if not element:
                element = _NAME_ELEMENT_HINTS.get(name.upper(), name[0].upper())
            atoms.append({
                "record": line[0:6].strip(),
                "serial": int(line[6:11]),
                "name": name,
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resseq": int(line[22:26]),
                "icode": line[26].strip(),
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "element": element,
            })
    if not atoms:
        raise SystemExit(f"no ATOM/HETATM records parsed from {path} — check the file")
    return atoms


def find_center(atoms, chain=None, resi=None, atom_name=None, hetero=None):
    if hetero:
        hits = [a for a in atoms if a["record"] == "HETATM" and a["resname"] == hetero
                and (chain is None or a["chain"] == chain)]
        if not hits:
            raise SystemExit(f"no HETATM with resname={hetero} found"
                              + (f" in chain {chain}" if chain else ""))
        if len(hits) > 1:
            print(f"WARNING: {len(hits)} {hetero} atoms found; using the first "
                  f"(chain {hits[0]['chain']}, resseq {hits[0]['resseq']}). "
                  f"Pass --chain to disambiguate.", file=sys.stderr)
        return hits[0]["xyz"]
    if chain is None or resi is None:
        raise SystemExit("provide --chain and --resi, or --hetero")
    name = atom_name or "CA"
    hits = [a for a in atoms if a["chain"] == chain and a["resseq"] == resi and a["name"] == name]
    if not hits:
        raise SystemExit(f"atom {name} not found at chain={chain} resi={resi} — "
                          f"check numbering against the PDB file directly")
    return hits[0]["xyz"]


def select_residues(atoms, center, radius):
    by_res = {}
    for a in atoms:
        by_res.setdefault((a["chain"], a["resseq"], a["icode"]), []).append(a)
    selected = {}
    for key, ats in by_res.items():
        if any(np.linalg.norm(a["xyz"] - center) <= radius for a in ats):
            selected[key] = ats
    return selected


def _atom(ats, name):
    for a in ats:
        if a["name"] == name:
            return a
    return None


def cap_dangling_bonds(selected):
    """For each selected polymer residue, if its sequence neighbor (same
    chain, resseq +/-1) is NOT selected, replace the missing backbone bond
    with a capping H placed along the original bond direction."""
    caps = []
    for (chain, resseq, icode), ats in selected.items():
        if (chain, resseq - 1, icode) not in selected:
            n_atom = _atom(ats, "N")
            ca_atom = _atom(ats, "CA")
            if n_atom is not None and ca_atom is not None:
                direction = n_atom["xyz"] - ca_atom["xyz"]
                norm = np.linalg.norm(direction)
                if norm > 1e-6:
                    h_pos = n_atom["xyz"] + (direction / norm) * (CAP_BOND_LENGTH - norm) \
                        if norm < CAP_BOND_LENGTH else n_atom["xyz"] + (direction / norm) * CAP_BOND_LENGTH
                    caps.append(("H", h_pos))
        if (chain, resseq + 1, icode) not in selected:
            c_atom = _atom(ats, "C")
            ca_atom = _atom(ats, "CA")
            if c_atom is not None and ca_atom is not None:
                direction = c_atom["xyz"] - ca_atom["xyz"]
                norm = np.linalg.norm(direction)
                if norm > 1e-6:
                    h_pos = c_atom["xyz"] + (direction / norm) * CAP_BOND_LENGTH
                    caps.append(("H", h_pos))
    return caps


def write_xyz(selected, caps, out_path, include_hetero=True, comment=""):
    rows = []
    for ats in selected.values():
        for a in ats:
            if a["record"] == "HETATM" and not include_hetero:
                continue
            rows.append((a["element"], a["xyz"]))
    rows.extend(caps)
    with open(out_path, "w") as f:
        f.write(f"{len(rows)}\n{comment}\n")
        for element, xyz in rows:
            f.write(f"{element:<3s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}\n")
    return rows


def suggest_avas(rows):
    """Print a starting --avas suggestion from the elements actually present.
    This is a starting point, not a verdict — confirm against what's
    chemically relevant before submitting."""
    elements = sorted({e for e, _ in rows if e not in ("H", "C")})
    hints = {"ZN": "Zn 3d", "MG": "Mg 3s", "FE": "Fe 3d", "MN": "Mn 3d",
             "S": "S 3p", "N": "N 2p", "O": "O 2p"}
    parts = [hints[e] for e in elements if e in hints]
    return ", ".join(parts) if parts else "(no metal/heteroatom in range — review manually)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb", help="local structure file, .pdb or .cif")
    src.add_argument("--pdb-id", help="4-character PDB ID to download from RCSB (cached)")
    src.add_argument("--target", help="mutation key from targets.json (e.g. TP53_C275F) — "
                     "resolves both the PDB ID and, when the target names a specific "
                     "residue, --resi. You must still supply --chain: targets.json "
                     "stores no chain, and guessing one would be a silent structural "
                     "assumption.")
    ap.add_argument("--cache-dir", default="./pdb_cache",
                    help="where downloaded structures are kept (default ./pdb_cache)")
    ap.add_argument("--chain", help="chain ID for the center residue")
    ap.add_argument("--resi", type=int, help="residue number for the center")
    ap.add_argument("--atom", default="CA", help="atom name within that residue to center on (default CA)")
    ap.add_argument("--hetero", help="alternative: center on a HETATM by resname, e.g. ZN")
    ap.add_argument("--radius", type=float, default=6.0, help="selection radius in Angstrom (default 6.0)")
    ap.add_argument("--no-hetero", action="store_true", help="exclude HETATM records (waters, ions, ligands) from the output")
    ap.add_argument("--out", required=True, help="output .xyz path")
    args = ap.parse_args()

    # Resolve the structure source before anything else touches args.pdb.
    if args.target:
        pdb_id, resi = resolve_target(args.target)
        print(f"targets.json: {args.target} -> PDB {pdb_id}"
              + (f", residue {resi}" if resi else ", no single-residue anchor"))
        if resi is not None and args.resi is None and not args.hetero:
            args.resi = resi
        if args.resi is None and not args.hetero:
            raise SystemExit(
                f"{args.target} names no specific residue (a generic loss-of-function "
                f"target has no single-residue anchor), so there is no center to "
                f"extract around. Pass --resi <n> or --hetero <RESNAME> explicitly.")
        if not args.chain and not args.hetero:
            raise SystemExit(
                "--target resolves the structure and residue but NOT the chain — "
                "targets.json does not store one. Pass --chain explicitly (open the "
                "structure and confirm which chain carries this residue).")
        args.pdb = fetch_structure(pdb_id, args.cache_dir)
    elif args.pdb_id:
        args.pdb = fetch_structure(args.pdb_id, args.cache_dir)

    atoms = parse_structure(args.pdb)
    center = find_center(atoms, chain=args.chain, resi=args.resi, atom_name=args.atom, hetero=args.hetero)
    selected = select_residues(atoms, center, args.radius)
    caps = cap_dangling_bonds(selected)

    n_heavy = sum(1 for ats in selected.values() for a in ats
                  if not (a["record"] == "HETATM" and args.no_hetero))
    comment = (f"QM cluster: {args.pdb} center=({args.chain or ''}{args.resi or ''}"
               f"{args.hetero or ''}) radius={args.radius}A residues={len(selected)} "
               f"heavy_atoms={n_heavy} caps={len(caps)}")
    rows = write_xyz(selected, caps, args.out, include_hetero=not args.no_hetero, comment=comment)

    print(f"selected {len(selected)} residues, {n_heavy} heavy atoms, {len(caps)} capping H")
    print(f"wrote {args.out}")
    print(f"suggested --avas: \"{suggest_avas(rows)}\"")
    print("REVIEW BEFORE SUBMITTING: confirm the radius captured what you intend, "
          "check for atoms with only partial sidechains that still needed capping "
          "(rare with whole-residue selection, but verify near proline/glycine), "
          "and confirm protonation state / net charge separately — this script "
          "does not add or remove hydrogens on existing heavy atoms.")


if __name__ == "__main__":
    main()
