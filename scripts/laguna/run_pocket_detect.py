#!/usr/bin/env python3
"""Geometric pocket detection (fpocket) -- the second axis of SOLANGE's
classification, independent of the electronic-structure/DMRG side.

Context (2026-10-07): SOLANGE computes the FIRST axis -- is this active
site's electronic structure classically tractable (Class B/C) or does it
need quantum hardware (Class A). That answer is already in hand or close
to it for most of this project's targets. The SECOND, independent axis --
does a 3D binding pocket even exist for a drug-like ligand to occupy -- had
NO computational tool in this codebase before this script; the existing
3D viewer (PDBMolViewer.tsx, NGL.js) only renders structure, it does not
detect cavities (confirmed by grep: the only "pocket" reference in the
whole frontend is a cosmetic NGL selection string for highlighting nearby
residues, not a geometric pocket-finding algorithm).

Verified against a known case before trusting this on anything new (DP1,
"verify, don't trust"): run on TP53 2OCJ (the project's own anchor
structure), this tool found 59 candidate cavities, 58 of them scoring
near-zero druggability (<0.1) -- consistent with the literature claim that
TP53's DNA-binding domain has no conventional druggable pocket (it's a
transcription factor, not an enzyme). The one high-scoring exception
(Druggability Score 0.774) turned out to be a crystal-contact artifact --
spans two separate chains (A and C) of the same asymmetric unit, not a
real binding site -- and fpocket's own combined Score (not the isolated
druggability sub-score) correctly flagged it as the worst-ranked pocket of
the 59 (-0.266). This script reports BOTH scores and flags multi-chain
pockets explicitly, rather than reporting the single Druggability Score
number alone, which would have been actively misleading in this exact
case.

Setup (fpocket has no apt/pip package; build from source once):
    git clone --depth 1 https://github.com/Discngine/fpocket.git
    cd fpocket && make
    # binary lands at fpocket/bin/fpocket -- point --fpocket-bin at it,
    # or put it on PATH

Usage
-----
    python run_pocket_detect.py --pdb 2OCJ
    python run_pocket_detect.py --pdb 2OCJ --near-residue C:275     # TP53 C275F site
    python run_pocket_detect.py --local-file struct/8GS8.pdb --fpocket-bin /path/to/fpocket

Output: <key>_pockets.json -- every pocket's full descriptor set, each one
flagged single_chain=True/False, plus a top-level best_single_chain_pocket
summary (None if no single-chain pocket clears a minimal bar). Does NOT
pick a winner for you past that flag -- a human (or a further-downstream
gate, per this project's own DP5/DP6 discipline) still has to read the
result, the same way every other classification in this pipeline refuses
to self-declare a grade it isn't entitled to.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys


def fetch_pdb(pdb_id, cache="struct"):
    """Plain-PDB-format fetch -- fpocket's bundled mmCIF reader exists but
    the plain .pdb path is the most-tested one; RCSB serves both for most
    entries. Falls back to converting a cached .cif via gemmi (same pattern
    as run_gate2_sdhb.py's cif_to_pdb()) if no legacy .pdb exists for a
    cryo-EM entry."""
    import urllib.request
    os.makedirs(cache, exist_ok=True)
    pdb_path = os.path.join(cache, f"{pdb_id}.pdb")
    if os.path.exists(pdb_path):
        return pdb_path
    try:
        print(f"[fetch] downloading {pdb_path} from RCSB...")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", pdb_path)
        return pdb_path
    except Exception as e:
        print(f"[fetch] plain .pdb fetch failed ({e}) -- falling back to .cif conversion "
              f"(expected for large/cryo-EM entries like 8GS8/9KC4)")
    cif_path = os.path.join(cache, f"{pdb_id}.cif")
    if not os.path.exists(cif_path):
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.cif", cif_path)
    import gemmi
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    st.write_pdb(pdb_path)
    print(f"[fetch] converted {cif_path} -> {pdb_path} via gemmi")
    return pdb_path


def find_fpocket_bin(explicit):
    if explicit:
        if not os.path.exists(explicit):
            sys.exit(f"--fpocket-bin {explicit} does not exist")
        return explicit
    found = shutil.which("fpocket")
    if found:
        return found
    sys.exit("fpocket binary not found on PATH and --fpocket-bin not given -- "
             "see this script's own module docstring for the one-time build steps "
             "(git clone + make, no apt/pip package exists for it)")


def run_fpocket(fpocket_bin, pdb_path, out_dir):
    """Runs fpocket; it writes its own <name>_out/ directory next to the
    input file (not controllable via a flag in fpocket 4.x), so this copies
    the input into out_dir first to keep output location predictable."""
    os.makedirs(out_dir, exist_ok=True)
    local_pdb = os.path.join(out_dir, os.path.basename(pdb_path))
    if os.path.abspath(local_pdb) != os.path.abspath(pdb_path):
        shutil.copy(pdb_path, local_pdb)
    print(f"[fpocket] running on {local_pdb}...")
    r = subprocess.run([fpocket_bin, "-f", os.path.basename(local_pdb)],
                        cwd=out_dir, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fpocket exited {r.returncode}\nstdout:\n{r.stdout}\nstderr:\n{r.stderr}")
    stem = os.path.splitext(os.path.basename(local_pdb))[0]
    info_path = os.path.join(out_dir, f"{stem}_out", f"{stem}_info.txt")
    atm_dir = os.path.join(out_dir, f"{stem}_out", "pockets")
    if not os.path.exists(info_path):
        sys.exit(f"fpocket ran but {info_path} was not created -- check stdout/stderr above")
    return info_path, atm_dir, stem


def parse_info_txt(info_path):
    """Parses fpocket's <name>_info.txt into a list of dicts, one per
    pocket, preserving every descriptor fpocket reports (not just the two
    this script cares about most) -- a future question about e.g. Volume
    or Hydrophobicity score shouldn't need re-parsing logic written again."""
    pockets = []
    cur = None
    for line in open(info_path):
        m = re.match(r"Pocket\s+(\d+)\s*:", line)
        if m:
            if cur:
                pockets.append(cur)
            cur = {"pocket_id": int(m.group(1))}
            continue
        m = re.match(r"\s*(.+?)\s*:\s*\t?\s*([-\d.]+)\s*$", line)
        if m and cur is not None:
            key = re.sub(r"[^a-z0-9]+", "_", m.group(1).strip().lower()).strip("_")
            try:
                cur[key] = float(m.group(2))
            except ValueError:
                pass
    if cur:
        pockets.append(cur)
    return pockets


def chains_in_pocket(atm_dir, pocket_id):
    """Reads pocket<N>_atm.pdb (the atoms fpocket assigned to this pocket)
    and returns the set of chain IDs involved. A pocket spanning more than
    one chain in a multi-copy asymmetric unit is almost always a crystal
    contact, not a real binding site -- see this script's own module
    docstring for the TP53 2OCJ case that motivated this check."""
    path = os.path.join(atm_dir, f"pocket{pocket_id}_atm.pdb")
    if not os.path.exists(path):
        return None
    chains = set()
    for line in open(path):
        if line[:6] in ("ATOM  ", "HETATM") and len(line) > 21:
            ch = line[21].strip()
            if ch:
                chains.add(ch)
    return chains


def residue_near(atm_dir, pocket_id, target_chain, target_resnum, radius_atoms=True):
    """Cheap proximity check: does this pocket's own atom set include the
    target residue, or (if not an exact hit) report the pocket's own
    residue range on the target chain for a human to eyeball distance from.
    Not a substitute for an actual distance calculation if precision matters
    -- good enough to triage 50+ pockets down to the ones worth a closer
    look near a specific mutation site."""
    path = os.path.join(atm_dir, f"pocket{pocket_id}_atm.pdb")
    if not os.path.exists(path):
        return False, []
    hit = False
    resnums_on_chain = set()
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM") or len(line) < 26:
            continue
        ch = line[21].strip()
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        if ch == target_chain:
            resnums_on_chain.add(resnum)
            if resnum == target_resnum:
                hit = True
    return hit, sorted(resnums_on_chain)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--pdb", help="PDB ID to fetch (e.g. 2OCJ)")
    g.add_argument("--local-file", help="path to an already-downloaded .pdb/.cif")
    ap.add_argument("--near-residue", default=None,
                     help="chain:resnum (e.g. C:275) -- if given, flags which pockets "
                          "include or sit on the same chain near this residue. Optional; "
                          "omit to just see the whole protein's pocket landscape.")
    ap.add_argument("--fpocket-bin", default=None,
                     help="path to the fpocket binary; omit to search PATH")
    ap.add_argument("--min-druggability", type=float, default=0.5,
                     help="fpocket's own conventional druggable/undruggable cutoff")
    ap.add_argument("--out-dir", default="pocket_detect_out")
    a = ap.parse_args()

    fpocket_bin = find_fpocket_bin(a.fpocket_bin)

    if a.pdb:
        pdb_path = fetch_pdb(a.pdb)
        key = a.pdb
    else:
        pdb_path = a.local_file
        key = os.path.splitext(os.path.basename(pdb_path))[0]

    info_path, atm_dir, stem = run_fpocket(fpocket_bin, pdb_path, a.out_dir)
    pockets = parse_info_txt(info_path)
    print(f"[fpocket] {len(pockets)} candidate pocket(s) found")

    target_chain = target_resnum = None
    if a.near_residue:
        target_chain, target_resnum = a.near_residue.split(":")
        target_resnum = int(target_resnum)

    for p in pockets:
        chains = chains_in_pocket(atm_dir, p["pocket_id"])
        p["chains"] = sorted(chains) if chains else []
        p["single_chain"] = len(p["chains"]) == 1
        if not p["single_chain"] and chains:
            print(f"[fpocket] NOTE: pocket {p['pocket_id']} spans chains {p['chains']} -- "
                  f"likely a crystal-contact artifact (see module docstring), not flagged "
                  f"as a real candidate below regardless of its druggability score")
        if target_chain:
            hit, resnums = residue_near(atm_dir, p["pocket_id"], target_chain, target_resnum)
            p["includes_target_residue"] = hit
            p["target_chain_resnum_range"] = (min(resnums), max(resnums)) if resnums else None

    # Rank candidates: single-chain only (excludes the exact artifact pattern
    # found on 2OCJ), above the conventional druggability cutoff, by fpocket's
    # own combined score (not the druggability sub-score alone -- see module
    # docstring for why that alone would have been misleading).
    candidates = [p for p in pockets
                  if p.get("single_chain") and p.get("druggability_score", 0) >= a.min_druggability]
    candidates.sort(key=lambda p: p.get("score", -999), reverse=True)
    best = candidates[0] if candidates else None

    if a.near_residue:
        near_candidates = [p for p in candidates if p.get("includes_target_residue")]
        print(f"\n[summary] {len(candidates)} single-chain pocket(s) above druggability "
              f"{a.min_druggability}; {len(near_candidates)} of those actually include "
              f"residue {a.near_residue}")
    else:
        print(f"\n[summary] {len(candidates)} single-chain pocket(s) above druggability "
              f"{a.min_druggability} (out of {len(pockets)} total candidates)")

    if best:
        print(f"[summary] best candidate: pocket {best['pocket_id']}, "
              f"druggability={best.get('druggability_score')}, score={best.get('score')}, "
              f"volume={best.get('volume')}, chains={best['chains']}")
    else:
        print(f"[summary] NO pocket clears the single-chain + druggability>={a.min_druggability} "
              f"bar -- consistent with a genuinely non-druggable target by this criterion alone "
              f"(not a tool failure; this is exactly what TP53 2OCJ itself returns)")

    out = f"{key}_pockets.json"
    json.dump(dict(
        key=key, pdb_source=pdb_path, n_pockets_total=len(pockets),
        near_residue=a.near_residue, min_druggability=a.min_druggability,
        best_single_chain_pocket=best, all_pockets=pockets,
    ), open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
