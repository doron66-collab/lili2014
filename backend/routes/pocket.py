"""Geometric pocket detection (fpocket) — the SECOND, independent axis of
SOLANGE's target assessment, alongside the existing electronic-structure/DMRG
classification (Class A/B/C). A mutation can be Class B (classically
tractable) and still have no 3D cavity for a drug-like ligand to occupy —
this endpoint answers that second question, on demand, for a user-selected
structure.

Logic mirrors scripts/laguna/run_pocket_detect.py (built and verified first
on TP53 2OCJ — see that script's own module docstring for the crystal-contact
artifact finding this endpoint also flags). Duplicated here rather than
imported across the backend/scripts boundary: Render's deploy only ships
whatever the service's root directory holds, and that boundary has already
drifted once elsewhere in this codebase (see pdb.py's MUTATION_PDB_MAP
comment) — a direct cross-directory import would be the same risk again,
silently fragile to how the two are packaged.

fpocket itself has no apt/pip package; render.yaml's buildCommand compiles it
from source into backend/bin/fpocket. If that build step fails (e.g. no
gcc/make in Render's Python image), every call here returns 503 rather than
crashing — this is a real, not-yet-fully-verified deploy step (see CLAUDE.md
Standing Tasks discipline: flag it, don't pretend it's certain).
"""
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

from fastapi import APIRouter, HTTPException

router = APIRouter()

_BACKEND_DIR = Path(__file__).parent.parent
_CACHE_DIR = _BACKEND_DIR / "pocket_cache"
_STRUCT_CACHE = _BACKEND_DIR / "pocket_struct_cache"


def _find_fpocket_bin():
    candidates = [
        os.environ.get("FPOCKET_BIN"),
        str(_BACKEND_DIR / "bin" / "fpocket"),
        shutil.which("fpocket"),
    ]
    for c in candidates:
        if c and os.path.exists(c) and os.access(c, os.X_OK):
            return c
    return None


def _fetch_pdb(pdb_id: str) -> str:
    """Plain-.pdb fetch only (no gemmi/.cif fallback — see module docstring
    for why that's a deliberately smaller scope than the Laguna script)."""
    import urllib.request
    _STRUCT_CACHE.mkdir(parents=True, exist_ok=True)
    pdb_path = _STRUCT_CACHE / f"{pdb_id.upper()}.pdb"
    if pdb_path.exists():
        return str(pdb_path)
    try:
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", pdb_path)
        return str(pdb_path)
    except Exception as e:
        raise HTTPException(
            status_code=502,
            detail=f"could not fetch {pdb_id.upper()}.pdb from RCSB (no legacy PDB-format file for "
                   f"large/cryo-EM entries is not handled by this endpoint yet): {e}",
        )


def _run_fpocket(fpocket_bin: str, pdb_path: str, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    local_pdb = out_dir / os.path.basename(pdb_path)
    if not local_pdb.exists():
        shutil.copy(pdb_path, local_pdb)
    try:
        r = subprocess.run(
            [fpocket_bin, "-f", local_pdb.name],
            cwd=out_dir, capture_output=True, text=True, timeout=120,
        )
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="fpocket did not finish within 120s")
    if r.returncode != 0:
        raise HTTPException(status_code=500, detail=f"fpocket exited {r.returncode}: {r.stderr[-800:]}")
    stem = local_pdb.stem
    info_path = out_dir / f"{stem}_out" / f"{stem}_info.txt"
    atm_dir = out_dir / f"{stem}_out" / "pockets"
    if not info_path.exists():
        raise HTTPException(status_code=500, detail="fpocket ran but produced no _info.txt")
    return info_path, atm_dir, stem


def _parse_info_txt(info_path: Path):
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


def _pocket_atoms(atm_dir: Path, pocket_id: int):
    """Returns (chain_set, centroid_xyz, atom_count) for one pocket's own
    lining-atom file. The centroid + a volume-derived radius is what the 3D
    viewer renders as a sphere — cheaper and more robust than exposing
    fpocket's raw alpha-sphere vertex file, and sufficient to show WHERE a
    candidate cavity sits on the structure."""
    path = atm_dir / f"pocket{pocket_id}_atm.pdb"
    if not path.exists():
        return set(), None, 0
    chains = set()
    xs, ys, zs = [], [], []
    for line in open(path):
        if line[:6] in ("ATOM  ", "HETATM") and len(line) > 54:
            ch = line[21].strip()
            if ch:
                chains.add(ch)
            try:
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
            except ValueError:
                continue
    if not xs:
        return chains, None, 0
    centroid = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    return chains, centroid, len(xs)


def _residue_near(atm_dir: Path, pocket_id: int, target_chain: str, target_resnum: int):
    path = atm_dir / f"pocket{pocket_id}_atm.pdb"
    if not path.exists():
        return False
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM") or len(line) < 26:
            continue
        if line[21].strip() != target_chain:
            continue
        try:
            if int(line[22:26].strip()) == target_resnum:
                return True
        except ValueError:
            continue
    return False


@router.get("/detect")
def detect_pockets(
    pdb: str,
    near_residue: str | None = None,
    min_druggability: float = 0.5,
):
    """On-demand geometric pocket detection for a user-selected PDB entry —
    the interactive entry point backing the "FIND POCKETS" action in the 3D
    viewer (NSCLCViewer.tsx), triggered by the user's own choice of target,
    not a fixed list.

    near_residue: "chain:resnum" (e.g. "C:275") — flags which pockets
    include that residue and reports it in the response; omit to see the
    whole structure's pocket landscape.
    """
    pdb_id = pdb.strip().upper()
    if not re.match(r"^[0-9][A-Z0-9]{3}$", pdb_id):
        raise HTTPException(status_code=400, detail=f"{pdb_id!r} doesn't look like a PDB ID (e.g. 2OCJ)")

    target_chain = target_resnum = None
    if near_residue:
        try:
            target_chain, resnum_s = near_residue.split(":")
            target_resnum = int(resnum_s)
        except ValueError:
            raise HTTPException(status_code=400, detail="near_residue must be 'chain:resnum', e.g. C:275")

    cache_key = f"{pdb_id}_{near_residue or 'all'}_{min_druggability}".replace(":", "-")
    cache_path = _CACHE_DIR / f"{cache_key}.json"
    if cache_path.exists():
        import json
        return json.loads(cache_path.read_text())

    fpocket_bin = _find_fpocket_bin()
    if not fpocket_bin:
        raise HTTPException(
            status_code=503,
            detail="pocket detection unavailable on this deploy — fpocket binary not found "
                   "(backend/bin/fpocket missing; check the Render build log for the fpocket "
                   "compile step in render.yaml)",
        )

    pdb_path = _fetch_pdb(pdb_id)
    out_dir = _CACHE_DIR / "_runs" / pdb_id
    info_path, atm_dir, _stem = _run_fpocket(fpocket_bin, pdb_path, out_dir)
    pockets = _parse_info_txt(info_path)

    for p in pockets:
        chains, centroid, _n = _pocket_atoms(atm_dir, p["pocket_id"])
        p["chains"] = sorted(chains)
        p["single_chain"] = len(chains) == 1
        if centroid:
            p["x"], p["y"], p["z"] = centroid
            vol = p.get("volume") or 0
            p["radius"] = round((3 * vol / (4 * math.pi)) ** (1 / 3), 2) if vol > 0 else 4.0
        if target_chain:
            p["includes_target_residue"] = _residue_near(atm_dir, p["pocket_id"], target_chain, target_resnum)

    candidates = [p for p in pockets if p.get("single_chain") and p.get("druggability_score", 0) >= min_druggability]
    candidates.sort(key=lambda p: p.get("score", -999), reverse=True)
    best = candidates[0] if candidates else None

    result = dict(
        pdb=pdb_id,
        near_residue=near_residue,
        min_druggability=min_druggability,
        n_pockets_total=len(pockets),
        n_single_chain_druggable=len(candidates),
        best_single_chain_pocket=best,
        all_pockets=pockets,
    )

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    import json
    cache_path.write_text(json.dumps(result))
    return result
