"""
PDB proxy — fetches real protein structure data from RCSB PDB REST API.
Caches responses for 24 h to avoid hammering the external API.
"""
import time
from functools import lru_cache

import httpx
from fastapi import APIRouter, HTTPException

router = APIRouter()

RCSB_BASE = "https://data.rcsb.org/rest/v1/core"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_DOWNLOAD = "https://files.rcsb.org/download"

# PDB IDs relevant to the five NSCLC mutations in this study
MUTATION_PDB_MAP = {
    "TP53_C275F":  "2OCJ",
    "TP53_Y220C":  "2VUK",
    "KEAP1_LOF":   "2FLU",
    "STK11_LKB1":  "2QK7",
    "CDKN2A_P16":  "2A5E",
}


@lru_cache(maxsize=32)
def _cached_fetch(url: str, _ttl_bucket: int) -> dict:
    """Fetch JSON from URL; _ttl_bucket buckets time into 24-h windows for cache TTL."""
    with httpx.Client(timeout=15.0) as client:
        r = client.get(url)
        r.raise_for_status()
        return r.json()


def fetch_pdb_entry(pdb_id: str) -> dict:
    ttl = int(time.time()) // 86400
    return _cached_fetch(f"{RCSB_BASE}/entry/{pdb_id.upper()}", ttl)


def fetch_pdb_polymer_entities(pdb_id: str) -> list:
    ttl = int(time.time()) // 86400
    data = _cached_fetch(f"{RCSB_BASE}/entry/{pdb_id.upper()}/polymer_entities", ttl)
    return data if isinstance(data, list) else []


@router.get("/{pdb_id}")
async def get_pdb_entry(pdb_id: str):
    """Return RCSB PDB entry metadata for a given PDB ID."""
    try:
        data = fetch_pdb_entry(pdb_id)
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code,
                            detail=f"RCSB PDB returned {e.response.status_code} for {pdb_id}")
    except httpx.RequestError as e:
        raise HTTPException(status_code=502, detail=f"Could not reach RCSB PDB: {e}")

    entry = data.get("entry", {})
    struct = data.get("struct", {})
    exptl = data.get("exptl", [{}])[0] if data.get("exptl") else {}
    cell = data.get("cell", {})
    return {
        "pdb_id":     pdb_id.upper(),
        "title":      struct.get("title", ""),
        "method":     exptl.get("method", ""),
        "resolution": data.get("refine", [{}])[0].get("ls_d_res_high") if data.get("refine") else None,
        "deposition_date": entry.get("rcsb_entry_info", {}).get("deposition_date"),
        "polymer_entity_count": data.get("rcsb_entry_info", {}).get("polymer_entity_count_protein", 0),
        "atom_count":  data.get("rcsb_entry_info", {}).get("deposited_atom_count", 0),
        "source_organism": (
            data.get("rcsb_entry_info", {}).get("rcsb_source_organism", [{}])[0]
            .get("ncbi_scientific_name", "")
            if data.get("rcsb_entry_info", {}).get("rcsb_source_organism") else ""
        ),
        "cell_dimensions": {
            "a": cell.get("length_a"), "b": cell.get("length_b"), "c": cell.get("length_c"),
            "alpha": cell.get("angle_alpha"), "beta": cell.get("angle_beta"), "gamma": cell.get("angle_gamma"),
        } if cell else None,
        "_source": "RCSB PDB REST API v1",
    }


@router.get("/{pdb_id}/structure-url")
async def get_structure_url(pdb_id: str, format: str = "cif"):
    """Return the download URL for a PDB structure file (cif or pdb)."""
    fmt = format.lower()
    if fmt not in ("cif", "pdb", "mmtf"):
        raise HTTPException(status_code=400, detail="format must be cif, pdb, or mmtf")
    ext = {"cif": ".cif", "pdb": ".pdb", "mmtf": ".mmtf"}[fmt]
    return {
        "pdb_id": pdb_id.upper(),
        "format": fmt,
        "url":    f"{PDB_DOWNLOAD}/{pdb_id.upper()}{ext}",
    }


@router.get("/mutation/{mutation_id}")
async def get_mutation_pdb(mutation_id: str):
    """Return PDB entry for a mutation ID used in this study."""
    pdb_id = MUTATION_PDB_MAP.get(mutation_id)
    if not pdb_id:
        raise HTTPException(status_code=404,
                            detail=f"Unknown mutation: {mutation_id}. "
                                   f"Valid: {list(MUTATION_PDB_MAP.keys())}")
    return await get_pdb_entry(pdb_id)


@router.get("/")
async def list_study_pdbs():
    """List all PDB structures used in this study."""
    return {
        "mutations": [
            {"mutation_id": k, "pdb_id": v,
             "structure_url": f"{PDB_DOWNLOAD}/{v}.cif"}
            for k, v in MUTATION_PDB_MAP.items()
        ]
    }
