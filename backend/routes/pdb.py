"""
PDB proxy — fetches real protein structure data from RCSB PDB REST API.
Also proxies UniProt gene lookup and AlphaFold predictions.
Caches responses for 24 h to avoid hammering external APIs.
"""
import time
from functools import lru_cache

import httpx
from fastapi import APIRouter, HTTPException

router = APIRouter()

RCSB_BASE      = "https://data.rcsb.org/rest/v1/core"
RCSB_SEARCH    = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_DOWNLOAD   = "https://files.rcsb.org/download"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API  = "https://alphafold.ebi.ac.uk/api/prediction"

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


# ── UniProt gene lookup ────────────────────────────────────────────────────

@lru_cache(maxsize=64)
def _uniprot_gene(gene_symbol: str, _ttl: int) -> dict:
    params = {
        "query":  f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,protein_name,gene_names,xref_pdb",
        "format": "json",
        "size":   "1",
    }
    with httpx.Client(timeout=15.0) as client:
        r = client.get(UNIPROT_SEARCH, params=params)
        r.raise_for_status()
        return r.json()


@router.get("/lookup/gene/{gene_symbol}")
async def lookup_gene(gene_symbol: str):
    """Look up a human gene in UniProt — returns UniProt accession + linked PDB IDs."""
    ttl = int(time.time()) // 86400
    try:
        data = _uniprot_gene(gene_symbol.upper(), ttl)
    except httpx.RequestError as e:
        raise HTTPException(status_code=502, detail=f"UniProt unreachable: {e}")

    results = data.get("results", [])
    if not results:
        return {"gene": gene_symbol, "found": False, "uniprot_id": None,
                "protein_name": None, "pdb_ids": [], "alphafold_url": None}

    entry      = results[0]
    accession  = entry.get("primaryAccession", "")
    prot_name  = (entry.get("proteinDescription", {})
                       .get("recommendedName", {})
                       .get("fullName", {})
                       .get("value", ""))
    pdb_ids = [x["id"] for x in entry.get("uniProtKBCrossReferences", [])
               if x.get("database") == "PDB"]

    return {
        "gene":         gene_symbol.upper(),
        "found":        True,
        "uniprot_id":   accession,
        "protein_name": prot_name,
        "pdb_ids":      pdb_ids[:10],
        "rcsb_url":     f"https://www.rcsb.org/search?query={accession}",
        "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{accession}",
        "alphafold_model_url": f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb",
        "_source": "UniProt REST API",
    }


# ── AlphaFold lookup ───────────────────────────────────────────────────────

@lru_cache(maxsize=64)
def _alphafold_fetch(uniprot_id: str, _ttl: int) -> list:
    with httpx.Client(timeout=15.0) as client:
        r = client.get(f"{ALPHAFOLD_API}/{uniprot_id}")
        if r.status_code == 404:
            return []
        r.raise_for_status()
        return r.json()


@router.get("/alphafold/{uniprot_id}")
async def get_alphafold(uniprot_id: str):
    """Return AlphaFold prediction metadata for a UniProt accession."""
    ttl = int(time.time()) // 86400
    try:
        data = _alphafold_fetch(uniprot_id.upper(), ttl)
    except httpx.RequestError as e:
        raise HTTPException(status_code=502, detail=f"AlphaFold DB unreachable: {e}")

    if not data:
        return {"uniprot_id": uniprot_id, "available": False}

    pred = data[0]
    return {
        "uniprot_id":  uniprot_id,
        "available":   True,
        "entry_id":    pred.get("entryId"),
        "pLDDT":       pred.get("globalMetricValue"),
        "organism":    pred.get("organismScientificName"),
        "model_url":   pred.get("pdbUrl"),
        "view_url":    f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}",
        "_source":     "AlphaFold DB · EMBL-EBI",
    }
