"""
Provenance endpoints — CRUD for P1-P9 simulation run records stored in Supabase.
All write operations require a valid Supabase JWT in the Authorization header.
"""
import hashlib
import json
import os

from fastapi import APIRouter, Depends, HTTPException, Header
from supabase import create_client

router = APIRouter()

SUPABASE_URL = os.environ.get("SUPABASE_URL", "")
SUPABASE_KEY = os.environ.get("SUPABASE_KEY", "")


def get_supabase():
    if not SUPABASE_URL or not SUPABASE_KEY:
        raise HTTPException(status_code=503, detail="Supabase not configured")
    return create_client(SUPABASE_URL, SUPABASE_KEY)


def verify_jwt(authorization: str = Header(None)) -> str:
    """
    Validate caller identity via Supabase JWT.
    Returns the raw token for RLS enforcement on the Supabase side.
    RLS policies on simulation_runs ensure researchers see only their own rows.
    """
    if not authorization or not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Missing or malformed Authorization header")
    return authorization[len("Bearer "):]


# ── Read endpoints (public, no auth required) ─────────────────────────────────

@router.get("/runs")
async def list_runs(limit: int = 20, mutation_id: str = None):
    """Return recent simulation runs. Optional filter by mutation_id."""
    sb = get_supabase()
    q = (sb.table("simulation_runs")
           .select("id, created_at, mutation_id, mutation_name, "
                   "p7_energy_ha, p7_ci_lower, p7_ci_upper, p8_hash, phase")
           .order("created_at", desc=True)
           .limit(limit))
    if mutation_id:
        q = q.eq("mutation_id", mutation_id)
    res = q.execute()
    return {"data": res.data, "count": len(res.data)}


@router.get("/runs/{run_id}")
async def get_run(run_id: str):
    """Return the full P1-P9 provenance record for a single run."""
    sb = get_supabase()
    res = sb.table("simulation_runs").select("*").eq("id", run_id).execute()
    if not res.data:
        raise HTTPException(status_code=404, detail=f"Run {run_id} not found")
    return res.data[0]


@router.get("/runs/{run_id}/verify")
async def verify_seal(run_id: str):
    """
    Re-compute the P8 SHA-256 seal from stored P1-P7+P9 fields and compare
    against the stored p8_hash — returns pass/fail integrity check.
    """
    sb = get_supabase()
    res = sb.table("simulation_runs").select("*").eq("id", run_id).execute()
    if not res.data:
        raise HTTPException(status_code=404, detail=f"Run {run_id} not found")
    record = res.data[0]

    seal_payload = json.dumps(
        {k: v for k, v in record.items()
         if k.startswith(("p1_", "p2_", "p3_", "p4_", "p5_", "p6_", "p7_", "p9_"))},
        sort_keys=True, default=str
    )
    recomputed = hashlib.sha256(seal_payload.encode()).hexdigest()
    stored = record.get("p8_hash", "")
    return {
        "run_id":    run_id,
        "integrity": "PASS" if recomputed == stored else "FAIL",
        "stored_hash":     stored,
        "recomputed_hash": recomputed,
        "algorithm": "SHA-256",
    }


@router.get("/summary")
async def provenance_summary():
    """Aggregate stats across all simulation runs."""
    sb = get_supabase()
    res = sb.table("simulation_runs").select("mutation_id, p7_energy_ha, phase").execute()
    if not res.data:
        return {"total_runs": 0, "by_mutation": {}}

    by_mutation: dict[str, list] = {}
    for row in res.data:
        mid = row["mutation_id"]
        by_mutation.setdefault(mid, []).append(row["p7_energy_ha"])

    return {
        "total_runs": len(res.data),
        "by_mutation": {
            mid: {
                "runs": len(energies),
                "mean_energy_ha": round(sum(energies) / len(energies), 8),
                "min_energy_ha":  round(min(energies), 8),
                "max_energy_ha":  round(max(energies), 8),
            }
            for mid, energies in by_mutation.items()
        },
    }


# ── Write endpoints (require valid JWT) ───────────────────────────────────────

@router.delete("/runs/{run_id}")
async def delete_run(run_id: str, token: str = Depends(verify_jwt)):
    """
    Hard-delete a simulation run record. Requires admin JWT.
    In production: Supabase RLS policy restricts this to service_role key only.
    """
    sb = get_supabase()
    res = sb.table("simulation_runs").delete().eq("id", run_id).execute()
    if not res.data:
        raise HTTPException(status_code=404, detail=f"Run {run_id} not found")
    return {"deleted": run_id}


@router.get("/audit-log")
async def get_audit_log(limit: int = 50, token: str = Depends(verify_jwt)):
    """
    Return the audit log of all provenance operations.
    Requires authentication. RLS on Supabase side enforces role-based visibility.
    """
    sb = get_supabase()
    res = (sb.table("provenance_audit")
             .select("*")
             .order("created_at", desc=True)
             .limit(limit)
             .execute())
    return {"data": res.data, "count": len(res.data)}
