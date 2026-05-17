"""
Admin endpoints — access log, user list, all simulation runs.
Requires service_role Supabase key on the backend; frontend must send
a valid Supabase JWT belonging to a user with role='admin' in users_profile.
"""
import os
from datetime import datetime, timezone

from fastapi import APIRouter, Header, HTTPException
from supabase import create_client

router = APIRouter()

SUPABASE_URL = os.environ.get("SUPABASE_URL", "")
SUPABASE_KEY = os.environ.get("SUPABASE_KEY", "")


def get_supabase():
    if not SUPABASE_URL or not SUPABASE_KEY:
        return None
    return create_client(SUPABASE_URL, SUPABASE_KEY)


def _require_admin(authorization: str | None, sb) -> str:
    """Verify the JWT belongs to a user with role='admin'. Returns user_id."""
    if not authorization or not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Authorization header required")
    import base64, json as _json
    try:
        token = authorization[7:]
        payload_b64 = token.split(".")[1]
        payload_b64 += "=" * (-len(payload_b64) % 4)
        payload = _json.loads(base64.urlsafe_b64decode(payload_b64))
        user_id = payload.get("sub")
        if not user_id:
            raise HTTPException(status_code=401, detail="Invalid token")
    except HTTPException:
        raise
    except Exception:
        raise HTTPException(status_code=401, detail="Invalid token")

    # Check users_profile role
    if sb:
        res = sb.table("users_profile").select("role").eq("id", user_id).single().execute()
        if not res.data or res.data.get("role") != "admin":
            raise HTTPException(status_code=403, detail="Admin role required")
    return user_id


@router.get("/runs")
async def admin_all_runs(limit: int = 50, authorization: str | None = Header(None)):
    """All simulation runs across all users, newest first."""
    sb = get_supabase()
    _require_admin(authorization, sb)
    res = (sb.table("simulation_runs")
             .select("id, created_at, user_id, mutation_id, mutation_name, p7_energy_ha, p8_hash, phase")
             .order("created_at", desc=True)
             .limit(limit)
             .execute())
    return {"data": res.data, "count": len(res.data)}


@router.get("/users")
async def admin_users(authorization: str | None = Header(None)):
    """All user profiles."""
    sb = get_supabase()
    _require_admin(authorization, sb)
    res = sb.table("users_profile").select("*").order("created_at", desc=True).execute()
    return {"data": res.data, "count": len(res.data)}


@router.get("/access-log")
async def admin_access_log(limit: int = 100, authorization: str | None = Header(None)):
    """Provenance audit log — simulation runs as audit trail."""
    sb = get_supabase()
    _require_admin(authorization, sb)
    res = (sb.table("simulation_runs")
             .select("id, created_at, user_id, mutation_name, mutation_id, p7_energy_ha, p8_hash, phase")
             .order("created_at", desc=True)
             .limit(limit)
             .execute())
    rows = res.data or []
    audit = [
        {
            "created_at": r["created_at"],
            "action": f"VQE simulation — {r['mutation_name']} ({r['mutation_id']})",
            "result": "success",
            "detail": f"Energy: {r['p7_energy_ha']} Hₐ · Phase: {r['phase']} · Seal: {(r['p8_hash'] or '')[:16]}…",
            "user_id": r["user_id"],
            "run_id": r["id"],
        }
        for r in rows
    ]
    return {"data": audit, "count": len(audit)}


@router.get("/stats")
async def admin_stats(authorization: str | None = Header(None)):
    """Summary stats for the admin dashboard."""
    sb = get_supabase()
    _require_admin(authorization, sb)
    runs  = sb.table("simulation_runs").select("id, mutation_id, created_at, user_id").execute()
    users = sb.table("users_profile").select("id, role").execute()
    data  = runs.data or []
    by_mutation = {}
    for r in data:
        m = r["mutation_id"]
        by_mutation[m] = by_mutation.get(m, 0) + 1
    return {
        "total_runs":   len(data),
        "total_users":  len(users.data or []),
        "by_mutation":  by_mutation,
        "last_run_at":  data[0]["created_at"] if data else None,
    }
