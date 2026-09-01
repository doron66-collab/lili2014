"""
Gate 1 — structural resolvability, live in SOLANGE.

Three gates decide whether a target's computed evidence means anything (see
routes/gate2.py's own module docstring for the full three-gate framing):
  Gate 1  Is there a residue to anchor an active space to at all? (this module)
  Gate 2  Which magnitude does the target's THERAPEUTIC mechanism actually
          require, and who decided that? (routes/gate2.py)
  Gate 3  Which computational method is trustworthy at that magnitude? — this
          is not a separate check: it IS the DMRG/SHCI classifier's own A/B/C
          verdict (Rung 3), not a gate that stands apart from a rung the way
          Gate 1 and Gate 2 do.

SINGLE SOURCE OF TRUTH, same rule as everywhere else in this codebase:
targets.json's structure_caveat field is authoritative for every mutation it
already covers (the four core NSCLC targets and their named point mutations).
This module does NOT duplicate that data into a second table — it reads
targets.json first, exactly like backend/routes/simulate.py already does.

WHAT THIS MODULE ACTUALLY OWNS: gate1_checks, a cache for mutations targets.json
has never seen — an NGS report can surface any gene/mutation, and re-deciding
"is this residue even resolvable" from scratch every time a name resurfaces
(this session's own PDB lookups, done by hand, are the reason this cache
exists at all) is wasted work once it has been decided once.

END-OF-DAY PROMOTION: gate1_checks is a WORKING cache, not a permanent store.
scripts/laguna/promote_gate1_checks.py merges every accumulated row into
targets.json (creating or updating the matching mutation's structure_caveat)
and then empties gate1_checks, so targets.json stays the one place a
structural verdict lives long-term and the cache never silently drifts into
being a second source of truth. Run it manually, end of day, reviewing the
diff before committing — the same "verify, don't trust" posture applied to
every other write into targets.json this project makes.

STATUS: like Gate 2, this does not itself block dispatch anywhere yet — no
live endpoint currently calls check_target() before running a job. It exists
so the verdict is recorded and reusable, not so a checkbox somewhere reads it.
Automating the actual PDB lookup (fetch structure, check residue coverage)
is future work — POST /check today only records a verdict a human already
worked out, it does not compute one.
"""
import logging
from datetime import datetime, timezone

from fastapi import APIRouter, Body, HTTPException, Header

from routes.simulate import get_supabase, _uid_from_auth, _TARGETS

router = APIRouter()

_GATE1_CHECK_COLUMNS = frozenset({
    "target", "pdb_id", "chain", "resi", "resolved", "reason", "source",
    "checked_by", "updated_at",
})


def _from_targets_json(target: str) -> dict | None:
    """targets.json is authoritative where it speaks. Mutation-level entry first
    (per-mutation structure_caveat, e.g. TP53_C275F); gene-level fallback for
    the LOF-style keys whose caveat lives on the gene (e.g. CDKN2A_LOF's real
    entry is under genes.CDKN2A, not mutations.CDKN2A_LOF — see targets.json's
    own _class_granularity note on why some genes carry no per-mutation split)."""
    mutations = _TARGETS.get("mutations") or {}
    entry = mutations.get(target)
    if entry and entry.get("structure_caveat"):
        return {
            "target": target, "resolved": entry["structure_caveat"].startswith("VERIFIED"),
            "reason": entry["structure_caveat"], "source": "targets.json (mutations)",
            "pdb_id": entry.get("pdb"),
        }
    gene = target.split("_")[0]
    genes = _TARGETS.get("genes") or {}
    gentry = genes.get(gene)
    if gentry and gentry.get("structure_caveat"):
        return {
            "target": target, "resolved": gentry["structure_caveat"].startswith("VERIFIED"),
            "reason": gentry["structure_caveat"], "source": "targets.json (genes." + gene + ")",
            "pdb_id": None,
        }
    return None


@router.get("/check/{target}")
async def check_target(target: str):
    """Look up target's Gate 1 verdict: targets.json first (authoritative,
    permanent), then the gate1_checks cache (working, pending promotion).
    Absent from both = never checked at all — distinct from 'checked and
    resolved', and returned as such rather than defaulting to either."""
    known = _from_targets_json(target)
    if known:
        return {"record": known}
    sb = get_supabase()
    if not sb:
        return {"record": None, "db": "not_configured"}
    try:
        res = sb.table("gate1_checks").select("*").eq("target", target).execute()
        if not res.data:
            return {"record": None}
        row = res.data[0]
        row["source"] = row.get("source") or "gate1_checks (pending promotion)"
        return {"record": row}
    except Exception as e:
        return {"record": None, "error": str(e)}


@router.get("/list")
async def list_pending(limit: int = 200):
    """Every row currently sitting in the working cache, awaiting the next
    promote_gate1_checks.py run. Not the full Gate 1 picture — targets.json's
    entries never pass through here, since they never needed the cache."""
    sb = get_supabase()
    if not sb:
        return {"records": [], "db": "not_configured"}
    try:
        res = (sb.table("gate1_checks").select("*")
                 .order("updated_at", desc=True).limit(limit).execute())
        return {"records": res.data or []}
    except Exception as e:
        return {"records": [], "error": str(e)}


@router.post("/check")
async def record_check(payload: dict = Body(...), authorization: str | None = Header(None)):
    """Record a structural-resolvability verdict for a mutation targets.json
    does not yet cover. Requires auth, same reasoning as Gate 2's upsert: this
    IS the sign-off action, recording who decided and when. Refuses to write
    over a target targets.json already covers -- that file is edited only via
    promote_gate1_checks.py (or by hand, reviewed), never by this endpoint,
    so there is exactly one path that can change the permanent record."""
    uid = _uid_from_auth(authorization)
    target = (payload or {}).get("target")
    if not target:
        raise HTTPException(400, "missing target")
    if _from_targets_json(target):
        raise HTTPException(
            409, f"'{target}' is already recorded in targets.json — the single "
                 f"source of truth. Edit it there (and re-run "
                 f"scripts/laguna/verify_consistency.py), not through this cache.")
    if "resolved" not in payload or "reason" not in payload:
        raise HTTPException(400, "missing 'resolved' (bool) or 'reason' (why)")
    sb = get_supabase()
    if not sb:
        return {"saved": False, "db": "not_configured"}
    row = {k: v for k, v in payload.items() if k in _GATE1_CHECK_COLUMNS}
    row["target"] = target
    row["checked_by"] = uid
    row["updated_at"] = datetime.now(timezone.utc).isoformat()
    try:
        sb.table("gate1_checks").upsert(row, on_conflict="target").execute()
        logging.info("Gate 1 check recorded for %s (resolved=%s) by %s",
                     target, payload.get("resolved"), uid)
        return {"saved": True, "record": row}
    except Exception as e:
        logging.error("Gate 1 check save failed: %s", e)
        return {"saved": False, "error": str(e),
                "hint": "run the gate1_checks migration — see RUN_GUIDE.md"}


@router.delete("/check/{target}")
async def delete_check(target: str, authorization: str | None = Header(None)):
    """Remove one row from the working cache — called by
    scripts/laguna/promote_gate1_checks.py right after that row's verdict has
    been merged into targets.json, never before. Not a general-purpose delete:
    there is nothing else this cache is for once a row is promoted."""
    _uid_from_auth(authorization)
    sb = get_supabase()
    if not sb:
        return {"deleted": False, "db": "not_configured"}
    try:
        sb.table("gate1_checks").delete().eq("target", target).execute()
        return {"deleted": True}
    except Exception as e:
        return {"deleted": False, "error": str(e)}
