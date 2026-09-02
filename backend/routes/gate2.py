"""
Gate 2 — mechanism category & chemist sign-off, live in SOLANGE (not yet gating).

Three gates decide whether a target's computed evidence means anything:
  Gate 1  Is there a residue to anchor an active space to at all?
          (STRUCTURALLY_UNRESOLVED, frontend + backend)
  Gate 2  Which magnitude does the target's THERAPEUTIC mechanism actually
          require, and who decided that? (this module)
  Gate 3  Which computational method is trustworthy at that magnitude?
          (the DMRG/SHCI classifiers, §06.i)

Gate 2 answers a different question from Gate 3, not an earlier version of the
same one. Gate 3 (DMRG/SHCI) decides whether a target's GROUND-STATE energy is
classically tractable — that question, and the routine NGS-driven screening
that depends on it, needs no chemist in the loop and this module does not gate
it. Gate 2 only becomes load-bearing the moment someone wants a MECHANISM-
SPECIFIC quantity for actual drug design — a covalent activation barrier
(ΔG‡), a metal center's spin-state gap — which is a deliberate escalation a
researcher makes for a specific promising target, not a step every screened
mutation passes through.

STATUS: tracking and sign-off only. No live SOLANGE endpoint currently checks
gate3_allowed() before running anything — because the tool that would consume
it (a mechanism-specific energy calculator, e.g. a covalent ΔG‡ pipeline) does
not exist yet. Wiring this module's records into that tool's dispatch path,
once it exists, is what turns this from documentation into an actual gate.

THE MAPPING THIS MODULE OWNS (mechanism category -> what quantity matters,
what pair of states it is a difference between): a taxonomy problem, argued
out on four real cases (RhoA/covalent, SDHB/metal-spin, DNMT3A/interface,
TP53 stabilizer/local-comparison). It does not require a chemist to re-derive
per target once the target IS placed in a category.

THE JUDGMENT THIS MODULE CANNOT OWN, and does not pretend to (each field
below is a REQUIRED, externally-supplied value — there is no default):
  - which mechanism category a given target actually belongs in
  - which AO/orbital criterion is chemically right for a given site
  - oxidation state, spin state, and charge for a metal-containing site
  - whether a local calculation is even meaningful for the question at hand

Promoted from the standalone gate2_requirements.py (same taxonomy, same
Gate2Record shape) so the record persists in SOLANGE and is visible in the
Orchestration tab, instead of living only in a terminal session.
"""
import logging
from datetime import datetime, timezone

from fastapi import APIRouter, Body, HTTPException, Header

from routes.simulate import get_supabase, _uid_from_auth, _require_dispatch_allowed

router = APIRouter()

MECHANISM_CATEGORIES = {
    "covalent_reactive_cysteine": dict(
        question="Will a covalent warhead react with this residue, and how fast?",
        quantity="activation barrier (deltaG double-dagger) along the reaction "
                 "coordinate for nucleophilic attack",
        state_pair=("reactant complex", "transition state"),
        flags=["needs_ts_search", "needs_protonation_signoff"],
        example="RhoA Cys16 (or a stabilizing cysteine elsewhere on a "
                "conformational mutant — the reactive site need not be the "
                "mutated residue)",
    ),
    "metal_redox_center": dict(
        question="Does this metal-containing site have near-degenerate spin "
                 "states that classical single-reference methods cannot "
                 "resolve?",
        quantity="energy gap between candidate spin states of the metal "
                 "cluster",
        state_pair=("candidate spin state A", "candidate spin state B"),
        flags=["needs_metal_signoff", "needs_spin_signoff",
               "needs_oxidation_state_signoff"],
        example="SDHB C101Y — classify the wild-type Fe-S cluster site, not "
                "the (likely unassembled) mutant",
    ),
    "protein_interface_disruption": dict(
        question="Does the mutation change a protein-protein contact enough "
                 "to matter energetically, or is this a conformational/"
                 "oligomerization effect outside local electronic structure "
                 "entirely?",
        quantity="UNRESOLVED — may not be a local-electronic-structure "
                 "question at all; interface stability is usually MD/FEP "
                 "territory",
        state_pair=None,
        flags=["needs_gate2_scope_signoff"],
        example="DNMT3A R882 — included as a sample-expansion target "
                "precisely because Gate 2 is weak here, not strong",
    ),
    "structural_stabilizer_local_comparison": dict(
        question="Has the local electronic environment around the mutation "
                 "changed enough to matter for a stabilizer strategy?",
        quantity="local energy difference between wild-type and mutant "
                 "clusters at the same site, same method, same active space",
        state_pair=("wild-type cluster", "mutant cluster"),
        flags=["needs_protonation_signoff"],
        example="TP53 C275F stabilizer question — narrower and more "
                "defensible than a general druggability claim",
    ),
    "catalytic_loss_of_function": dict(
        question="Has catalysis been lost by disrupting local bonding at the "
                 "active site (as opposed to folding, expression, or an "
                 "unrelated interface)?",
        quantity="local energy/entanglement of the wild-type catalytic site",
        state_pair=("wild-type active site", None),
        flags=["needs_gate2_scope_signoff"],
        example="the eight archived Class B targets — this was the implicit "
                "category all along, never stated as one",
    ),
}

_GATE2_COLUMNS = frozenset({
    "target", "category", "metal_present", "spin_assigned_by",
    "oxidation_state_assigned_by", "protonation_assigned_by",
    "scope_signoff_by", "ts_search_configured", "source",
    "updated_by", "updated_at",
})


def _missing_requirements(rec: dict) -> list[str]:
    """Same logic as the standalone gate2_requirements.py's
    Gate2Record.missing_requirements(), operating on a plain dict so it can be
    applied to whatever a Supabase row hands back."""
    category = rec.get("category")
    if category not in MECHANISM_CATEGORIES:
        return [f"unknown category '{category}' — must be one of "
                f"{list(MECHANISM_CATEGORIES)}"]
    spec = MECHANISM_CATEGORIES[category]
    missing = []
    if not rec.get("source"):
        missing.append("no source recorded for WHY this target is in "
                        f"category '{category}' — who decided, and on "
                        "what basis")
    flags = spec["flags"]
    if "needs_metal_signoff" in flags and not rec.get("metal_present"):
        missing.append("category requires a metal-present determination; "
                        "none recorded")
    if "needs_spin_signoff" in flags and not rec.get("spin_assigned_by"):
        missing.append("spin state not assigned by anyone — this is "
                        "exactly the --dmrg-spin blocker: do not default "
                        "to spin=0 by omission")
    if "needs_oxidation_state_signoff" in flags and not rec.get("oxidation_state_assigned_by"):
        missing.append("oxidation state / cluster charge not assigned by "
                        "anyone")
    if "needs_protonation_signoff" in flags and not rec.get("protonation_assigned_by"):
        missing.append("protonation state not explicitly signed off — "
                        "automatic pdb2pqr/PROPKA output is not sufficient "
                        "for this category per protonate.py's own "
                        "disclaimer on metal-bound/reactive sites")
    if "needs_ts_search" in flags and not rec.get("ts_search_configured"):
        missing.append("transition-state search not configured — this "
                        "category is not a single-point calculation")
    if "needs_gate2_scope_signoff" in flags and not rec.get("scope_signoff_by"):
        missing.append("this category's own quantity is marked UNRESOLVED "
                        "or scope-limited in MECHANISM_CATEGORIES — an "
                        "explicit human decision is required to proceed at "
                        "all, not just to fill in a number")
    return missing


def _enrich(rec: dict) -> dict:
    missing = _missing_requirements(rec)
    return {**rec, "missing_requirements": missing, "gate3_allowed": len(missing) == 0}


@router.get("/categories")
async def list_categories():
    """The mechanism-category taxonomy itself — static, not per-target."""
    return {"categories": MECHANISM_CATEGORIES}


@router.get("/list")
async def list_records(limit: int = 100):
    """Every target's Gate 2 record, with missing_requirements/gate3_allowed
    computed fresh on every read (never stored — this must always reflect the
    CURRENT taxonomy, not whatever it evaluated to when last saved)."""
    sb = get_supabase()
    if not sb:
        return {"records": [], "db": "not_configured"}
    try:
        res = (sb.table("gate2_records").select("*")
                 .order("updated_at", desc=True).limit(limit).execute())
        return {"records": [_enrich(r) for r in (res.data or [])]}
    except Exception as e:
        return {"records": [], "error": str(e)}


@router.get("/record/{target}")
async def get_record(target: str):
    """One target's Gate 2 record. Absent = not yet categorised at all, which
    is itself meaningful (distinct from 'categorised but missing fields')."""
    sb = get_supabase()
    if not sb:
        return {"record": None, "db": "not_configured"}
    try:
        res = sb.table("gate2_records").select("*").eq("target", target).execute()
        if not res.data:
            return {"record": None}
        return {"record": _enrich(res.data[0])}
    except Exception as e:
        return {"record": None, "error": str(e)}


@router.post("/record")
async def upsert_record(payload: dict = Body(...), authorization: str | None = Header(None)):
    """Record (or update) a target's mechanism category and sign-offs.

    This is itself the sign-off action — it requires auth, same as any other
    action that spends or commits real judgment in SOLANGE. It does NOT
    validate that the signer is actually a qualified chemist (SOLANGE has no
    role for that today); it records WHO clicked save and WHEN, which is what
    'source' and 'updated_by' are for — verify, don't trust applies to the
    computation, not to identity verification this platform cannot perform.
    """
    uid = _uid_from_auth(authorization)
    target = (payload or {}).get("target")
    category = (payload or {}).get("category")
    if not target:
        raise HTTPException(400, "missing target")
    if category is not None and category not in MECHANISM_CATEGORIES:
        raise HTTPException(400, f"unknown category '{category}' — must be one of "
                                  f"{list(MECHANISM_CATEGORIES)}")
    sb = get_supabase()
    if not sb:
        return {"saved": False, "db": "not_configured"}
    row = {k: v for k, v in payload.items() if k in _GATE2_COLUMNS}
    row["target"] = target
    row["updated_by"] = uid
    row["updated_at"] = datetime.now(timezone.utc).isoformat()
    try:
        sb.table("gate2_records").upsert(row, on_conflict="target").execute()
        logging.info("Gate 2 record saved for %s (category=%s) by %s", target, category, uid)
        return {"saved": True, "record": _enrich(row)}
    except Exception as e:
        logging.error("Gate 2 record save failed: %s", e)
        return {"saved": False, "error": str(e),
                "hint": "run the gate2_records migration — see RUN_GUIDE.md"}


@router.post("/dispatch_custom_compound")
async def dispatch_custom_compound(payload: dict = Body(...),
                                    authorization: str | None = Header(None)):
    """The FIRST real gate: Gate 2 as a pre-dispatch check, not just tracking.

    Every other Gate 2 endpoint in this module only RECORDS a category and its
    sign-offs — nothing consumes gate3_allowed to stop anything (see the module
    docstring's STATUS line). This is the one exception, added because it was
    the concrete case that exposed the gap: a Gate-1-unresolvable, non-PDB,
    non-library target (e.g. an SDHB Fe-S-cluster proxy compound) had NO
    queue-dispatchable path at all before today — the only way to run one was
    to type a --geometry command by hand into a Laguna terminal, with no
    record of who decided the spin state, oxidation state, or AVAS criterion
    that made the run mean anything.

    This endpoint closes that gap the way Gate 2 was always meant to: it
    upserts the target's category + sign-off fields (same shape as /record,
    same table), computes missing_requirements against the CURRENT taxonomy,
    and refuses to queue anything if the list is non-empty. Only once every
    required field for that category is on record does this insert the actual
    hpc_dispatch row (job_type='dmrg', geometry mode — see
    solange_hpc.py's job_type=='dmrg' branch and simulate.py's dispatch_hpc).

    Uses _require_dispatch_allowed, the SAME executive-account block as every
    other path that spends real HPC/DMRG/QPU time — a custom-compound
    dispatch is not a lesser spend than any other.
    """
    sb = get_supabase()
    uid = _require_dispatch_allowed(authorization, sb)
    payload = payload or {}
    target = payload.get("target")
    category = payload.get("category")
    geometry = payload.get("geometry")
    avas = payload.get("avas")
    if not target:
        raise HTTPException(400, "missing target")
    if not geometry or not str(geometry).strip():
        raise HTTPException(400, "missing geometry (paste an .xyz file's contents)")
    if not avas:
        raise HTTPException(400, "missing avas (the AO criterion string, e.g. 'Fe 3d, S 3p')")
    if category not in MECHANISM_CATEGORIES:
        raise HTTPException(400, f"unknown category '{category}' — must be one of "
                                  f"{list(MECHANISM_CATEGORIES)}")
    if not sb:
        return {"queued": False, "db": "not_configured"}

    gate2_row = {k: v for k, v in payload.items() if k in _GATE2_COLUMNS}
    gate2_row["target"] = target
    gate2_row["category"] = category
    gate2_row["updated_by"] = uid
    gate2_row["updated_at"] = datetime.now(timezone.utc).isoformat()
    try:
        sb.table("gate2_records").upsert(gate2_row, on_conflict="target").execute()
    except Exception as e:
        return {"queued": False, "error": f"could not save the Gate 2 record: {e}",
                "hint": "run the gate2_records migration — see RUN_GUIDE.md"}

    missing = _missing_requirements(gate2_row)
    if missing:
        logging.info("Gate 2 blocked dispatch for %s (category=%s): %s",
                     target, category, missing)
        return {"queued": False, "blocked_by_gate2": True,
                "missing_requirements": missing,
                "note": "Fill in every field this category requires, then try again — "
                        "nothing was sent to the compute queue."}

    row = {
        "requested_by": uid, "status": "queued", "job_type": "dmrg",
        "key": target, "geometry": geometry, "avas": avas,
        "basis": payload.get("basis", "sto-3g"),
        "charge": int(payload.get("charge", 0)), "spin": int(payload.get("spin", 0)),
        "ncas": 0, "nelecas": 0,   # unused in geometry mode; AVAS decides the active space
        # bond_dims is NOT an hpc_dispatch column (checked against RUN_GUIDE.md's
        # migrations) — the agent falls back to its own default
        # ("250,500,1000,2000") when job.get("bond_dims") is absent, same as
        # every other geometry-mode job type. Do not add it here without a
        # migration first, or this insert breaks against a strict schema.
        "dmrg_scf": True, "dmrg_scf_maxm": int(payload.get("dmrg_scf_maxm", 250)),
    }
    try:
        res = sb.table("hpc_dispatch").insert(row).execute()
        did = (res.data or [{}])[0].get("id")
        logging.info("Gate 2 cleared dispatch for %s (category=%s) by %s — queued %s",
                     target, category, uid, did)
        return {"queued": True, "dispatch_id": did, "gate2_record": _enrich(gate2_row)}
    except Exception as e:
        return {"queued": False, "error": str(e),
                "hint": "run backend/migrations to create table hpc_dispatch"}
