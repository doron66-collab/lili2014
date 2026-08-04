"""
security_log.py — the append-only record of denied access attempts.

leon_audit (backend/routes/leon.py) answers "was this record's content real and
unaltered". This module answers the other half of DP1 — "was this caller allowed
to do what they just tried to do". A denial that goes unrecorded gives an
IT/Security reviewer nothing to look at, which was the actual gap: the platform
enforced two real authorization boundaries (executive accounts cannot dispatch
real HPC/DMRG/QPU work; only admins can reach the admin endpoints) but neither
rejection left any trace anywhere. The check existed; the evidence that it fired
did not.

Deliberately narrow. This logs REJECTIONS only, from the two enforcement points
that already exist in the codebase — not every successful call. A row for every
allowed request would be noise that drowns the one signal an IT/Security
reviewer actually needs: attempted overreach. If a third enforcement point is
added later, it should call log_denied() too rather than growing its own
logging shape.
"""
import logging
import uuid
from datetime import datetime, timezone


def log_denied(sb, *, event: str, uid: str | None, email: str | None,
               role: str | None, endpoint: str, detail: str = "") -> None:
    """Record one denied access attempt. Best-effort and non-fatal — mirrors
    leon.write_audit()'s posture exactly: a logging failure (e.g. the
    access_denied_log migration not yet run) must never turn an already-correct
    403 into a 500 for the caller who was correctly denied."""
    if not sb:
        return
    row = {
        "id": str(uuid.uuid4()),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "event": event,
        "uid": uid,
        "email": email,
        "role": role,
        "endpoint": endpoint,
        "detail": detail,
    }
    try:
        sb.table("access_denied_log").insert(row).execute()
    except Exception as e:
        logging.error("security_log.log_denied failed (non-fatal, denial still enforced): %s", e)
