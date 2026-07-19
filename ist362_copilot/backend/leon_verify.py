"""
leon_verify.py — live seal verification, "LEON" in the copilot.

SOLANGE's provenance records carry a P8 SHA-256 seal computed over the P1-P7 and
P9 fields. LEON's principle is "verify, don't trust": integrity is derived from
recomputing the seal, not from trusting the submitter. This module reproduces
the platform's exact sealing algorithm so the copilot can re-attest any record
locally and detect tampering.

Reverse-engineered and confirmed bit-for-bit against a real HPC record:
  payload = json.dumps({all p1..p7 and p9 fields EXCEPT p3_calibration_epoch},
                       sort_keys=True, default=str)   # ensure_ascii=True (default)
  seal    = sha256(payload).hexdigest()

Editing any sealed field (e.g. an energy) changes the payload and therefore the
hash, so a mismatch pinpoints tampering — exactly what a reviewer wants.
"""
from __future__ import annotations

import hashlib
import json
import re

# The calibration epoch is a per-run timestamp, not scientific content, and the
# platform excludes it from the seal. Everything else in P1-P7,P9 is included.
_SEAL_EXCLUDE = {"p3_calibration_epoch"}
_SEAL_KEY_RE = re.compile(r"^p[1-7]_")


def _is_sealed_key(key: str) -> bool:
    return (_SEAL_KEY_RE.match(key) or key.startswith("p9_")) and key not in _SEAL_EXCLUDE


def canonical_payload(record: dict) -> str:
    """Rebuild the exact byte string the platform hashes for the P8 seal."""
    body = {k: v for k, v in record.items() if _is_sealed_key(k)}
    return json.dumps(body, sort_keys=True, default=str)


def recompute_seal(record: dict) -> dict:
    """
    Re-attest a provenance record's P8 seal. Returns a verdict dict:
      verifiable=False  -> record isn't in the flat P1-P9 schema with a p8_hash
      verified=True/False, with recomputed vs stored hash, and (if the sealed
      payload is present) the specific fields that were tampered with.
    """
    stored = record.get("p8_hash")
    has_flat = any(_SEAL_KEY_RE.match(k) for k in record)

    if not has_flat or not stored:
        nested = record.get("provenance")
        nested_seal = nested.get("p8_seal") if isinstance(nested, dict) else None
        return {
            "verifiable": False,
            "reason": ("No reconstructable P1-P7,P9 seal fields with a p8_hash were "
                       "found. Live verification supports the flat P1-P9 schema."),
            "stored_hash": stored or nested_seal,
        }

    payload = canonical_payload(record)
    recomputed = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    verified = recomputed == stored

    result = {
        "verifiable": True,
        "verified": verified,
        "notarized_by": "LEON",
        "algorithm": record.get(
            "p8_algorithm", "SHA-256 over P1-P7,P9 (sort_keys, default=str)"),
        "recomputed_hash": recomputed,
        "stored_hash": stored,
        "sealed_fields": sum(1 for k in record if _is_sealed_key(k)),
    }

    # If the sealed snapshot is stored, pinpoint which fields changed.
    sealed_payload = record.get("p8_seal_payload")
    if not verified and sealed_payload:
        try:
            sealed = json.loads(sealed_payload)
            changed = [k for k in sealed if record.get(k) != sealed.get(k)]
            missing = [k for k in sealed if k not in record]
            added = [k for k in record if _is_sealed_key(k) and k not in sealed]
            result["tampered_fields"] = sorted(set(changed) | set(missing) | set(added))
        except Exception:
            pass
    return result
