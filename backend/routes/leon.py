"""
LEON — Lineage-Evidence Orchestration & Notarization.

The single notarization authority of the SOLANGE platform: the guardian that
re-verifies every simulation result crossing back into the platform, and refuses
any record whose lineage does not check out. Named in memory of Leon.

Design principle (DP1, §06.iv) — *verify, don't trust*: within a governed
pipeline, trust in a result is derived from reproducible cryptographic proof,
never from the identity of the boundary a result arrived through. LEON is the
component that makes the P1–P9 provenance chain of §06.iii self-enforcing rather
than declarative.

This module is the ONE canonical home of the P8 seal definition. Both the
ingestion path (routes.simulate ·/hpc/submit) and the query-time re-verification
path (routes.provenance ·/runs/{id}/verify) delegate here, so the seal can never
drift between "sealed at ingestion" and "re-verified later".
"""
import base64
import hashlib
import json
import logging
import os
from datetime import datetime, timezone

try:
    from cryptography.hazmat.primitives.asymmetric.ed25519 import (
        Ed25519PrivateKey, Ed25519PublicKey)
    _CRYPTO_AVAILABLE = True
except ImportError:  # pragma: no cover - only in an environment missing the dep
    _CRYPTO_AVAILABLE = False

NAME = "LEON"
FULL_NAME = "Lineage-Evidence Orchestration & Notarization"

# Result fields tamper-checked against the sealed payload at query time.
_TAMPER_FIELDS = ("p7_energy_ha", "p5_casscf_ref_ha", "p5_raw_energy")
_TAMPER_TOL = 1e-6
_CONSISTENCY_TOL = 1e-3


def _canonical_calibration_epoch(value):
    """Canonicalize p3_calibration_epoch to integer milliseconds since the Unix
    epoch (UTC), as a string.

    CORRECTED 2026-09-29 (Claude Science security review): this field used to be
    EXCLUDED from the P8 hash entirely, on the reasoning that a timestamptz's
    textual representation can be reformatted by the Postgres/Supabase round-trip
    (precision, 'Z' vs '+00:00', trailing zeros) without changing the represented
    instant -- which would make a naive re-hash spuriously fail. But the seal's
    own §11.10(a) criterion is defined as "consistent re-derivation conditional on
    the same calibration epoch" -- so the one field that criterion is conditioned
    on was the one field the seal did not cover, and could be changed silently
    with no detection. Fixed properly instead of left excluded: normalize to a
    fixed-precision representation (integer ms since epoch) that is invariant to
    reformatting but still changes if the actual instant changes, then include it
    like every other P1-P9 field.
    """
    if value is None:
        return None
    if isinstance(value, datetime):
        dt = value
    elif isinstance(value, str):
        s = value.strip()
        if s.endswith("Z"):
            s = s[:-1] + "+00:00"
        dt = datetime.fromisoformat(s)
    else:
        return str(value)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return str(int(dt.astimezone(timezone.utc).timestamp() * 1000))


def build_p8_payload(record: dict) -> str:
    """The exact canonical JSON string that the P8 seal hashes over (P1–P7 + P9)."""
    fields = {k: v for k, v in record.items()
              if k.startswith(("p1_", "p2_", "p3_", "p4_",
                                "p5_", "p6_", "p7_", "p9_"))}
    if "p3_calibration_epoch" in fields:
        fields["p3_calibration_epoch"] = _canonical_calibration_epoch(
            fields["p3_calibration_epoch"])
    return json.dumps(fields, sort_keys=True, default=str)


def build_p8_seal(record: dict) -> str:
    """SHA-256 hash of P1–P7 + P9 fields — the P8 integrity digest.

    NAMING, corrected 2026-09-29 (Claude Science): this is an integrity digest,
    not a signature. A hash reveals accidental or third-party tampering, but
    provides no non-repudiation and does not bind the record to an actor --
    anyone who can write to the store can recompute a valid hash after tampering.
    See sign_p8_seal()/verify_p8_signature() below for the actor-bound guarantee;
    this function's name is kept for backward compatibility with existing callers
    and stored data, but nothing in this module or its docs calls it a signature.
    """
    return hashlib.sha256(build_p8_payload(record).encode()).hexdigest()


def _signing_key():
    """LEON's Ed25519 private key, from the LEON_SIGNING_PRIVATE_KEY env var
    (a base64-encoded 32-byte seed). Returns None if unset or invalid -- a
    deployment without this configured still gets a valid integrity digest via
    build_p8_seal(), just not a non-repudiable signature. Generate a keypair with
    scripts/laguna/leon_keygen.py; the private key never leaves the deployment
    environment, and the public key is safe to publish (dissertation, docs).
    """
    if not _CRYPTO_AVAILABLE:
        return None
    raw = os.environ.get("LEON_SIGNING_PRIVATE_KEY")
    if not raw:
        return None
    try:
        return Ed25519PrivateKey.from_private_bytes(base64.b64decode(raw))
    except Exception as e:
        logging.warning("LEON_SIGNING_PRIVATE_KEY set but invalid: %s", e)
        return None


def sign_p8_seal(seal_hex: str) -> str | None:
    """Sign the P8 integrity digest with LEON's private key (Ed25519).

    Returns a base64-encoded signature, or None if no signing key is configured
    (soft-fail, matching this module's existing non-fatal philosophy -- a record
    without a signature still carries a valid, checkable integrity digest; it is
    simply not yet non-repudiable). Signing the digest, not the raw payload,
    keeps the signed quantity fixed-size and identical to what's already stored
    and displayed as p8_hash.
    """
    key = _signing_key()
    if key is None:
        return None
    return base64.b64encode(key.sign(seal_hex.encode())).decode()


def verify_p8_signature(seal_hex: str, signature_b64: str, public_key_b64: str) -> bool:
    """Verify an Ed25519 signature over a P8 integrity digest against a known
    public key. Returns False (not an exception) on any malformed input, missing
    dependency, or signature mismatch -- callers treat "not verified" uniformly
    regardless of cause, per this module's verify-don't-trust standard."""
    if not _CRYPTO_AVAILABLE or not signature_b64 or not public_key_b64:
        return False
    try:
        pub = Ed25519PublicKey.from_public_bytes(base64.b64decode(public_key_b64))
        pub.verify(base64.b64decode(signature_b64), seal_hex.encode())
        return True
    except Exception:
        return False


def notarize(prov: dict, jw: dict) -> dict:
    """Notarize an incoming external run at ingestion.

    LEON recomputes the P8 seal over the submitted provenance and re-checks the
    internal physics-consistency invariant (core + active == reference). It does
    NOT trust the submitter: `seal_ok=False` or `consistency_ok=False` means the
    caller must reject the record (never store it).

    Returns a verdict dict; the HTTP layer decides status codes.
    """
    submitted = prov.get("p8_hash")
    recomputed = build_p8_seal(prov)
    seal_ok = (submitted == recomputed)
    signature = sign_p8_seal(recomputed)  # None if no signing key configured

    consistency_ok = None
    ecore, ecas, eact = jw.get("ecore"), jw.get("e_casscf"), jw.get("e_active_exact")
    if None not in (ecore, ecas, eact):
        consistency_ok = abs((ecore + eact) - ecas) < _CONSISTENCY_TOL

    ok = seal_ok and (consistency_ok is not False)
    return {
        "notary": NAME,
        "ok": ok,
        "seal_ok": seal_ok,
        "consistency_ok": consistency_ok,
        "submitted_hash": submitted,
        "recomputed_hash": recomputed,
        "signature": signature,
        "signed": signature is not None,
    }


def reverify(record: dict) -> dict:
    """Re-attest a stored record's integrity on demand (query time).

    Robust path: the exact hashed JSON was stored verbatim in p8_seal_payload
    (text survives the DB round-trip byte-identically, unlike floats/timestamps).
    Re-hash it, then tamper-check the key result fields against the live columns.
    Legacy path (records sealed before p8_seal_payload existed): reconstruction is
    unreliable across the round-trip, so it is reported as LEGACY-UNVERIFIABLE.
    """
    stored = record.get("p8_hash", "")
    payload = record.get("p8_seal_payload")

    signature = record.get("p8_signature")
    public_key = os.environ.get("LEON_PUBLIC_KEY")
    # None (not False) when unsigned or unconfigured -- absence of a signature is
    # not itself a failure (older records, or a deployment not yet configured with
    # LEON_SIGNING_PRIVATE_KEY), only a mismatched signature is.
    signature_ok = verify_p8_signature(stored, signature, public_key) if signature else None

    if payload:
        recomputed = hashlib.sha256(payload.encode()).hexdigest()
        seal_ok = (recomputed == stored)
        tamper_ok, tamper_note = True, None
        try:
            pj = json.loads(payload)
            for f in _TAMPER_FIELDS:
                pv, cv = pj.get(f), record.get(f)
                if pv is not None and cv is not None and abs(float(pv) - float(cv)) > _TAMPER_TOL:
                    tamper_ok = False
                    tamper_note = f"{f}: sealed {pv} != stored {cv}"
        except Exception:
            pass
        ok = seal_ok and tamper_ok and (signature_ok is not False)
        return {
            "method": "sealed-payload", "notary": NAME,
            "integrity": "PASS" if ok else "FAIL",
            "seal_ok": seal_ok, "tamper_ok": tamper_ok, "note": tamper_note,
            "stored_hash": stored, "recomputed_hash": recomputed, "algorithm": "SHA-256",
            "signature_ok": signature_ok,
        }

    recomputed = hashlib.sha256(build_p8_payload(record).encode()).hexdigest()
    ok = (recomputed == stored) and (signature_ok is not False)
    return {
        "method": "legacy-reconstruction", "notary": NAME,
        "integrity": "PASS" if ok else "LEGACY-UNVERIFIABLE",
        "stored_hash": stored, "recomputed_hash": recomputed, "algorithm": "SHA-256",
        "signature_ok": signature_ok,
        "note": None if ok else "Pre-payload record; re-run to get a robustly verifiable seal.",
    }


def build_generic_payload(record: dict, exclude=frozenset()) -> str:
    """Canonical JSON for a NON-P1–P9 record (e.g. a DMRG A/B/C classification) —
    the whole dict, sorted keys, minus any self-referential seal fields. Used where
    a result doesn't fit the provenance schema but still needs LEON's guarantee:
    sealed at the source, rejected at ingestion if the seal doesn't recompute."""
    return json.dumps({k: v for k, v in record.items() if k not in exclude},
                      sort_keys=True, default=str)


def build_generic_seal(record: dict, exclude=frozenset()) -> str:
    return hashlib.sha256(build_generic_payload(record, exclude).encode()).hexdigest()


def notarize_generic(record: dict, hash_field: str, exclude=frozenset()) -> dict:
    """Notarize an incoming record that isn't a P1-P9 provenance record. Same
    verify-don't-trust contract as notarize(): recompute the seal from every field
    except the hash field itself, compare to what was submitted, and let the caller
    reject on mismatch. No physics-consistency check here — that's schema-specific
    to CASSCF/JW runs; generic records are sealed on completeness+non-tampering only.
    """
    submitted = record.get(hash_field)
    recomputed = build_generic_seal(record, exclude=exclude | {hash_field})
    seal_ok = (submitted == recomputed)
    return {
        "notary": NAME,
        "ok": seal_ok,
        "seal_ok": seal_ok,
        "submitted_hash": submitted,
        "recomputed_hash": recomputed,
    }


# The exact field set write_audit() hashes into each chain link -- verify_chain()
# must select the SAME subset from a row read back from the DB (which also
# carries id/created_at/seq/chain_hash, none of which were part of the original
# hash input), or every recomputed link would mismatch by construction.
_AUDIT_CHAIN_FIELDS = ("event", "run_id", "integrity", "seal_ok", "consistency_ok",
                       "method", "stored_hash", "recomputed_hash", "actor", "note")


def _chain_hash(prev_chain_hash: str, row: dict) -> str:
    """One link of the audit hash-chain: sha256(prev_chain_hash + canonical(row)).

    Added 2026-09-29 (Claude Science security review): "append-only" was enforced
    only by Postgres row-level-security policy (no update/delete policy defined) --
    real, but configuration, not construction. Whoever holds the service-role key
    can alter that policy. Hash-chaining each row to the previous one is
    tamper-EVIDENT BY CONSTRUCTION instead: deleting or editing any row breaks the
    chain at every row after it, detectably, independent of what access policy is
    in force at the time someone checks. Same mechanism already cited in §3.6 as
    QCIVET's (Yeniaras & Karimov 2026) hash-chained audit trail -- adopted here
    rather than left as an unused "structural neighbor" citation.
    """
    fields = {k: row.get(k) for k in _AUDIT_CHAIN_FIELDS}
    payload = json.dumps(fields, sort_keys=True, default=str)
    return hashlib.sha256(f"{prev_chain_hash or ''}{payload}".encode()).hexdigest()


def _last_chain_hash(sb):
    try:
        res = (sb.table("leon_audit").select("chain_hash")
               .order("seq", desc=True).limit(1).execute())
        return res.data[0]["chain_hash"] if res.data else None
    except Exception as e:
        logging.warning("LEON audit chain: could not read previous link: %s", e)
        return None


def verify_chain(sb, limit: int = None) -> dict:
    """Walk the audit trail in sequence order and recompute every chain_hash link.

    Returns {"ok": bool, "rows_checked": int, "break_at_seq": int|None} -- a break
    means every row from that seq onward could have been altered after the fact.
    This is the check that makes chain_hash meaningful: a hash nobody ever
    re-verifies is not tamper-evidence, just an unused column.
    """
    q = sb.table("leon_audit").select("*").order("seq")
    if limit:
        q = q.limit(limit)
    rows = q.execute().data or []
    prev = None
    for r in rows:
        stored_chain = r.get("chain_hash")
        expected = _chain_hash(prev, r)
        if expected != stored_chain:
            return {"ok": False, "rows_checked": len(rows), "break_at_seq": r.get("seq")}
        prev = stored_chain
    return {"ok": True, "rows_checked": len(rows), "break_at_seq": None}


def write_audit(sb, event: str, run_id, verdict: dict, actor: str = None, note: str = None):
    """Append one immutable record to LEON's audit trail (21 CFR §11.10(e)).

    Best-effort and non-fatal: an audit-write failure (e.g. the leon_audit migration
    not yet run) must never break the request it is auditing. The trail is
    append-only by table policy, AND hash-chained (see _chain_hash) so tampering
    is also detectable by construction, not only blocked by policy.
    """
    if sb is None:
        return
    row = {
        "event": event,
        "run_id": str(run_id) if run_id is not None else None,
        "integrity": verdict.get("integrity")
                     or ("PASS" if verdict.get("ok") else "REJECTED"),
        "seal_ok": verdict.get("seal_ok"),
        "consistency_ok": verdict.get("consistency_ok"),
        "method": verdict.get("method", "ingestion"),
        "stored_hash": verdict.get("stored_hash") or verdict.get("submitted_hash"),
        "recomputed_hash": verdict.get("recomputed_hash"),
        "actor": actor,
        "note": note or verdict.get("note"),
    }
    try:
        prev_chain_hash = _last_chain_hash(sb)
        row["chain_hash"] = _chain_hash(prev_chain_hash, row)
        sb.table("leon_audit").insert(row).execute()
    except Exception as e:
        logging.warning("LEON audit-write skipped (%s): %s", event, e)
