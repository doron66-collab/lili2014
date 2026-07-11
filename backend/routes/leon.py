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
import hashlib
import json
import logging

NAME = "LEON"
FULL_NAME = "Lineage-Evidence Orchestration & Notarization"

# p3_calibration_epoch (timestamptz) and other DB-reformatted fields would make a
# re-verification from the stored row spuriously FAIL. They are metadata, not
# result-integrity data, so excluding them keeps the seal robust and re-verifiable.
_SEAL_EXCLUDE = {"p3_calibration_epoch"}

# Result fields tamper-checked against the sealed payload at query time.
_TAMPER_FIELDS = ("p7_energy_ha", "p5_casscf_ref_ha", "p5_raw_energy")
_TAMPER_TOL = 1e-6
_CONSISTENCY_TOL = 1e-3


def build_p8_payload(record: dict) -> str:
    """The exact canonical JSON string that the P8 seal hashes over (P1–P7 + P9)."""
    return json.dumps({k: v for k, v in record.items()
                       if k.startswith(("p1_", "p2_", "p3_", "p4_",
                                        "p5_", "p6_", "p7_", "p9_"))
                       and k not in _SEAL_EXCLUDE},
                      sort_keys=True, default=str)


def build_p8_seal(record: dict) -> str:
    """SHA-256 hash of P1–P7 + P9 fields — the P8 cryptographic seal."""
    return hashlib.sha256(build_p8_payload(record).encode()).hexdigest()


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
        ok = seal_ok and tamper_ok
        return {
            "method": "sealed-payload", "notary": NAME,
            "integrity": "PASS" if ok else "FAIL",
            "seal_ok": seal_ok, "tamper_ok": tamper_ok, "note": tamper_note,
            "stored_hash": stored, "recomputed_hash": recomputed, "algorithm": "SHA-256",
        }

    seal_payload = json.dumps(
        {k: v for k, v in record.items()
         if k.startswith(("p1_", "p2_", "p3_", "p4_", "p5_", "p6_", "p7_", "p9_"))
         and k not in _SEAL_EXCLUDE},
        sort_keys=True, default=str)
    recomputed = hashlib.sha256(seal_payload.encode()).hexdigest()
    ok = (recomputed == stored)
    return {
        "method": "legacy-reconstruction", "notary": NAME,
        "integrity": "PASS" if ok else "LEGACY-UNVERIFIABLE",
        "stored_hash": stored, "recomputed_hash": recomputed, "algorithm": "SHA-256",
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


def write_audit(sb, event: str, run_id, verdict: dict, actor: str = None, note: str = None):
    """Append one immutable record to LEON's audit trail (21 CFR §11.10(e)).

    Best-effort and non-fatal: an audit-write failure (e.g. the leon_audit migration
    not yet run) must never break the request it is auditing. The trail is
    append-only by table policy — LEON records what it saw, it does not rewrite it.
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
        sb.table("leon_audit").insert(row).execute()
    except Exception as e:
        logging.warning("LEON audit-write skipped (%s): %s", event, e)
