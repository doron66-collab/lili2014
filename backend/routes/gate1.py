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

AUTOMATED PDB LOOKUP (POST /auto_check): fetches the real deposited structure
from RCSB and checks whether the given chain/residue actually has ATOM
records at that position — the same objective fact a human was checking by
hand for every existing structure_caveat entry (e.g. "1U6D resolves residues
322-609; G333 sits shortly inside that boundary"). It answers the mechanical,
checkable half of Gate 1 ("does this file cover this residue") — never the
judgment half (is the flexible-linker/disordered-region call correct, is the
chosen chain the biologically relevant one). The verdict this returns is NOT
saved automatically: POST /check (the human-reviewed save) still requires a
person to look at the result and click save, the same "verify, don't trust"
posture as every other write into this project's records. This mirrors, on
purpose, the residue-coverage check computational chemists already do by eye
in a structure viewer — it does not replace judgment about WHICH chain/PDB is
biologically the right one to check in the first place, only automates
re-deriving whether a chosen one covers a chosen residue.
"""
import logging
import re
from datetime import datetime, timezone

import httpx
from fastapi import APIRouter, Body, HTTPException, Header

from routes.simulate import get_supabase, _uid_from_auth, _TARGETS
from routes.pdb import MUTATION_PDB_MAP

router = APIRouter()

_RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

# ── Standard amino acid single<->three-letter codes ───────────────────────────
_AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
_AA_1TO3 = {v: k for k, v in _AA_3TO1.items()}
_AA_3_PATTERN = "|".join(_AA_3TO1)


def _parse_variant_notation(mutation: str) -> tuple[str, int, str] | None:
    """Extract (wt_1letter, resnum, mut_1letter) from a single-residue
    substitution written in any of the common forms an NGS report actually
    uses: 'R882H', 'p.R882H', 'Arg882His', 'p.Arg882His', 'p.(Arg882His)'.
    Returns None for anything this cannot confidently parse as a simple
    substitution -- frameshift ('p.Gln1118*'), deletion, free text ('Loss of
    function') -- rather than guess. A wrong guess here would feed a wrong
    residue number into a real RCSB lookup and silently mislabel Gate 1.
    """
    if not mutation:
        return None
    m = mutation.strip().removeprefix("p.").strip("()")
    # Three-letter form: Arg882His
    m3 = re.fullmatch(rf"({_AA_3_PATTERN})(\d+)({_AA_3_PATTERN})", m, re.IGNORECASE)
    if m3:
        wt3, resnum, mut3 = m3.group(1).upper(), int(m3.group(2)), m3.group(3).upper()
        return _AA_3TO1[wt3], resnum, _AA_3TO1[mut3]
    # Single-letter form: R882H
    m1 = re.fullmatch(r"([A-Za-z])(\d+)([A-Za-z])", m)
    if m1 and m1.group(1).upper() in _AA_1TO3 and m1.group(3).upper() in _AA_1TO3:
        return m1.group(1).upper(), int(m1.group(2)), m1.group(3).upper()
    return None

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


def _parse_atom_site_chain_residues(cif_text: str) -> dict:
    """Minimal mmCIF _atom_site loop parser — same approach and same caveats as
    scripts/laguna/build_qm_cluster.py's parse_cif() (one _atom_site loop,
    values on one line each — true for every RCSB-deposited structure), but
    lighter: only chain/resseq/resname/record-type are needed here, not
    coordinates. Kept as a separate, self-contained copy rather than an
    import — this module runs on Render, build_qm_cluster.py runs on Laguna,
    different deployments of the same repo.

    Returns {chain: {resseq: resname}} for ATOM records only (HETATM residues
    are not part of the polymer chain a mutation's residue numbering refers
    to, and would give a false "resolved" for e.g. a nearby ion or water)."""
    lines = cif_text.splitlines()
    tag_start = None
    for i, raw in enumerate(lines):
        if raw.strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].strip().startswith("_atom_site."):
            tag_start = i + 1
            break
    if tag_start is None:
        raise ValueError("no _atom_site loop found — not a valid RCSB mmCIF file")
    tags, j = [], tag_start
    while j < len(lines) and lines[j].strip().startswith("_atom_site."):
        tags.append(lines[j].strip().split(".", 1)[1])
        j += 1
    idx = {t: n for n, t in enumerate(tags)}
    required = ["group_PDB", "auth_comp_id", "auth_asym_id", "auth_seq_id"]
    missing = [t for t in required if t not in idx]
    if missing:
        raise ValueError(f"_atom_site loop is missing expected columns: {missing}")

    by_chain: dict[str, dict[int, str]] = {}
    while j < len(lines):
        line = lines[j].strip()
        j += 1
        if not line or line.startswith("#") or line.startswith("_") or line == "loop_":
            break
        parts = line.split()
        if len(parts) != len(tags):
            continue
        if parts[idx["group_PDB"]] != "ATOM":
            continue
        chain = parts[idx["auth_asym_id"]]
        try:
            resseq = int(parts[idx["auth_seq_id"]])
        except ValueError:
            continue
        by_chain.setdefault(chain, {})[resseq] = parts[idx["auth_comp_id"]]
    if not by_chain:
        raise ValueError("parsed the _atom_site loop but found zero ATOM rows")
    return by_chain


async def _known_answer(target: str) -> dict | None:
    """Check BOTH permanent stores before anything else runs — targets.json,
    then the gate1_checks working cache (a colleague may have already
    resolved this exact target earlier today). A new employee has no way to
    know which mutations are already answered; this is what makes 'type the
    variant, get an answer' true regardless of whether that answer already
    exists or still needs deriving."""
    known = _from_targets_json(target)
    if known:
        return known
    sb = get_supabase()
    if not sb:
        return None
    try:
        res = sb.table("gate1_checks").select("*").eq("target", target).execute()
        if res.data:
            row = res.data[0]
            row["source"] = row.get("source") or "gate1_checks (pending promotion)"
            return row
    except Exception:
        pass
    return None


@router.post("/resolve_variant")
async def resolve_variant(payload: dict = Body(...)):
    """Turn what an NGS report actually gives you (gene + mutation notation,
    e.g. {"gene": "DNMT3A", "mutation": "p.Arg882His"}) into the fields Gate 1
    needs, WITHOUT requiring the end user to know a PDB ID, a chain letter, or
    a three-letter amino acid code -- or even whether this exact mutation was
    already checked by someone else. Checks BOTH permanent stores (targets.json,
    gate1_checks) FIRST, under every target-key form this project actually
    uses (a parsed substitution key, and the gene's own {GENE}_LOF key for
    anything this parser cannot reduce to one residue, e.g. a frameshift) --
    only derives a NEW answer if nothing already covers it. Read-only beyond
    that lookup: no RCSB call here (that's /auto_check, next) -- this only
    parses the notation and looks up a PDB ID already known to this project
    (targets.json / routes.pdb.MUTATION_PDB_MAP).

    Deliberately does NOT guess a chain: chain resolution needs the actual
    structure (see /auto_check's chain-omitted mode, which checks every chain
    present and reports which one(s) actually have the residue)."""
    gene = str((payload or {}).get("gene") or "").strip().upper()
    mutation = str((payload or {}).get("mutation") or "").strip()
    if not gene:
        raise HTTPException(400, "missing gene")
    parsed = _parse_variant_notation(mutation)

    if not parsed:
        # Can't reduce this to one residue (frameshift, deletion, free text) --
        # but this project's own convention already has a home for exactly
        # that case: {GENE}_LOF. Check it before giving up.
        lof_target = f"{gene}_LOF"
        known = await _known_answer(lof_target)
        if known:
            return {"resolved": known["resolved"], "target": lof_target,
                     "reason": known["reason"], "source": known.get("source"),
                     "already_known": True,
                     "note": f"'{mutation}' could not be reduced to a single residue, but "
                              f"{lof_target} is already answered — this project's convention "
                              f"for exactly this case (frameshift/deletion/generic LOF)."}
        return {"resolved": False,
                "error": f"could not parse '{mutation}' as a single-residue substitution "
                         f"(e.g. 'R882H' or 'p.Arg882His'), and {lof_target} has no existing "
                         f"answer either -- frameshift, deletion, and loss-of-function variants "
                         f"have no single residue to anchor Gate 1 to at all, and need a "
                         f"human's structural judgment, not this parser."}

    wt1, resnum, mut1 = parsed
    target = f"{gene}_{wt1}{resnum}{mut1}"

    known = await _known_answer(target)
    if known:
        return {"resolved": known["resolved"], "target": target,
                 "reason": known["reason"], "source": known.get("source"),
                 "already_known": True}

    # PDB ID: targets.json (this exact target, then a same-gene entry as a
    # weaker fallback) before the older routes.pdb literal map.
    pdb_id = None
    mutations = _TARGETS.get("mutations") or {}
    if target in mutations and mutations[target].get("pdb"):
        pdb_id = mutations[target]["pdb"]
    else:
        for key, entry in mutations.items():
            if key.startswith(gene + "_") and entry.get("pdb"):
                pdb_id = entry["pdb"]
                break
    if not pdb_id:
        pdb_id = MUTATION_PDB_MAP.get(target) or next(
            (v for k, v in MUTATION_PDB_MAP.items() if k.startswith(gene + "_")), None)

    if not pdb_id:
        return {"resolved": False, "target": target, "resi": resnum,
                "expect_resname": _AA_1TO3[wt1],
                "error": f"no PDB structure known for gene {gene} anywhere in this project "
                         f"(targets.json, routes.pdb) -- a chemist needs to supply one manually; "
                         f"this cannot be auto-derived."}
    return {
        "resolved": True, "target": target, "pdb_id": pdb_id,
        "resi": resnum, "expect_resname": _AA_1TO3[wt1],
    }


@router.post("/auto_check")
async def auto_check(payload: dict = Body(...)):
    """Automated half of Gate 1: does RCSB's own deposited structure for
    pdb_id actually have ATOM records at chain/resi? Read-only against RCSB,
    writes nothing — the caller reviews this result and calls POST /check
    (authenticated) to actually record it, same as every other Gate 1/2
    sign-off in this project.

    chain is OPTIONAL: when omitted (the normal case coming from
    /resolve_variant, which cannot know the chain without fetching the
    structure), every chain present in the file is checked and any chain(s)
    that actually have the residue are reported -- never assumed to be
    chain A, which is only sometimes the biologically relevant one."""
    pdb_id = str((payload or {}).get("pdb_id") or "").strip().upper()
    chain = str((payload or {}).get("chain") or "").strip() or None
    resi = (payload or {}).get("resi")
    expect_resname = str((payload or {}).get("expect_resname") or "").strip().upper() or None
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise HTTPException(400, "pdb_id must be a 4-character RCSB ID")
    try:
        resi = int(resi)
    except (TypeError, ValueError):
        raise HTTPException(400, "resi must be an integer residue number")

    url = _RCSB_CIF_URL.format(pdb_id=pdb_id)
    try:
        async with httpx.AsyncClient(timeout=30) as client:
            r = await client.get(url)
        if r.status_code != 200:
            return {"resolved": None, "error": f"RCSB returned HTTP {r.status_code} for {pdb_id} — check the PDB ID"}
        by_chain = _parse_atom_site_chain_residues(r.text)
    except httpx.HTTPError as e:
        return {"resolved": None, "error": f"could not reach RCSB: {e}"}
    except ValueError as e:
        return {"resolved": None, "error": str(e)}

    chains_to_check = [chain] if chain else sorted(by_chain)
    missing_chain = chain and chain not in by_chain
    if missing_chain:
        return {
            "resolved": False,
            "reason": f"chain '{chain}' not found in {pdb_id} — chains present: {sorted(by_chain)}",
            "pdb_id": pdb_id, "chain": chain, "resi": resi,
        }

    matches, near_misses = [], []
    for c in chains_to_check:
        residues = by_chain.get(c, {})
        if resi in residues:
            found = residues[resi]
            mismatch = expect_resname is not None and found != expect_resname
            matches.append({"chain": c, "found_resname": found, "resname_mismatch": mismatch})
        elif residues:
            near_misses.append({"chain": c, "range": [min(residues), max(residues)]})

    if not matches:
        ranges = "; ".join(f"chain {m['chain']}: {m['range'][0]}-{m['range'][1]}" for m in near_misses) or "no chain has any ATOM records"
        return {
            "resolved": False,
            "reason": (f"residue {resi} not found in {pdb_id} in {'chain ' + chain if chain else 'any chain'} — "
                       f"observed ranges: {ranges}. Falls in a gap or outside range (missing density / "
                       f"disordered / unresolved region, or simply out of range)."),
            "pdb_id": pdb_id, "chain": chain, "resi": resi,
        }
    clean = [m for m in matches if not m["resname_mismatch"]]
    chosen = clean[0] if clean else matches[0]
    reason = (f"VERIFIED — residue {resi} has ATOM records in {pdb_id} chain {chosen['chain']}"
              + (f" (also present in chain(s) {', '.join(m['chain'] for m in matches if m is not chosen)})" if len(matches) > 1 else ""))
    if chosen["resname_mismatch"]:
        reason += (f" — WARNING: expected residue name {expect_resname} but the structure has "
                   f"{chosen['found_resname']} at this position; check numbering convention/isoform before trusting this")
    return {
        "resolved": True, "reason": reason, "pdb_id": pdb_id, "chain": chosen["chain"], "resi": resi,
        "found_resname": chosen["found_resname"], "resname_mismatch": chosen["resname_mismatch"],
        "all_matching_chains": [m["chain"] for m in matches],
    }


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
