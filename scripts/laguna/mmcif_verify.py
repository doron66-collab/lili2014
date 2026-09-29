"""Structure-level verification that a stated point mutation is actually present
in a deposited coordinate file.

Rationale (SOLANGE Gate 2): a prior stabilizer-comparison result was retracted
because both the "wild-type" and "mutant" inputs carried the wild-type residue
at the mutated position. Internal consistency checks could not see it. The only
check that can is reading the coordinates themselves, which is what this does.

Evidence used, in order of authority:
  1. _atom_site      -- the residue actually present at the mapped position.
  2. _struct_ref_seq_dif -- the depositor's own record of differences from the
                        reference sequence ("ENGINEERED MUTATION", "VARIANT", ...).
  3. _struct_ref_seq -- UniProt <-> author numbering alignment, so a UniProt
                        position is mapped rather than assumed equal.
No filename, title, or upstream record is trusted.
"""
import re

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
AA1 = {v: k for k, v in AA3.items() if len(k) == 3 and k not in ("MSE", "SEC", "PYL")}


# ---------------------------------------------------------------- mmCIF reader
def _tokenize(line):
    """Split one mmCIF data line, honouring single/double quotes."""
    out, i, n = [], 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c in "'\"":
            j = i + 1
            while j < n:
                if line[j] == c and (j + 1 >= n or line[j + 1] in " \t"):
                    break
                j += 1
            out.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and line[j] not in " \t":
                j += 1
            out.append(line[i:j])
            i = j
    return out


def read_cif(path):
    """Return {category: {tag: [values]}}. Loops and key-value items unified:
    every tag maps to a list, length 1 for non-loop items."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        raw = fh.read().splitlines()

    # fold multi-line semicolon text fields into single tokens
    lines, k = [], 0
    while k < len(lines) or k < len(raw):
        if k >= len(raw):
            break
        ln = raw[k]
        if ln.startswith(";"):
            buf = [ln[1:]]
            k += 1
            while k < len(raw) and not raw[k].startswith(";"):
                buf.append(raw[k])
                k += 1
            k += 1
            lines.append(("TEXT", "\n".join(buf).strip()))
        else:
            lines.append(("LINE", ln))
            k += 1

    data, i = {}, 0
    while i < len(lines):
        kind, ln = lines[i]
        if kind == "TEXT":
            i += 1
            continue
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("data_"):
            i += 1
            continue

        if s.lower() == "loop_":
            i += 1
            tags = []
            while i < len(lines) and lines[i][0] == "LINE" and lines[i][1].strip().startswith("_"):
                tags.append(lines[i][1].strip().split()[0])
                i += 1
            cat = tags[0].split(".")[0]
            cols = {t.split(".", 1)[1]: [] for t in tags}
            names = [t.split(".", 1)[1] for t in tags]
            row = []
            while i < len(lines):
                kind2, ln2 = lines[i]
                if kind2 == "TEXT":
                    row.append(ln2)
                    i += 1
                else:
                    st = ln2.strip()
                    if not st or st.startswith("#"):
                        i += 1
                        if row:
                            continue
                        break
                    if st.startswith("_") or st.lower() in ("loop_",) or st.startswith("data_"):
                        break
                    row.extend(_tokenize(ln2))
                    i += 1
                while len(row) >= len(names):
                    for nm, val in zip(names, row[:len(names)]):
                        cols[nm].append(val)
                    row = row[len(names):]
            store = data.setdefault(cat, {})
            for nm in names:
                store.setdefault(nm, []).extend(cols[nm])
            continue

        if s.startswith("_"):
            parts = _tokenize(s)
            tag = parts[0]
            cat, nm = tag.split(".", 1) if "." in tag else (tag, "value")
            if len(parts) > 1:
                val = " ".join(parts[1:])
            else:
                val = lines[i + 1][1] if i + 1 < len(lines) else "?"
                i += 1
            data.setdefault(cat, {}).setdefault(nm, []).append(val)
            i += 1
            continue
        i += 1
    return data


def rows(data, cat):
    """Category as a list of dicts."""
    d = data.get(cat)
    if not d:
        return []
    n = max(len(v) for v in d.values())
    return [{k: (v[i] if i < len(v) else "?") for k, v in d.items()} for i in range(n)]


# ------------------------------------------------------------------- structure
def structure_summary(data):
    def one(cat, tag, default="?"):
        return (data.get(cat, {}).get(tag) or [default])[0]
    res = one("_refine", "ls_d_res_high")
    if res in ("?", ".", None):
        res = one("_em_3d_reconstruction", "resolution")
    return {
        "method": one("_exptl", "method"),
        "resolution_A": res,
        "title": one("_struct", "title"),
    }


def observed_residues(data):
    """{(auth_asym_id, auth_seq_id): comp_id} for polymer atoms, plus het list."""
    obs, het = {}, []
    for r in rows(data, "_atom_site"):
        comp = r.get("label_comp_id", "?")
        grp = r.get("group_PDB", "ATOM")
        ch = r.get("auth_asym_id") or r.get("label_asym_id")
        sq = r.get("auth_seq_id") or r.get("label_seq_id")
        if grp == "HETATM" and comp not in AA3:
            het.append((ch, sq, comp))
            continue
        obs[(ch, sq)] = comp
    return obs, het


def ref_alignments(data, accession):
    """Alignment blocks for the entity matching a UniProt accession."""
    ref_ids = {r["id"] for r in rows(data, "_struct_ref")
               if (r.get("pdbx_db_accession", "").upper() == accession.upper()
                   or r.get("db_code", "").upper().startswith(accession.upper()))}
    blocks = []
    for r in rows(data, "_struct_ref_seq"):
        if r.get("ref_id") not in ref_ids:
            continue
        try:
            blocks.append({
                "chain": r.get("pdbx_strand_id"),
                "db_beg": int(r["db_align_beg"]),
                "db_end": int(r["db_align_end"]),
                "auth_beg": int(r["pdbx_auth_seq_align_beg"]),
                "auth_end": int(r["pdbx_auth_seq_align_end"]),
            })
        except (KeyError, ValueError):
            continue
    return blocks, ref_ids


def seq_difs(data, ref_ids):
    out = []
    for r in rows(data, "_struct_ref_seq_dif"):
        if ref_ids and r.get("_ref_id", r.get("ref_id")) not in ref_ids:
            continue
        out.append({
            "chain": r.get("pdbx_pdb_strand_id"),
            "auth_seq": r.get("pdbx_auth_seq_num"),
            "pdb_mon": r.get("mon_id"),
            "db_mon": r.get("db_mon_id"),
            "details": r.get("details", "?"),
        })
    return out


def map_uniprot_pos(blocks, pos):
    """UniProt position -> [(chain, auth_seq_id)] using alignment offsets."""
    hits = []
    for b in blocks:
        if b["db_beg"] <= pos <= b["db_end"]:
            hits.append((b["chain"], b["auth_beg"] + (pos - b["db_beg"])))
    return hits


def poly_seq_scheme(data):
    """(chain, seqres_index) -> (auth_seq_num, mon_id) over SEQRES, including
    residues present in the construct but not modelled."""
    out = {}
    for r in rows(data, "_pdbx_poly_seq_scheme"):
        ch = r.get("pdb_strand_id")
        try:
            sid = int(r.get("seq_id"))
        except (TypeError, ValueError):
            continue
        out[(ch, sid)] = (r.get("auth_seq_num"), r.get("mon_id"))
    return out


def all_difs(data, accession=None):
    """Every depositor-recorded difference from the reference sequence.

    This is where an engineered mutation is declared. Filtering is by
    accession only when the entry references UniProt; entries that cite
    GenBank (e.g. 2FLU's KEAP1 chain) are returned unfiltered.
    """
    out = []
    for r in rows(data, "_struct_ref_seq_dif"):
        acc = r.get("pdbx_seq_db_accession_code", "?")
        if accession and acc not in ("?", ".", None) and acc.upper() != accession.upper():
            if r.get("pdbx_seq_db_name", "").upper() == "UNP":
                continue
        out.append({
            "chain": r.get("pdbx_pdb_strand_id"),
            "auth_seq": r.get("pdbx_auth_seq_num"),
            "pdb_mon": r.get("mon_id"),
            "db_mon": r.get("db_mon_id"),
            "details": (r.get("details") or "?").strip(),
            "accession": acc,
        })
    return out


def sifts_blocks(sifts_json, pdb_id, accession):
    """SIFTS residue-level alignment blocks for one accession.

    Returns [{chain, unp_start, unp_end, seqres_start}] where seqres_start is
    the SEQRES index corresponding to unp_start. Authoritative because SIFTS
    re-aligns to the current canonical sequence, rather than trusting the
    numbering the depositor used at deposition time.
    """
    inner = (sifts_json.get(pdb_id.lower(), {}).get("UniProt", {}) or {})
    dat = inner.get(accession)
    if not dat:
        return []
    blocks = []
    for m in dat.get("mappings", []):
        try:
            blocks.append({
                "chain": m["chain_id"],
                "unp_start": int(m["unp_start"]),
                "unp_end": int(m["unp_end"]),
                "seqres_start": int(m["start"]["residue_number"]),
            })
        except (KeyError, TypeError, ValueError):
            continue
    return blocks


def construct_scope(data, sifts_json, pdb_id, accession, chain=None):
    """How much of the protein the entry actually contains, for this accession.

    Verifying that a mutation is present is not sufficient. An entry can carry
    the mutation and still be useless: the six "TP53 R175H" entries in the PDB
    (6VQO/6VR5/6VRM/6VRN/6W51/7RM4) are nine-residue peptides, HMTEVVRHC, bound
    in an MHC groove. HIS is genuinely there at the mapped position and the
    depositor records the substitution -- but there is no p53 fold, no DNA
    binding domain and no structural zinc, so no domain-level active space can
    be built from them. This function makes that visible.
    """
    blocks = sifts_blocks(sifts_json, pdb_id, accession)
    spans = sorted({(b["unp_start"], b["unp_end"]) for b in blocks})
    mapped = sum(hi - lo + 1 for lo, hi in spans)
    seqres = {}
    for r in rows(data, "_entity_poly"):
        ids = (r.get("pdbx_strand_id") or "").split(",")
        seq = (r.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
        for c in ids:
            if c.strip():
                seqres[c.strip()] = len(seq)
    return {"unp_spans": [list(s) for s in spans],
            "mapped_residues": mapped,
            "seqres_length_by_chain": seqres,
            "seqres_length_at_site": seqres.get(chain) if chain else None}


def verify_sifts(path, sifts_json, pdb_id, accession, wt1, pos, mut1,
                 require_span=None, min_construct_len=None):
    """Adjudicate presence of `mut1` at UniProt `pos`, mapping via SIFTS.

    `require_span=(lo, hi)`: the UniProt range the entry must cover for the
    result to be usable at domain level (e.g. (94, 292) for the TP53 DNA
    binding domain). `min_construct_len`: minimum mapped residue count.
    When either is given and not satisfied, a confirmed mutation is downgraded
    to MUTANT_CONFIRMED_FRAGMENT_ONLY -- the mutation is real, the context is
    not there.

    Three independent sources, reported separately so a disagreement is
    visible rather than averaged away:
      observed   -- the residue in _atom_site at the mapped author number
      difs       -- the depositor's own engineered-mutation records
      mapping    -- SIFTS UniProt <-> SEQRES <-> author numbering
    """
    data = read_cif(path)
    res = {**structure_summary(data), "pdb": pdb_id}
    obs, het = observed_residues(data)
    scheme = poly_seq_scheme(data)
    difs = all_difs(data, accession)
    blocks = sifts_blocks(sifts_json, pdb_id, accession)
    res["hetero"] = sorted({h[2] for h in het if h[2] != "HOH"})
    res["waters"] = any(h[2] == "HOH" for h in het)
    res["unp_coverage"] = sorted({(b["unp_start"], b["unp_end"]) for b in blocks})
    res["all_engineered"] = [d for d in difs if "engineered" in d["details"].lower()]

    if not blocks:
        res["verdict"] = "NO_MAPPING"
        res["detail"] = f"SIFTS reports no mapping from {pdb_id} to {accession}."
        return res

    mapped = {}
    for b in blocks:
        if b["unp_start"] <= pos <= b["unp_end"]:
            sid = b["seqres_start"] + (pos - b["unp_start"])
            auth, seqres_mon = scheme.get((b["chain"], sid), (None, None))
            mapped[b["chain"]] = {"seqres_index": sid, "auth_seq": auth,
                                  "seqres_residue": seqres_mon}
    if not mapped:
        res["verdict"] = "SITE_OUTSIDE_CONSTRUCT"
        res["detail"] = (f"UniProt {pos} is outside the SIFTS-mapped coverage "
                         f"{res['unp_coverage']} of {accession} in {pdb_id}.")
        return res

    verdicts, seen = set(), {}
    for ch, m in mapped.items():
        comp = obs.get((ch, str(m["auth_seq"]))) if m["auth_seq"] else None
        # SEQRES tells us what the construct contains even when unmodelled
        construct = m["seqres_residue"]
        seen[f"{ch}/{m['auth_seq']}"] = {"modelled": comp, "seqres": construct}
        ref = comp or construct
        if comp is None and construct is None:
            verdicts.add("SITE_NOT_IN_SCHEME")
        elif comp is None:
            verdicts.add("SITE_NOT_MODELLED")
        elif AA3.get(ref) == mut1:
            verdicts.add("MUTANT_CONFIRMED")
        elif AA3.get(ref) == wt1:
            verdicts.add("WILD_TYPE_AT_SITE")
        else:
            verdicts.add("THIRD_RESIDUE")
    res["observed"] = seen
    res["difs_at_site"] = [d for d in difs
                           if d["auth_seq"] in {str(m["auth_seq"]) for m in mapped.values()}]
    for v in ("MUTANT_CONFIRMED", "THIRD_RESIDUE", "WILD_TYPE_AT_SITE",
              "SITE_NOT_MODELLED", "SITE_NOT_IN_SCHEME"):
        if v in verdicts:
            res["verdict"] = v
            break
    if len(verdicts) > 1:
        res["verdict"] += "_MIXED"
    res["detail"] = (f"{accession} {wt1}{pos}{mut1} -> " +
                     "; ".join(f"chain {k}: modelled={v['modelled']}, seqres={v['seqres']}"
                               for k, v in seen.items()) +
                     f". Depositor records at site: {res['difs_at_site'] or 'none'}.")

    # construct scope: is the surrounding domain actually present?
    site_chain = next(iter(mapped), None)
    scope = construct_scope(data, sifts_json, pdb_id, accession, chain=site_chain)
    res["construct_scope"] = scope
    if require_span or min_construct_len:
        lo_hi = require_span or (None, None)
        covers = any(b["unp_start"] <= lo_hi[0] and b["unp_end"] >= lo_hi[1]
                     for b in blocks) if require_span else True
        long_enough = (scope["mapped_residues"] >= min_construct_len
                       if min_construct_len else True)
        res["covers_required_span"] = covers
        res["construct_long_enough"] = long_enough
        if res["verdict"].startswith("MUTANT_CONFIRMED") and not (covers and long_enough):
            res["verdict"] = "MUTANT_CONFIRMED_FRAGMENT_ONLY"
            res["detail"] += (
                f" SCOPE: the entry maps only {scope['mapped_residues']} residues of "
                f"{accession} ({scope['unp_spans']}); required span {require_span}, "
                f"minimum length {min_construct_len}. The mutation is genuinely "
                f"present but the surrounding domain is not, so no domain-level "
                f"active space can be built from this entry.")
    return res


def verify(path, accession, wt1, pos, mut1):
    """Adjudicate whether `mut1` is present at UniProt `pos`.

    Returns a dict with a verdict and the evidence it rests on.
    """
    data = read_cif(path)
    summ = structure_summary(data)
    obs, het = observed_residues(data)
    blocks, ref_ids = ref_alignments(data, accession)
    difs = seq_difs(data, ref_ids)

    res = {**summ, "accession_matched": bool(blocks),
           "hetero": sorted({h[2] for h in het}),
           "chains_aligned": sorted({b["chain"] for b in blocks}),
           "construct_ranges": sorted({(b["db_beg"], b["db_end"]) for b in blocks}),
           "depositor_difs_at_site": [], "observed": {}, "auth_numbering": {}}

    if not blocks:
        accs = sorted({r.get("pdbx_db_accession", "?") for r in rows(data, "_struct_ref")})
        res["verdict"] = "WRONG_PROTEIN"
        res["detail"] = (f"No _struct_ref entity matches {accession}; entry references {accs}. "
                         "The cited structure is not this gene product.")
        return res

    mapped = map_uniprot_pos(blocks, pos)
    if not mapped:
        res["verdict"] = "SITE_OUTSIDE_CONSTRUCT"
        res["detail"] = (f"UniProt position {pos} lies outside every aligned block "
                         f"{res['construct_ranges']}; the crystallised construct does not contain it.")
        return res

    seen, verdicts = {}, set()
    for ch, auth in mapped:
        comp = obs.get((ch, str(auth)))
        seen[f"{ch}{auth}"] = comp or "NOT_MODELLED"
        res["auth_numbering"][ch] = auth
        if comp is None:
            verdicts.add("SITE_NOT_MODELLED")
        elif AA3.get(comp) == mut1:
            verdicts.add("MUTANT_CONFIRMED")
        elif AA3.get(comp) == wt1:
            verdicts.add("WILD_TYPE_AT_SITE")
        else:
            verdicts.add("UNEXPECTED_RESIDUE")
    res["observed"] = seen
    res["depositor_difs_at_site"] = [
        d for d in difs if d["auth_seq"] in {str(a) for _, a in mapped}]

    for v in ("MUTANT_CONFIRMED", "WILD_TYPE_AT_SITE", "UNEXPECTED_RESIDUE",
              "SITE_NOT_MODELLED"):
        if v in verdicts:
            res["verdict"] = v
            break
    if len(verdicts) > 1:
        res["verdict"] += "_MIXED_CHAINS"

    obs_names = sorted({c for c in seen.values()})
    res["detail"] = (f"UniProt {wt1}{pos}{mut1} maps to author numbering "
                     f"{res['auth_numbering']}; residue(s) present: {obs_names}. "
                     f"Depositor difference records at that position: "
                     f"{res['depositor_difs_at_site'] or 'none'}.")
    return res
