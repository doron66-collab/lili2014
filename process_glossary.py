#!/usr/bin/env python3
"""
Process dissertation_revised.html:
1. Add id="glos-SLUG" to each <tr> in the appendix-a glossary
2. Wrap acronym occurrences with <a href="#glos-SLUG" class="glos-link"> links
3. Add CSS for .glos-link
"""

import re
import sys
from collections import defaultdict

FILE_PATH = "/home/user/lili2014/public/dissertation_revised.html"

# ── Tier 1: link ALL occurrences ──────────────────────────────────────────────
TIER1 = [
    "SqDRIFT", "QPU", "VQE", "BQP", "QEC", "ZNE", "PEC", "CASSCF", "DMRG",
    "RCSB", "ADAPT-GQE", "EWF-SQD", "GMLP", "RBAC", "JWT", "FEDS", "NSCLC",
    "TP53", "KEAP1", "STK11", "CDKN2A", "KRAS", "EGFR", "ALK", "IND", "FCI",
    "FEP+", "NVQLink", "HIPAA", "CCSD(T)", "CUDA-Q", "OSF", "GARD", "Qiskit",
]

# ── Tier 2: link FIRST occurrence per <section> only ─────────────────────────
TIER2 = [
    "FDA", "HPC", "AI", "ML", "API", "DFT", "PDB", "NISQ", "CSA", "IRB", "QC", "DSR",
]

CSS_RULE = (
    ".glos-link { color: inherit; text-decoration: none; "
    "border-bottom: 1px dotted #999; } "
    ".glos-link:hover { border-bottom-color: currentColor; }"
)


def slug(acronym: str) -> str:
    """Convert acronym to HTML-safe ID slug."""
    return re.sub(r"[^A-Za-z0-9]", "-", acronym)


def make_pattern(acronym: str) -> re.Pattern:
    """
    Build a regex pattern for the given acronym.
    Handles special characters: (, ), +, -, etc.
    Uses \b where possible (before alphanumeric start, after alphanumeric end).
    """
    escaped = re.escape(acronym)
    # Left boundary: always \b before first char (all start with alnum)
    pat = r"\b" + escaped
    # Right boundary: if last char is alphanumeric → \b; else → (?![A-Za-z0-9])
    if acronym[-1].isalnum():
        pat += r"\b"
    else:
        pat += r"(?![A-Za-z0-9])"
    return re.compile(pat)


def build_replacement(acronym: str, glos_id: str, text: str, pattern: re.Pattern,
                      counts: dict) -> str:
    """Replace all occurrences in `text` and update counts."""
    def repl(m):
        counts[acronym] += 1
        return f'<a href="#{glos_id}" class="glos-link">{m.group(0)}</a>'
    return pattern.sub(repl, text)


def main():
    with open(FILE_PATH, "r", encoding="utf-8") as f:
        content = f.read()

    # ── Step 1: find section boundaries ──────────────────────────────────────
    # appendix-a section: from id="appendix-a" to next </section> at same depth
    # We'll use a simple approach: find the start tag and the matching </section>

    # Find appendix-a section start
    appendix_start_match = re.search(r'<section[^>]*id="appendix-a"[^>]*>', content)
    if not appendix_start_match:
        print("ERROR: Could not find appendix-a section", file=sys.stderr)
        sys.exit(1)
    appendix_start = appendix_start_match.start()

    # Find the matching </section> by counting nested sections
    def find_matching_section_end(text, start_of_open_tag):
        """Given position of opening <section...>, find the matching </section>."""
        depth = 0
        pos = start_of_open_tag
        while pos < len(text):
            open_m = re.search(r'<section\b', text[pos:])
            close_m = re.search(r'</section>', text[pos:])
            if close_m is None:
                break
            if open_m and open_m.start() < close_m.start():
                depth += 1
                pos += open_m.start() + len("<section")
            else:
                depth -= 1
                close_pos = pos + close_m.start()
                pos = close_pos + len("</section>")
                if depth == 0:
                    return pos  # position after </section>
        return len(text)

    appendix_end = find_matching_section_end(content, appendix_start)
    appendix_section = content[appendix_start:appendix_end]

    # Find terminology/reader's guide section
    terminology_start_match = re.search(r'<section[^>]*id="terminology"[^>]*>', content)
    if not terminology_start_match:
        print("ERROR: Could not find terminology section", file=sys.stderr)
        sys.exit(1)
    terminology_start = terminology_start_match.start()
    terminology_end = find_matching_section_end(content, terminology_start)

    print(f"Appendix-a section: chars {appendix_start}–{appendix_end}")
    print(f"Terminology section: chars {terminology_start}–{terminology_end}")

    # ── Step 2: add id="glos-SLUG" to each <tr> in appendix-a tbody ──────────
    # Extract the glossary acronyms as we go.
    # Pattern: <tr>\n  <td><strong>ACRONYM</strong></td>
    # We need to find the glossary acronyms from the table cells.

    # First, let's find all <tr> tags in the appendix-a tbody and add IDs.
    # The table structure: <tbody> ... <tr>\n<td><strong>ACRONYM</strong></td>...

    # Extract acronym from first <td> of each <tr>
    # Pattern in tbody: <tr>\s*<td><strong>ACRONYM</strong>
    tr_pattern = re.compile(
        r'(<tr>)(\s*<td><strong>)(.*?)(</strong></td>)',
        re.DOTALL
    )

    glossary_ids = {}  # acronym -> glos-id
    acronym_slugs = {}  # acronym -> slug (for known acronyms)

    def process_appendix_tr(m):
        acronym_raw = m.group(3).strip()
        glos_id = "glos-" + slug(acronym_raw)
        glossary_ids[acronym_raw] = glos_id
        return f'<tr id="{glos_id}">{m.group(2)}{acronym_raw}{m.group(4)}'

    new_appendix_section = tr_pattern.sub(process_appendix_tr, appendix_section)

    print(f"\nFound {len(glossary_ids)} glossary entries:")
    for acr, gid in sorted(glossary_ids.items()):
        print(f"  {acr!r} → #{gid}")

    # Update content with new appendix section
    content = content[:appendix_start] + new_appendix_section + content[appendix_end:]
    # Recalculate end positions since lengths may have changed
    # Re-find the section boundaries in updated content
    appendix_start_match2 = re.search(r'<section[^>]*id="appendix-a"[^>]*>', content)
    appendix_start2 = appendix_start_match2.start()
    appendix_end2 = find_matching_section_end(content, appendix_start2)

    terminology_start_match2 = re.search(r'<section[^>]*id="terminology"[^>]*>', content)
    terminology_start2 = terminology_start_match2.start()
    terminology_end2 = find_matching_section_end(content, terminology_start2)

    print(f"\nAfter ID insertion:")
    print(f"  Appendix-a: {appendix_start2}–{appendix_end2}")
    print(f"  Terminology: {terminology_start2}–{terminology_end2}")

    # ── Step 3: build patterns for all tier1+tier2 acronyms ──────────────────
    # Only keep acronyms that exist in the glossary (map to glos- IDs)
    # For acronyms not in glossary, we still need to link them if they appear

    def get_glos_id(acronym):
        """Get the glossary ID for an acronym, trying exact match first then case-insensitive."""
        if acronym in glossary_ids:
            return glossary_ids[acronym]
        # Try to find by matching the slug
        target = slug(acronym)
        for acr, gid in glossary_ids.items():
            if slug(acr) == target:
                return gid
        # Return a computed ID even if not in glossary (graceful fallback)
        return "glos-" + slug(acronym)

    tier1_info = [(acr, get_glos_id(acr), make_pattern(acr)) for acr in TIER1]
    tier2_info = [(acr, get_glos_id(acr), make_pattern(acr)) for acr in TIER2]

    print("\nTier 1 patterns (all occurrences):")
    for acr, gid, pat in tier1_info:
        in_glos = "✓" if acr in glossary_ids else "?"
        print(f"  {in_glos} {acr!r} → #{gid}  pattern={pat.pattern!r}")

    print("\nTier 2 patterns (first per section):")
    for acr, gid, pat in tier2_info:
        in_glos = "✓" if acr in glossary_ids else "?"
        print(f"  {in_glos} {acr!r} → #{gid}  pattern={pat.pattern!r}")

    # ── Step 4: split content into body + excluded zones ─────────────────────
    # We need to process only text outside:
    #   - appendix-a section
    #   - terminology section
    #   - <code>...</code>, <a ...>...</a>, <th>...</th>, <strong> inside glossary
    #
    # Strategy: split the full content into segments:
    #   "processable" vs "skip"
    # Then within processable segments, split on HTML tags and process text nodes.

    # Build exclusion zones (character ranges to skip)
    def find_excluded_zones(text):
        """Return list of (start, end) ranges that should NOT be processed."""
        zones = []

        # Excluded sections (appendix-a and terminology)
        zones.append((appendix_start2, appendix_end2))
        zones.append((terminology_start2, terminology_end2))

        # Find <code>...</code> tags (anywhere in document)
        for m in re.finditer(r'<code\b[^>]*>.*?</code>', text, re.DOTALL):
            zones.append((m.start(), m.end()))

        # Find existing <a ...>...</a> tags
        for m in re.finditer(r'<a\b[^>]*>.*?</a>', text, re.DOTALL):
            zones.append((m.start(), m.end()))

        # Sort and merge overlapping zones
        zones.sort()
        merged = []
        for start, end in zones:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append([start, end])
        return merged

    # ── Step 5: process the body text ────────────────────────────────────────
    # We'll use a tag-splitting approach on processable text segments.
    # For Tier 2, we track first occurrence per section.

    counts = defaultdict(int)

    # We process the content in chunks, skipping excluded zones.
    # For the processable parts, we further split by HTML tags and
    # only modify text nodes (content between > and <).

    # Tier 2 state: set of acronyms already linked in current section
    tier2_seen_in_section = set()

    def process_text_node(text_node: str, in_section: bool) -> str:
        """
        Process a raw text node (between HTML tags).
        Apply Tier 1 replacements (all) and Tier 2 replacements (first per section).
        """
        # Skip if empty or just whitespace
        if not text_node.strip():
            return text_node

        result = text_node

        # Tier 1: all occurrences
        for acr, gid, pat in tier1_info:
            def make_repl(a, g):
                def repl(m):
                    counts[a] += 1
                    return f'<a href="#{g}" class="glos-link">{m.group(0)}</a>'
                return repl
            result = pat.sub(make_repl(acr, gid), result)

        # Tier 2: first occurrence per section
        for acr, gid, pat in tier2_info:
            if acr not in tier2_seen_in_section:
                def make_repl2(a, g):
                    def repl(m):
                        tier2_seen_in_section.add(a)
                        counts[a] += 1
                        return f'<a href="#{g}" class="glos-link">{m.group(0)}</a>'
                    return repl
                new_result = pat.sub(make_repl2(acr, gid), result)
                # Only mark as seen if a replacement actually happened
                if new_result != result:
                    result = new_result

        return result

    # Split content into processable and excluded zones
    excluded_zones = find_excluded_zones(content)
    print(f"\nExcluded zones: {len(excluded_zones)}")

    # We'll process character by character through the content,
    # but more efficiently by working on spans between excluded zones.

    # Build list of processable spans
    processable_spans = []
    prev_end = 0
    for start, end in excluded_zones:
        if prev_end < start:
            processable_spans.append((prev_end, start, True))   # processable
        processable_spans.append((start, end, False))           # excluded
        prev_end = end
    if prev_end < len(content):
        processable_spans.append((prev_end, len(content), True))  # processable

    # Now process each processable span
    # Split each span into HTML tokens (tags vs text) and process text nodes
    # Track section boundaries for Tier 2 reset

    tag_pattern = re.compile(r'(<[^>]+>|<!--.*?-->)', re.DOTALL)

    result_parts = []

    for span_start, span_end, is_processable in processable_spans:
        chunk = content[span_start:span_end]

        if not is_processable:
            result_parts.append(chunk)
            continue

        # Split into tokens: tags and text nodes
        tokens = tag_pattern.split(chunk)
        processed_tokens = []

        for token in tokens:
            if not token:
                continue
            if token.startswith('<') and token.endswith('>'):
                # It's an HTML tag
                processed_tokens.append(token)
                # Check for section boundary (reset Tier 2 tracking)
                if re.match(r'<section\b', token):
                    tier2_seen_in_section.clear()
            else:
                # It's a text node - process it
                processed_tokens.append(process_text_node(token, True))

        result_parts.append("".join(processed_tokens))

    new_content = "".join(result_parts)

    # ── Step 6: add CSS rule ──────────────────────────────────────────────────
    # Add before first </style>
    css_injection = f"\n  {CSS_RULE}\n"
    new_content = new_content.replace("</style>", css_injection + "</style>", 1)

    # ── Step 7: write back ────────────────────────────────────────────────────
    with open(FILE_PATH, "w", encoding="utf-8") as f:
        f.write(new_content)

    print("\n=== REPLACEMENT SUMMARY ===")
    all_acr = TIER1 + TIER2
    total = 0
    for acr in all_acr:
        n = counts[acr]
        total += n
        tier = "T1" if acr in TIER1 else "T2"
        print(f"  [{tier}] {acr}: {n}")
    print(f"\nTotal replacements: {total}")
    print(f"File written to: {FILE_PATH}")


if __name__ == "__main__":
    main()
