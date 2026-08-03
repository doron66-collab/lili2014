"""
Evidence grading — what standing does a number have, and what may be said about it.

THE PROBLEM THIS EXISTS TO PREVENT
SOLANGE computes real quantities over stand-in systems: CAS(2e,2o)/STO-3G on
toluene where a Phe sidechain belongs, methanethiol where a Cys does. Those are
correct calculations of what they actually describe — two electrons of a model
compound — and they are exactly what a pipeline demonstration needs. They become
misconduct the moment a report lets them read as a characterisation of a protein
active site. The distance between "we computed this" and "this is true of the
target" is the whole distance between a demonstration and a claim.

This module is the single place that decides which side of that line a given set
of inputs falls on, so no report, table, export, or recommendation can quietly
answer it differently.

THE DISSERTATION'S POSITION, MADE ENFORCEABLE
§07 already states that operational stability and scientific validity are separate
axes: the platform's claim is that a governed pipeline runs end-to-end, NOT that a
computed energy is chemically correct. That has been a declaration. Here it becomes
a gate: grade_evidence() derives a grade from the inputs themselves and refuses to
accept a caller's assertion about its own standing. Verify, don't trust — the same
principle LEON applies to records, applied to claims.

This is an Information Systems contribution, and the grading framework IS the
artifact. The architecture, the provenance discipline, and this vocabulary are what
is being defended. The chemistry is not: no target in this platform is presently
validated to a standard that would support a chemical or biological finding, and
every grade below reference_exemplar says so in the output rather than in a footnote.
"""

# ── Grades, weakest first ──────────────────────────────────────────────────────
DEMONSTRATION_ONLY  = "demonstration_only"
TEMPLATE_CONFORMANT = "template_conformant"
REFERENCE_EXEMPLAR  = "reference_exemplar"

GRADE_RANK = {DEMONSTRATION_ONLY: 0, TEMPLATE_CONFORMANT: 1, REFERENCE_EXEMPLAR: 2}

GRADES = {
    DEMONSTRATION_ONLY: {
        "label": "Demonstration only — proof of concept",
        "means": "Values are placeholders or computed over stand-in systems. They "
                 "demonstrate that the pipeline runs end-to-end; they say nothing "
                 "about this target's chemistry.",
        "may_support": "Claims about the ARCHITECTURE: that a job was governed, "
                       "sealed, verified, and reproducible.",
        "may_not_support": "Any chemical, biological, or clinical statement about "
                           "the target — including whether it is druggable, whether "
                           "quantum simulation is warranted for it, or what its "
                           "energy is.",
    },
    TEMPLATE_CONFORMANT: {
        "label": "Template-conformant — measured at target scale, not externally validated",
        "means": "A computation covering this target's own active space exists and is "
                 "sealed, but has not been checked against an independent published "
                 "reference. It is internally consistent, not externally corroborated.",
        "may_support": "Architecture claims, plus internal comparisons across runs of "
                       "this platform.",
        "may_not_support": "Findings presented as chemically established. An "
                           "unvalidated calculation at the right scale is still an "
                           "unvalidated calculation.",
    },
    REFERENCE_EXEMPLAR: {
        "label": "Reference exemplar — reproduces a published, peer-reviewed result",
        "means": "The pipeline reproduces an independently published reference value "
                 "for this system, within a stated tolerance. The chemistry was "
                 "validated by its original authors and peer review; this platform "
                 "demonstrates it can recover it.",
        "may_support": "That the pipeline computes correctly for systems of this class "
                       "— the strongest claim this platform is designed to make.",
        "may_not_support": "Novel chemical findings. Reproducing a known answer "
                           "validates the instrument, not a discovery.",
    },
}

# Minimal basis sets: adequate for demonstrating a pipeline, never for a result that
# will be reported as chemistry. Listed explicitly so the judgement is auditable
# rather than buried in a heuristic.
MINIMAL_BASIS_SETS = {"sto-3g", "sto3g", "sto-6g", "minao", "min"}

# A calculation must cover at least this fraction of the target's own active space
# before it can speak for the target rather than for a fragment of it.
SITE_COVERAGE_MIN = 0.9

# The one sentence that must accompany any demonstration-grade output, wherever it
# appears — UI, export, API, or report. Defined once so it cannot drift between them.
POC_DISCLAIMER = (
    "PROOF OF CONCEPT — the values shown are not validated chemistry. They exist to "
    "demonstrate that the orchestration, provenance, and verification pipeline works "
    "end-to-end. Accurate values for this target have not yet been computed, and no "
    "chemical or biological conclusion may be drawn from what is shown."
)


def _f(v):
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def grade_evidence(ev: dict) -> dict:
    """Derive an evidence grade from the inputs themselves.

    Deliberately does NOT accept a declared grade: a caller asserting it is
    validated is precisely the claim that needs checking. Every promotion above
    demonstration_only must be earned by a specific, checkable property, and every
    reason a target failed to reach the next grade is returned so the path to
    upgrading it is explicit rather than mysterious.
    """
    blockers = []          # what stands between this target and the next grade up
    earned = []            # what it did establish

    site_e   = _f(ev.get("active_electrons"))
    scope_e  = _f(ev.get("computed_scope_e"))      # active electrons actually computed
    basis    = str(ev.get("basis_set") or "").strip().lower()
    model    = ev.get("model_compound")
    target   = ev.get("target_ref")
    ref_src  = ev.get("reference_source")          # DOI/citation for published values
    ref_delta_mha = _f(ev.get("reference_delta_mha"))
    ref_tol_mha   = _f(ev.get("reference_tolerance_mha")) or 1.6
    geom_src = ev.get("geometry_source")

    # ── Does the calculation describe the target, or a stand-in for it? ────────
    covers_site = None
    if site_e and scope_e:
        covers_site = scope_e >= SITE_COVERAGE_MIN * site_e
        if covers_site:
            earned.append(f"Computed over {scope_e:.0f}e, covering this target's ~{site_e:.0f}e site.")
        else:
            blockers.append(
                f"Computed over {scope_e:.0f}e against a ~{site_e:.0f}e site — it characterises "
                f"a fragment, not this target.")
    else:
        blockers.append("The computed active space was not stated, so it cannot be shown "
                        "to cover this target's site.")

    # A model compound is a stand-in by definition. Naming one is honest; it is also
    # disqualifying for any statement about the target it stands in for.
    if model and target and str(model).strip().lower() not in str(target).strip().lower():
        blockers.append(
            f"Computed on model compound '{model}', a stand-in — not on the target system itself.")

    if basis and basis in MINIMAL_BASIS_SETS:
        blockers.append(
            f"Basis set '{basis}' is minimal: adequate to demonstrate the pipeline, not to "
            f"support a reported chemical value.")
    elif basis:
        earned.append(f"Basis set '{basis}' is beyond minimal.")
    else:
        blockers.append("No basis set stated, so basis adequacy is unverifiable.")

    if geom_src:
        earned.append(f"Geometry provenance stated: {geom_src}.")
    else:
        blockers.append("No geometry source stated (experimental structure, PDB entry, "
                        "or documented construction).")

    # ── Promotion to reference_exemplar requires an EXTERNAL check ─────────────
    externally_validated = False
    if ref_src and ref_delta_mha is not None:
        if abs(ref_delta_mha) <= ref_tol_mha:
            externally_validated = True
            earned.append(
                f"Reproduces the published reference ({ref_src}) to {abs(ref_delta_mha):.2f} mHa, "
                f"within the {ref_tol_mha} mHa tolerance.")
        else:
            blockers.append(
                f"Differs from the published reference ({ref_src}) by {abs(ref_delta_mha):.2f} mHa, "
                f"outside the {ref_tol_mha} mHa tolerance.")
    else:
        blockers.append("No independent published reference value to check against.")

    # ── Decide ─────────────────────────────────────────────────────────────────
    if externally_validated and covers_site and not any(
            b.startswith(("Basis set", "Computed on model")) for b in blockers):
        grade = REFERENCE_EXEMPLAR
    elif covers_site:
        grade = TEMPLATE_CONFORMANT
    else:
        grade = DEMONSTRATION_ONLY

    info = GRADES[grade]
    out = {
        "evidence_grade": grade,
        "grade_label": info["label"],
        "what_this_means": info["means"],
        "may_support": info["may_support"],
        "may_not_support": info["may_not_support"],
        "established": earned,
        "blocking_upgrade": blockers,
    }
    # The disclaimer rides with anything not externally validated — the two lower
    # grades are both "not validated chemistry", they differ only in scale.
    if grade != REFERENCE_EXEMPLAR:
        out["disclaimer"] = POC_DISCLAIMER
    return out


def may_state_chemical_finding(grade: str) -> bool:
    """The gate a report must pass before presenting anything as a chemical result.
    Only an externally validated grade clears it, and even then only for the system
    that was validated."""
    return GRADE_RANK.get(grade, 0) >= GRADE_RANK[REFERENCE_EXEMPLAR]
