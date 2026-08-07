"""
Computational Decision Gateway — recommend WHERE a problem should run, and say why.

SOLANGE's execution ladder (Rungs 1-4) answers "run this here". This module answers
the question that comes before it: *should* this run there at all, and what is the
evidence? Its purpose is as much to withhold expensive resources as to grant them —
the most valuable recommendation it makes is often "do not buy quantum time yet".

WHAT IT DECIDES ON
The routing rule is the one the dissertation already argues for in §06.i: an
evidence-based classification, never a size heuristic. Three signals, in order of
authority:

  1. Active-electron count against the exact-classical wall (18e). Below it, exact
     classical methods reach chemical accuracy and no amount of quantum hardware
     improves the answer — only the bill.
  2. S_max, the maximum bipartite entanglement actually measured by a DMRG sweep.
     This is the load-bearing signal: it is measured evidence about the physics,
     not an assumption about the size. Above S_HARD, DMRG cannot capture the state
     at a practical bond dimension, and quantum simulation becomes the only route.
  3. Correlation energy against chemical accuracy — how far mean-field falls short.
     Modulates confidence; does not by itself pick a route.

WHEN EVIDENCE IS ABSENT the Gateway does NOT guess a route. It recommends running
the cheap classical classifier (DMRG) first, and says so explicitly. Committing
quantum time without an entanglement measurement at the target's own scale is
speculation, and the Gateway names it as such rather than dressing it as a
recommendation.

SCOPE GATE — the same rule the UI applies (solangeDmrgCoversSite): a DMRG verdict
is authoritative for a target only if the run's active space actually covered that
target's site. A CAS(6,4) validation run says nothing about a 44e site, and is
reported as an evidence GAP rather than silently treated as a measurement.

IP BOUNDARY (deployment)
This module is the protected core. Its inputs are scalar evidence and an opaque
target label — never a structure, a sequence, or a patient record. That is what
lets a customer run every byte of sensitive data inside their own perimeter and
still consult the Gateway: the hybrid on-prem model asks it "given these numbers,
where should this go?" and receives a route plus a rationale. The customer sees
the recommendation and the reasoning behind it; the scoring function that produced
them stays here. No black box for the scientist, no source disclosure for us.
"""
import logging
from fastapi import APIRouter, Body, HTTPException

from routes import evidence as ev_grade

router = APIRouter()

# ── Thresholds — kept identical to scripts/laguna/solange_dmrg.py, which is the
# component that actually measures these quantities. If they ever diverge, the
# Gateway would recommend a route the classifier disagrees with. ────────────────
CHEM_ACC_MHA = 1.6      # 1 kcal/mol — the accuracy the whole tractability argument is phrased in
# Exact diagonalisation builds every determinant in the active space, so its
# cost is combinatorial in the size of that space and carries no dependence on
# entanglement whatever. For CAS(n,n) at singlet spin the determinant count is
# C(n, n/2)^2:  16e -> 1.7e8,  18e -> 2.4e9,  20e -> 3.4e10,  22e -> 5.0e11.
# Eighteen is where routine exact diagonalisation stops being routine — a few
# billion determinants is reachable with effort, tens of billions is not.
#
# It is a practical convention, not a derived constant, and it is NOT "the
# classical wall": DMRG carries the classical frontier far past it (FeMoco at
# CAS(113,76)). What this boundary marks is where the classifier must START
# ASKING about entanglement. Below it the question is moot, because exact
# diagonalisation settles the space regardless of how entangled the state is.
# Above it that shortcut is gone and S_max decides.
EXACT_WALL_E = 18   # active electrons up to which exact diagonalisation is routine
# Max bipartite entanglement above which DMRG stops being practical AT
# PRACTICAL_M. Calibrated, not derived: N2 at equilibrium measures S_max = 0.42
# (weakly correlated) and stretched N2 measures 2.86 (strongly correlated), and
# 1.5 sits between them. A policy of this pipeline, not a constant of nature -
# it marks where this pipeline's committed bond dimension runs out.
#
# Deliberately NOT replaced by a value derived from M ~ C*e^S. Calibration
# (results/calibration/N2_CAS10-16_ccpvdz_2026-08-05.json) measured that
# relation directly and found a growth exponent of 0.31, not 1: the exponential
# form is an upper bound, following from a flat entanglement spectrum, and real
# spectra decay. Both its prefactor and its exponent depend on the chemistry and
# on the size of the active space, so a fit on N2 does not transfer to a
# 44-electron site.
S_HARD       = 1.5

# The bond dimension this pipeline is willing to spend. Not derived: a stated
# commitment, and the value that gives "practical" its meaning everywhere below.
# For scale, cost runs as M^3, so 2000 is 512x the M=250 the DMRG-SCF solver
# uses by default, and roughly what a single node reaches in hours; the largest
# published chemistry DMRG runs go to M=6000 on thousands of cores for weeks
# (FeMoco, CAS(113,76)). Calibration on an N2 stretch series at CAS(10,16)
# found the requirement in the 160-256 range, an order of magnitude under this
# ceiling - so it is generous for small systems and untested against large ones.
# Raising it does not make a target classically tractable; it changes what this
# pipeline is prepared to pay before it says so.
PRACTICAL_M  = 2000

# Largest usable width on the quantum hardware actually available to this project.
# IBM Heron r2/r3 are 156-qubit devices; a target needing more than this cannot be
# executed today regardless of how well it scores, and the Gateway must say so
# rather than recommending something that cannot run.
HERON_QUBITS = 156

# Scope gate: a DMRG verdict reclassifies a target only if its active space covered
# that target's site. Sized as a fraction because GENE_MAP site counts are stated
# approximately ("~44e").
SITE_COVERAGE_MIN = 0.9

# ── Provenance discipline ──────────────────────────────────────────────────────
# A number's VALUE and a number's STANDING are different things, and conflating
# them is how a demonstration becomes a claim. Every quantity this Gateway routes
# on is therefore graded by where it came from, and the recommendation's confidence
# is capped by the weakest input that actually decided the route — never merely by
# whether the inputs were present.
#
# This matters concretely and today: the platform's stored correlation energies are
# computed over CAS(2e,2o)/STO-3G *model compounds* (toluene standing in for a Phe
# sidechain, methanethiol for Cys). They are real computations and they are correct
# about what they measured — two electrons of a proxy molecule. They are not a
# characterisation of a 44-electron protein active site, and the Gateway must never
# let one be read as the other. Likewise, GENE_MAP's active-electron counts are
# curated literature estimates, not measurements, and the route can hinge on them.
CONFIDENCE_RANK = {"low": 0, "medium": 1, "high": 2}

PROVENANCE_CEILING = {
    "measured_at_site": "high",    # a real run whose active space covered this target
    "literature":       "medium",  # curated from published chemistry — an expectation
    "unstated":         "medium",  # caller did not say; do not assume it is gold
    "measured_proxy":   "low",     # a real computation, but of a stand-in system
    "heuristic":        "low",     # a size guess with no chemistry behind it
}

# Within this fraction of a threshold, plausible estimate error could flip the
# route, so the call is borderline regardless of how good the number's source is.
BORDERLINE_FRAC = 0.15


def _cap(confidence, ceiling):
    """Lower `confidence` to `ceiling` if it currently claims more than that."""
    return confidence if CONFIDENCE_RANK[confidence] <= CONFIDENCE_RANK[ceiling] else ceiling


def _borderline(value, threshold):
    """True when a value sits close enough to a threshold that the route is fragile."""
    try:
        return abs(float(value) - float(threshold)) <= BORDERLINE_FRAC * float(threshold)
    except (TypeError, ValueError):
        return False


def _covers_site(measured_e, site_e):
    """True when a measurement's active space is commensurate with the target site.
    Unknown scope is not second-guessed — only a demonstrably undersized run is
    held back from carrying authority."""
    try:
        m, s = float(measured_e), float(site_e)
    except (TypeError, ValueError):
        return True
    if m <= 0 or s <= 0:
        return True
    return m >= SITE_COVERAGE_MIN * s


ROUTES = {
    "exact_classical": {
        "label": "Exact classical (CCSD(T) / FCI)",
        "runs_on": "CPU — laptop or a single Laguna node",
        "rung": 1,
        "spends_quantum_time": False,
    },
    "dmrg": {
        "label": "DMRG / tensor-network",
        "runs_on": "Laguna largemem (USC CARC)",
        "rung": 3,
        "spends_quantum_time": False,
    },
    "dmrg_first": {
        "label": "DMRG classifier first — evidence before commitment",
        "runs_on": "Laguna largemem (USC CARC)",
        "rung": 3,
        "spends_quantum_time": False,
    },
    "qpu": {
        "label": "Quantum hardware",
        "runs_on": "IBM Heron (external QPU)",
        "rung": 4,
        "spends_quantum_time": True,
    },
}


def recommend(ev: dict) -> dict:
    """Pure function: scalar evidence in, route + rationale out. No I/O, no state,
    no data access — deliberately, so it can be reasoned about and tested in
    isolation, and so the IP boundary above is structural rather than a promise."""
    rationale = []
    gaps = []
    # (ceiling, reason) pairs — each one an upper bound the final confidence may not
    # exceed. Collected as the evidence is examined, applied once at the end, and
    # reported, so a caller can see WHY a recommendation is only as strong as it is.
    caps = []

    # ── Read the evidence, tolerating absence rather than assuming a default ────
    def num(key):
        v = ev.get(key)
        try:
            return None if v is None else float(v)
        except (TypeError, ValueError):
            return None

    active_e = num("active_electrons")
    if active_e is None:
        raise HTTPException(400, "active_electrons is required — it is the one piece "
                                 "of evidence with no sensible default")

    active_orb = num("active_orbitals")
    qubits = num("full_qubits")
    if qubits is None and active_orb is not None:
        qubits = 2 * active_orb          # Jordan-Wigner: one qubit per spin-orbital
    if qubits is None:
        qubits = 2 * active_e            # the platform's own convention in GENE_MAP

    s_max        = num("s_max")
    s_max_scope  = num("s_max_scope_e")
    corr_mha     = num("correlation_mha")
    corr_scope   = num("correlation_scope_e")
    corr_basis   = ev.get("correlation_basis_set")
    prior_class  = ev.get("prior_class")
    prior_basis  = ev.get("prior_basis")

    # How the active-electron count itself was arrived at. The 18e wall is a hard
    # route decision, and it is taken against THIS number — so if the number is a
    # curated estimate rather than a measurement, the recommendation inherits that.
    ae_basis = str(ev.get("active_electrons_basis") or "unstated").strip().lower()
    if ae_basis not in PROVENANCE_CEILING:
        ae_basis = "unstated"

    # ── Is the entanglement measurement, if any, actually about THIS target? ────
    s_max_usable = s_max is not None
    if s_max_usable and s_max_scope is not None and not _covers_site(s_max_scope, active_e):
        s_max_usable = False
        gaps.append(
            f"An entanglement measurement exists but covered only {s_max_scope:.0f}e "
            f"against this target's ~{active_e:.0f}e site. It characterises that fragment, "
            f"not this target, so it is not used to route.")

    # ── 1. The exact-classical wall — the cheapest possible correct answer ──────
    if active_e <= EXACT_WALL_E:
        route = "exact_classical"
        quantum_value = False
        confidence = "high"
        decisive = f"active-electron count ({ae_basis})"
        caps.append((PROVENANCE_CEILING[ae_basis],
                     f"the route turns on the active-electron count, which is {ae_basis}"))
        rationale.append(
            f"{active_e:.0f} active electrons is at or below the {EXACT_WALL_E}e "
            f"exact-classical wall: CCSD(T)/FCI reaches chemical accuracy here. "
            f"Quantum hardware would add cost and queue time without adding accuracy.")
        if _borderline(active_e, EXACT_WALL_E):
            caps.append(("medium", f"{active_e:.0f}e sits close to the {EXACT_WALL_E}e wall"))
            rationale.append(
                f"This call is borderline: at {active_e:.0f}e the target is within "
                f"{int(BORDERLINE_FRAC*100)}% of the {EXACT_WALL_E}e wall, so a modest revision of "
                f"the active-space definition could move it off this route.")

    # ── 2. Measured entanglement at commensurate scale — the authoritative signal ─
    elif s_max_usable:
        decisive = "measured S_max at site scale"
        if s_max_scope is None:
            # Scope unstated: not second-guessed for routing (see _covers_site), but
            # it cannot then carry full confidence either.
            caps.append(("medium", "the S_max measurement's active space was not stated"))
        if s_max > S_HARD:
            route = "qpu"
            quantum_value = True
            # NOT "high": this is an absence claim (DMRG cannot represent the state),
            # relative to one classical method family (MPS/tensor-network). It does not
            # by itself rule out selected CI, FCIQMC, non-MPS tensor-network geometries,
            # or neural quantum states, whose sensitivity to entanglement differs. The
            # dissertation's own §06.i already states the classification is "conditional
            # on the current state of classical methods... and revisable as those
            # methods advance" — this asymmetry is that hedge, reflected in the number
            # a caller actually sees rather than only in prose. The s_max <= S_HARD
            # branch below stays "high": DMRG succeeding is a positive, demonstrated
            # result, not dependent on which other methods were or weren't tried.
            confidence = "medium"
            rationale.append(
                f"Measured S_max = {s_max:.2f} exceeds {S_HARD}: the state carries more "
                f"entanglement than DMRG can represent at a practical bond dimension "
                f"(M ≤ {PRACTICAL_M}). This is single-method evidence — DMRG specifically, "
                f"not classical computation in general — so confidence is capped at medium "
                f"pending a cross-check against a non-tensor-network classical method.")
        else:
            route = "dmrg"
            quantum_value = False
            confidence = "high"
            rationale.append(
                f"Measured S_max = {s_max:.2f} is at or below {S_HARD}: DMRG converges at a "
                f"practical bond dimension, so a classical tensor-network run reaches "
                f"chemical accuracy. Quantum time is not warranted for this target.")
        if _borderline(s_max, S_HARD):
            caps.append(("medium", f"S_max {s_max:.2f} sits close to the {S_HARD} threshold"))
            rationale.append(
                f"This call is borderline: S_max is within {int(BORDERLINE_FRAC*100)}% of the "
                f"{S_HARD} threshold, so it separates the classical and quantum routes only "
                f"narrowly. Treat the route as provisional and prefer re-measuring before "
                f"committing significant resource either way.")

    # ── 3. No usable measurement — recommend earning the evidence, not guessing ──
    else:
        route = "dmrg_first"
        quantum_value = None
        confidence = "low"
        decisive = "absence of any usable measurement"
        gaps.append("No entanglement measurement at this target's own scale.")
        rationale.append(
            f"{active_e:.0f}e is past the {EXACT_WALL_E}e exact-classical wall, but no "
            f"entanglement measurement covering this target exists yet. Routing to quantum "
            f"on size alone would be speculation: run the classical DMRG classifier first "
            f"(no quantum time) to establish whether classical methods actually fail here.")
        if prior_class:
            rationale.append(
                f"A prior Class {prior_class} label exists"
                + (f" ({prior_basis})" if prior_basis else "")
                + ", but a literature or heuristic label is an expectation, not a "
                  "measurement, and is not sufficient grounds to commit quantum time.")

    # ── Correlation energy — same scope discipline as S_max, for the same reason ─
    # A correlation energy computed over a CAS(2e,2o) model compound is a true
    # statement about two electrons of a stand-in molecule. Saying "mean-field is
    # insufficient HERE" on that basis silently promotes it into a claim about a
    # 44-electron protein site it never touched. Scope is checked, and the sentence
    # changes to match what was actually measured.
    if corr_mha is not None:
        ratio = abs(corr_mha) / CHEM_ACC_MHA
        corr_covers = _covers_site(corr_scope, active_e) if corr_scope is not None else None
        basis_note = f", {corr_basis}" if corr_basis else ""
        if corr_covers is False:
            rationale.append(
                f"Correlation energy of {abs(corr_mha):.1f} mHa ({ratio:.0f}x chemical accuracy) was "
                f"computed over a {corr_scope:.0f}e model space{basis_note} — not this target's "
                f"~{active_e:.0f}e site. It shows that model chemistry is correlated; it does NOT "
                f"establish that mean-field is insufficient for this target.")
            gaps.append(
                f"Correlation energy is proxy-scale ({corr_scope:.0f}e of a model compound vs a "
                f"~{active_e:.0f}e site) and is not evidence about this target.")
            caps.append(("low", "correlation energy is proxy-scale"))
        elif corr_covers is True:
            rationale.append(
                f"Correlation energy is {abs(corr_mha):.1f} mHa — {ratio:.0f}x chemical accuracy "
                f"({CHEM_ACC_MHA} mHa) — measured over this target's own active space{basis_note}. "
                f"Mean-field (RHF/DFT) is insufficient here, so a correlated treatment is required "
                f"regardless of which route is taken.")
        else:
            rationale.append(
                f"A correlation energy of {abs(corr_mha):.1f} mHa ({ratio:.0f}x chemical accuracy) was "
                f"supplied without stating the active space it was computed over{basis_note}, so it "
                f"cannot be attributed to this target. Treated as indicative only.")
            gaps.append("Correlation energy supplied without a scope — provenance unverifiable.")
            caps.append(("medium", "correlation energy scope unstated"))
    else:
        gaps.append("No correlation energy available — mean-field adequacy is unverified.")

    # ── Hardware feasibility overlay — an unrunnable recommendation is not one ──
    feasible = True
    if route == "qpu" and qubits > HERON_QUBITS:
        feasible = False
        confidence = "high"      # we are highly confident it CANNOT run, which is itself useful
        rationale.append(
            f"However: this target needs ~{qubits:.0f} qubits and the largest device available "
            f"is {HERON_QUBITS}q (IBM Heron). It is not executable on today's hardware — the "
            f"quantum route is correct in principle but blocked in practice.")

    # ── Apply the confidence ceilings ──────────────────────────────────────────
    # Done last and in one place: a recommendation may not present itself as more
    # certain than the weakest input that decided it. The reasons are returned
    # alongside, so "medium" is never an unexplained hedge.
    stated_confidence = confidence
    limits = []
    for ceiling, reason in caps:
        if CONFIDENCE_RANK[ceiling] < CONFIDENCE_RANK[confidence]:
            limits.append(f"capped to {ceiling} — {reason}")
        confidence = _cap(confidence, ceiling)

    # Always present, never empty: what a chemist must know before treating any of
    # this as a finding rather than a routing suggestion.
    limitations = list(limits)
    limitations.append(
        "This is a routing recommendation, not a chemical result. It states which "
        "computational method the available evidence justifies — not what the answer is.")
    if ae_basis != "measured_at_site":
        limitations.append(
            f"The active-space size ({active_e:.0f}e) is {ae_basis}, not a measurement of this "
            f"target's site; every threshold test above was taken against that number.")

    # ── Evidence grade — what may be SAID about any of this ────────────────────
    # Separate axis from the route. The route answers "where should this run"; the
    # grade answers "what standing does the answer have". A recommendation can be
    # perfectly sound as routing and still rest on demonstration-grade inputs, and
    # the reader must be told so in the same breath, not in a footnote.
    grade = ev_grade.grade_evidence({
        "target_ref":        ev.get("target_ref"),
        "active_electrons":  active_e,
        "computed_scope_e":  s_max_scope if s_max_scope is not None else corr_scope,
        "basis_set":         corr_basis or ev.get("basis_set"),
        "model_compound":    ev.get("model_compound"),
        "geometry_source":   ev.get("geometry_source"),
        "reference_source":  ev.get("reference_source"),
        "reference_delta_mha":     ev.get("reference_delta_mha"),
        "reference_tolerance_mha": ev.get("reference_tolerance_mha"),
    })
    if not ev_grade.may_state_chemical_finding(grade["evidence_grade"]):
        limitations.insert(0, grade["disclaimer"])

    info = ROUTES[route]
    return {
        "target_ref": ev.get("target_ref"),
        "recommended_route": route,
        "route_label": info["label"],
        "runs_on": info["runs_on"],
        "ladder_rung": info["rung"],
        "spends_quantum_time": info["spends_quantum_time"],
        "executable_today": feasible,
        "quantum_adds_value": quantum_value,     # None = not yet determinable
        "confidence": confidence,
        # What the branch alone would have claimed, before provenance capped it.
        # Exposed so the gap between the two is visible rather than hidden.
        "confidence_before_provenance": stated_confidence,
        "decided_by_evidence": decisive,
        "qubits_required": int(qubits),
        "rationale": rationale,
        "evidence_gaps": gaps,
        "limitations": limitations,
        "evidence": grade,
        "may_state_chemical_finding": ev_grade.may_state_chemical_finding(grade["evidence_grade"]),
        # Real quantum time is metered and paid for. Anything that spends it is
        # gated on a human saying yes — the Gateway recommends, it never commits.
        "requires_human_approval": info["spends_quantum_time"],
        "decided_by": "SOLANGE Computational Decision Gateway",
        "basis": f"§06.i evidence-based classification · thresholds: "
                 f"{EXACT_WALL_E}e wall, S_max {S_HARD}, {CHEM_ACC_MHA} mHa chemical accuracy",
    }


@router.post("/recommend")
async def recommend_route(payload: dict = Body(...)):
    """Recommend a computational route from evidence, with the reasoning attached.

    Accepts scalar evidence only — no structures, no sequences, no records — so a
    customer can call it from inside their own perimeter without exporting anything
    sensitive. Returns the route, whether it is executable today, whether quantum
    would add value, the confidence, and an explicit list of evidence gaps.

    Read-only and side-effect free: nothing is stored, nothing is dispatched. A
    recommendation is not an instruction — routes that spend real quantum time come
    back flagged requires_human_approval.
    """
    if not isinstance(payload, dict):
        raise HTTPException(400, "body must be a JSON object of evidence fields")
    out = recommend(payload)
    logging.info("Gateway recommended %s for %s (confidence=%s, quantum_value=%s)",
                 out["recommended_route"], out.get("target_ref"),
                 out["confidence"], out["quantum_adds_value"])
    return out


@router.get("/thresholds")
async def gateway_thresholds():
    """The decision thresholds, published deliberately.

    Publishing the thresholds is not the same as publishing the scoring function:
    a reviewer (or a regulator) can see exactly what line was drawn and check a
    verdict against it, which is what makes a recommendation auditable rather than
    oracular. This is the anti-black-box surface of the Gateway."""
    return {
        "exact_classical_wall_e": EXACT_WALL_E,
        "s_max_hard": S_HARD,
        "chemical_accuracy_mha": CHEM_ACC_MHA,
        "practical_bond_dimension": PRACTICAL_M,
        "hardware_qubit_ceiling": HERON_QUBITS,
        "site_coverage_min_fraction": SITE_COVERAGE_MIN,
        "note": "Thresholds match scripts/laguna/solange_dmrg.py, the component that "
                "measures these quantities, so a Gateway recommendation can never "
                "contradict the classifier's own verdict.",
    }
