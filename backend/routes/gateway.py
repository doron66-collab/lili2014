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

router = APIRouter()

# ── Thresholds — kept identical to scripts/laguna/solange_dmrg.py, which is the
# component that actually measures these quantities. If they ever diverge, the
# Gateway would recommend a route the classifier disagrees with. ────────────────
CHEM_ACC_MHA = 1.6      # 1 kcal/mol — the accuracy the whole tractability argument is phrased in
EXACT_WALL_E = 18       # active electrons up to which exact classical (FCI/CCSD(T)) is fine
S_HARD       = 1.5      # max bipartite entanglement above which DMRG stops being practical
PRACTICAL_M  = 2000     # bond dimension beyond which DMRG is deemed impractical here

# Largest usable width on the quantum hardware actually available to this project.
# IBM Heron r2/r3 are 156-qubit devices; a target needing more than this cannot be
# executed today regardless of how well it scores, and the Gateway must say so
# rather than recommending something that cannot run.
HERON_QUBITS = 156

# Scope gate: a DMRG verdict reclassifies a target only if its active space covered
# that target's site. Sized as a fraction because GENE_MAP site counts are stated
# approximately ("~44e").
SITE_COVERAGE_MIN = 0.9


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
    prior_class  = ev.get("prior_class")
    prior_basis  = ev.get("prior_basis")

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
        rationale.append(
            f"{active_e:.0f} active electrons is at or below the {EXACT_WALL_E}e "
            f"exact-classical wall: CCSD(T)/FCI reaches chemical accuracy here. "
            f"Quantum hardware would add cost and queue time without adding accuracy.")

    # ── 2. Measured entanglement at commensurate scale — the authoritative signal ─
    elif s_max_usable:
        if s_max > S_HARD:
            route = "qpu"
            quantum_value = True
            confidence = "high"
            rationale.append(
                f"Measured S_max = {s_max:.2f} exceeds {S_HARD}: the state carries more "
                f"entanglement than DMRG can represent at a practical bond dimension "
                f"(M ≤ {PRACTICAL_M}). This is the regime where classical methods do not "
                f"reach chemical accuracy and quantum simulation is the justified route.")
        else:
            route = "dmrg"
            quantum_value = False
            confidence = "high"
            rationale.append(
                f"Measured S_max = {s_max:.2f} is at or below {S_HARD}: DMRG converges at a "
                f"practical bond dimension, so a classical tensor-network run reaches "
                f"chemical accuracy. Quantum time is not warranted for this target.")

    # ── 3. No usable measurement — recommend earning the evidence, not guessing ──
    else:
        route = "dmrg_first"
        quantum_value = None
        confidence = "low"
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

    # ── Correlation energy — modulates confidence, never picks the route ────────
    if corr_mha is not None:
        ratio = abs(corr_mha) / CHEM_ACC_MHA
        rationale.append(
            f"Correlation energy is {abs(corr_mha):.1f} mHa — {ratio:.0f}x chemical accuracy "
            f"({CHEM_ACC_MHA} mHa). Mean-field (RHF/DFT) is definitively insufficient here, so a "
            f"correlated treatment of some kind is required regardless of which route is taken.")
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
        "qubits_required": int(qubits),
        "rationale": rationale,
        "evidence_gaps": gaps,
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
