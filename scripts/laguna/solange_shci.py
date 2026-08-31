#!/usr/bin/env python3
"""
solange_shci.py — independent SHCI A/B/C classifier for SOLANGE, optionally
cross-validated against a DMRG classification.

SHCI (Selected CI with perturbative correction, via the Dice implementation)
is a genuinely independent classical method — its failure mechanism is
determinant-sparsity growth, unrelated to DMRG's bond-dimension growth across
a one-dimensional chain (block2's MPS representation). Subordinating it to
"only ever checks a DMRG verdict" would waste exactly the capability that
makes it useful: SHCI can reach its own Class A verdict on a target DMRG was
never run against, not only confirm or contest an existing one.

CLASSIFICATION RULE (own convergence signal, own resource axis):
run SHCI at successively tighter epsilon (the heat-bath selection threshold —
looser epsilon keeps fewer determinants, tighter keeps more). If the energy
stops changing by more than chemical accuracy (1.6 mHa) at a PRACTICAL epsilon,
SHCI has delivered -> classically tractable (B). If it is still changing at the
tightest epsilon this run was willing to pay for, SHCI has NOT delivered ->
quantum-necessary (A), by the same "DMRG-specific evidence is not a
cross-method exclusion" logic solange_dmrg.py's own classify() states for
itself, now read the other way round for this solver.

HONEST SCOPE, stated once rather than buried: this is a ONE-SIGNAL classifier
(convergence only). DMRG's own classify() reads TWO signals — convergence AND
S_max (entanglement) — because a small system can be exactly DMRG-solvable yet
still reveal strong correlation through S_max alone. SHCI has no equivalent
"resource is about to blow up" signal wired in yet (the natural candidate is
the PT2 perturbative correction size, or the variational determinant count,
at the tightest epsilon — neither is read from Dice's output here, because
guessing the wrong parsing rule against output this session has not actually
seen would risk exactly the kind of confident-wrong result this project's own
evidence discipline exists to catch). Until that second signal is added and
verified against a real run, a Class B verdict from THIS script means only
"SHCI's own energy converged" — not "SHCI also confirms low correlation."

PRACTICAL_EPS below is a REASONED GUESS, not a derived or calibrated value —
stated with the same honesty already applied to DMRG's own PRACTICAL_M and
S_HARD (solange_dmrg.py, both flagged in their own source comments as
commitments pending calibration, not derivations). Expect it to be revised
once real chemist-supplied cases exist to calibrate against.

TWO INDEPENDENT USES, both real, neither privileged over the other:
  1. Standalone classification: omit --dmrg-classification-id. SHCI reaches
     its own A/B/C verdict on whatever active space you give it — including a
     target DMRG has never touched.
  2. Cross-validation: pass --dmrg-classification-id <uuid>. The backend
     ADDITIONALLY fetches that DMRG record's own stored energy and computes
     delta_mha/agreement against it (DP1, verify-don't-trust: never taken from
     this script's own claim). This does not replace SHCI's own classification
     above — both are recorded on the same submission.

USAGE
  # standalone — no DMRG record needed or referenced
  python3 solange_shci.py --geometry cluster.xyz --charge 1 --basis sto-3g \\
      --avas "N 2p, O 2p" --key SDHB_C101Y \\
      --dice-scripts ~/lili2014/Dice/scripts --sweep-eps 1e-2,1e-3,5e-4,1e-4 --submit

  # cross-validated against an existing DMRG classification
  python3 solange_shci.py --geometry cluster.xyz --charge 1 --basis sto-3g \\
      --avas "N 2p, O 2p" --key DNMT3A_R882 \\
      --dmrg-classification-id <uuid-from-dmrg-submit-response> \\
      --dice-scripts ~/lili2014/Dice/scripts --sweep-eps 1e-3,5e-4,1e-4 --submit
"""
import argparse
import hashlib
import json
import sys
import time
import uuid
from pathlib import Path

EXACT_WALL_E = 18       # same wall solange_dmrg.py uses -- method-independent:
                         # below it, exact diagonalisation already suffices.
CHEM_ACC_MHA = 1.6
# REASONED GUESS, not calibrated -- see module docstring. Mirrors PRACTICAL_M's
# epistemic status exactly: a stated commitment about what this pipeline is
# willing to pay for, not a derived property of SHCI or of chemistry.
PRACTICAL_EPS = 1e-4


def classify_shci(active_electrons, energies):
    """Map SHCI convergence behaviour to an A/B/C class with an explicit
    rationale. energies: [(epsilon, E), ...], loosest epsilon first, tightest
    last. ONE signal only -- see module docstring's Honest Scope section."""
    dE_final = abs(energies[-1][1] - energies[-2][1]) * 1000 if len(energies) >= 2 else None
    converged = (dE_final is not None) and (dE_final < CHEM_ACC_MHA) \
                and (energies[-1][0] <= PRACTICAL_EPS)
    dE_str = "n/a" if dE_final is None else f"{dE_final:.2f}"

    if active_electrons <= EXACT_WALL_E:
        return "C", (f"{active_electrons}e <= {EXACT_WALL_E}e exact-classical wall — "
                     f"CCSD(T)/FCI reaches chemical accuracy; classical sufficient.")
    if converged:
        return "B", (f"SHCI reaches chemical accuracy at practical epsilon="
                     f"{energies[-1][0]:g} (ΔE={dE_str} mHa) — classical (SHCI) delivers; "
                     f"quantum-advantaged, not necessary. CONVERGENCE-ONLY VERDICT: no "
                     f"determinant-count/PT2-correction hardness signal is read yet — see "
                     f"this script's own module docstring, Honest Scope.")
    return "A", (f"SHCI not at chemical accuracy by practical epsilon={energies[-1][0]:g} "
                f"(ΔE={dE_str} mHa) — quantum-necessary. CONVERGENCE-ONLY VERDICT: this is "
                f"SHCI-specific evidence, not a cross-method exclusion (the same reading "
                f"solange_dmrg.py's own classify() applies to a DMRG Class A) — a classical "
                f"method other than SHCI or DMRG could still succeed here.")


# ── LEON seal (self-contained — mirrors backend/routes/leon.py's generic seal
# bit-for-bit and solange_dmrg.py's own _seal_payload, deliberately duplicated
# rather than imported: this script must run standalone on Laguna with no
# dependency on the SOLANGE backend package). LEON re-verifies this at
# ingestion; a mismatch is REJECTED, not stored.
def _seal_payload(record, exclude):
    return json.dumps({k: v for k, v in record.items() if k not in exclude},
                      sort_keys=True, default=str)


def _submit_shci(api, out):
    """POST the sealed SHCI classification to SOLANGE. Best-effort: prints the
    outcome but never raises — a submit failure must not discard the local JSON
    already written to --out."""
    import urllib.request
    try:
        url = api.rstrip("/") + "/api/simulate/hpc/shci/submit"
        body = json.dumps(out, default=str).encode()
        req = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            resp = json.loads(r.read().decode())
        print(f"SUBMITTED → {url}")
        print(f"  status={resp.get('status')}  seal_ok={resp.get('seal_ok')}  "
              f"db={resp.get('db_status')}  run_id={resp.get('run_id')}")
        if resp.get("delta_mha") is not None:
            print(f"  delta_mha={resp.get('delta_mha')}  agreement={resp.get('agreement')}"
                  f"  (computed server-side against the referenced DMRG record)")
    except Exception as e:
        print(f"  SUBMIT FAILED (result is still safe locally in --out): {e}", file=sys.stderr)
        if hasattr(e, "read"):
            try:
                print(f"  server said: {e.read().decode()}", file=sys.stderr)
            except Exception:
                pass


def read_xyz(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    if lines and lines[0].strip().isdigit():
        lines = lines[2:]
    return "\n".join(lines)


def detect_hardware(n_threads):
    """Best-effort hostname/CPU tag — Dice/SHCI is CPU (+MPI) only, no GPU path."""
    import platform
    return f"CPU · {platform.node()} · {n_threads} threads"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geometry", required=True)
    ap.add_argument("--avas", required=True)
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--threshold", type=float, default=0.2)
    ap.add_argument("--max-memory", type=int, default=16000)
    ap.add_argument("--df-auxbasis", default="def2-universal-jkfit")
    ap.add_argument("--dice-scripts", required=True,
                    help="path to Dice/scripts (contains shci.py)")
    ap.add_argument("--sweep-eps", default="1e-2,1e-3,5e-4,1e-4",
                    help="comma-separated epsilon schedule, LOOSEST first, TIGHTEST last — "
                         "each value is run as its own independent SHCI solve, so the "
                         "energy at each is recorded for the convergence check (unlike a "
                         "single Dice call given the whole schedule, which returns only "
                         "the final energy)")
    ap.add_argument("--num-thrds", type=int, default=8)
    ap.add_argument("--key", required=True, help="target key, e.g. DNMT3A_R882")
    ap.add_argument("--dmrg-classification-id", default=None,
                    help="OPTIONAL: id of a DMRG record (from its /hpc/dmrg/submit "
                         "response) to additionally cross-validate against. SHCI reaches "
                         "its own A/B/C verdict either way — this only adds a delta_mha/"
                         "agreement comparison on top, computed server-side.")
    ap.add_argument("--out", default=".")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    default=None, metavar="API_URL",
                    help="POST the sealed result to SOLANGE (default URL if given bare)")
    args = ap.parse_args()

    sys.path.insert(0, str(Path(args.dice_scripts).expanduser()))
    try:
        from shci import SHCI
    except ImportError as e:
        sys.exit(f"could not import SHCI from {args.dice_scripts}: {e}\n"
                 f"  check --dice-scripts points at the directory containing shci.py")

    from pyscf import gto, scf, ao2mo, mcscf
    from pyscf.mcscf import avas

    aolabels = [s.strip() for s in args.avas.split(",")]
    mol = gto.M(atom=read_xyz(args.geometry), basis=args.basis, charge=args.charge,
                spin=args.spin, verbose=0, max_memory=args.max_memory)
    print(f"cluster: {mol.natm} atoms, {mol.nao} basis functions, "
          f"{mol.nelectron} electrons, charge={args.charge}, spin={args.spin}")

    mf = scf.RHF(mol).density_fit(auxbasis=args.df_auxbasis)
    mf.max_memory = args.max_memory
    mf.kernel()
    if not mf.converged:
        sys.exit("RHF did not converge.")
    print(f"RHF E = {mf.e_tot:.8f} Ha")

    ncas, nelec, mo = avas.avas(mf, aolabels, threshold=args.threshold)
    ncas, nelec = int(ncas), int(nelec)
    occ = nelec // 2
    print(f"AVAS(threshold={args.threshold}) -> CAS({nelec},{ncas})  "
          f"{occ} occupied, {ncas - occ} virtual")

    # Same effective-Hamiltonian construction the DMRG path uses (solange_hpc.py's
    # run_casscf / integrals_from_geometry) -- fixed AVAS orbitals, no
    # optimisation -- so a later cross-validation compares the SAME orbitals a
    # --casci DMRG run used, not a DMRG-SCF-optimised set.
    mc = mcscf.CASCI(mf, ncas, nelec)
    h1e, ecore = mc.get_h1eff(mo_coeff=mo)
    h2e = ao2mo.restore(1, mc.get_h2eff(mo), ncas)
    ecore = float(ecore)

    eps_schedule = [float(x) for x in args.sweep_eps.split(",")]
    if eps_schedule != sorted(eps_schedule, reverse=True):
        sys.exit("--sweep-eps must be given loosest-first, tightest-last (descending values)")
    nelec_ab = ((nelec + args.spin) // 2, (nelec - args.spin) // 2)

    # Run EACH epsilon as its OWN independent solve (not one Dice call given the
    # whole schedule, which returns only the final energy) -- classify_shci()
    # needs the energy AT EACH step to read a convergence trend, the same way
    # solange_dmrg.py's dmrg_energies records one (M, E) pair per bond dimension.
    print(f"SHCI epsilon schedule: {eps_schedule} (run independently, one solve each)")
    energies = []
    t0 = time.time()
    for eps in eps_schedule:
        solver = SHCI(mol)
        solver.sweep_iter = [0]
        solver.sweep_epsilon = [eps]
        solver.num_thrds = args.num_thrds
        t_step = time.time()
        result = solver.kernel(h1e, h2e, ncas, nelec_ab, ecore=ecore)
        # pyscf FCI-solver convention returns (e, civec); some external solvers
        # return e alone -- handle both without guessing which silently.
        e_step = float(result[0] if isinstance(result, (tuple, list)) else result)
        dt = time.time() - t_step
        print(f"  epsilon={eps:g}  E={e_step:.8f} Ha  [{dt:.1f}s]")
        energies.append((eps, e_step))
    elapsed_s = round(time.time() - t0, 1)
    e_shci = energies[-1][1]
    print(f"SHCI E (tightest epsilon) = {e_shci:.8f} Ha")
    print(f"elapsed {elapsed_s}s")

    bqp_class, class_rationale = classify_shci(nelec, energies)
    print("-" * 68)
    print(f"CLASS {bqp_class}")
    print(f"  {class_rationale}")

    out = {
        "id": str(uuid.uuid4()),
        "key": args.key,
        "dmrg_classification_id": args.dmrg_classification_id,  # None -> standalone
        "ncas": ncas, "nelec": nelec,
        "e_shci": e_shci,
        "shci_energies": energies,  # [[eps, E], ...] -- own convergence trace, mirrors dmrg_energies
        "bqp_class": bqp_class, "class_rationale": class_rationale,
        "sweep_eps": ",".join(str(e) for e in eps_schedule),
        "method": "SHCI (Dice, semistochastic heat-bath CI), independent classification",
        "elapsed_s": elapsed_s,
        "provenance_source": "HPC/Laguna (SHCI classifier)",
        "hardware": detect_hardware(args.num_thrds),
    }
    # Seal at source (LEON re-verifies at ingestion — a mismatch is rejected, not
    # trusted). e_dmrg_ref/delta_mha/agreement are NOT computed here on purpose,
    # even when --dmrg-classification-id is given — the backend computes them
    # from that DMRG record's own current stored value (DP1), so this script
    # cannot assert an agreement the server hasn't checked.
    out["shci_seal_payload"] = _seal_payload(out, exclude={"shci_hash", "shci_seal_payload"})
    out["shci_hash"] = hashlib.sha256(out["shci_seal_payload"].encode()).hexdigest()

    Path(args.out).mkdir(parents=True, exist_ok=True)
    p = Path(args.out) / f"shci_class_{args.key}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"WROTE {p}")
    if args.submit:
        _submit_shci(args.submit, out)
    print("=" * 68)
    print("NOTE: SHCI is classical (CPU/MPI), NOT quantum — same as DMRG.")
    if args.dmrg_classification_id:
        print("This run ALSO cross-validates the named DMRG record's own energy at the")
        print("SAME active space (server-side check) -- it does not decide the active")
        print("space, the charge, or the AVAS criterion, which must already match it.")
    else:
        print("Standalone classification -- no DMRG record referenced or required.")


if __name__ == "__main__":
    main()
