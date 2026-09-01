#!/usr/bin/env python3
"""
solange_screen_and_classify.py — THE CLASSIFIER: PDB to a two-method decision,
NO terminal interaction from the platform's end user.

This is what the "New Target from PDB" block in the Orchestration tab
dispatches (job_type='screen_classify'). Two parts, matching the design of
that block:

  PART 1 — cluster acquisition. If (pdb_id, chain, resi) was already built in
  a previous run, the cached cluster (geometry/AVAS/charge/spin/active space)
  is reused via GET /api/simulate/cluster/lookup — no re-protonation, no
  re-carving, no re-probing. Otherwise this chains protonate.py -> carve at
  each radius via build_qm_cluster.py + avas_probe.py (shrinking the radius
  until the active space fits --max-orbitals) -> saves the accepted cluster
  via POST /api/simulate/cluster/save so the NEXT run of this same site skips
  straight to Part 2.

  PART 2 — classify and wait for both methods to agree. DMRG runs first
  (run_dmrg.sh --submit); SHCI then runs against the SAME active space,
  cross-validated against that DMRG record (solange_shci.py
  --dmrg-classification-id --submit). Both are mandatory, not optional (pass
  --skip-shci only for local testing without a Dice build) — the classifier's
  decision is the two-method outcome, not one method's opinion. The final
  line printed is that outcome: classes from both methods, and whether their
  energies agree, in the same event this endpoint's caller (the agent) already
  streams live.

WHAT THIS DOES NOT DO, on purpose — same limits screen_target.py already
stated, now inherited rather than re-litigated:
  - Does not choose the AVAS chemical criterion for you. A generic default is
    used unless --avas overrides it, and that default has already been shown
    (in this project's own investigation) to pick up unintended atoms on a
    protein-heavy cluster — chemist question 1, unresolved.
  - Does not decide protonation pH, or which oxidation/spin state a
    metal-containing site is in. --charge here is a STARTING POINT taken from
    build_qm_cluster.py's own printed suggestion, the same way screen_target.py
    already used it — not a chemist's sign-off, and Gate 2
    (gate2_requirements.py) is not consulted here at all.
  - Does not retry a failed DMRG/SHCI step with different parameters. A
    failure at any step stops the job and reports why; it does not guess a
    fix and try again.

USAGE (what the agent actually runs — a human could also run this directly)
  python3 solange_screen_and_classify.py --pdb-id 8TDR --chain A --resi 882 \\
      --expect-resname ARG --key DNMT3A_R882 \\
      --radii 5.0,4.0,3.5,3.0 --max-orbitals 45 --bond-dims 250,500,1000 \\
      --submit https://qcaihpc-simulation-api.onrender.com
"""
import argparse
import json
import re
import subprocess
import sys
import urllib.request
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
_RUN_TAG = uuid.uuid4().hex[:8]  # unique per invocation — see --scratch note below


def run(cmd, label):
    print(f"  $ {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr, file=sys.stderr)
        sys.exit(f"[{label}] failed with exit code {r.returncode}")
    return r.stdout


def run_streaming(cmd, label):
    """Like run(), but for the long steps (DMRG, SHCI) — subprocess.run's
    capture_output=True BUFFERS everything until the child exits, so a 7+
    minute DMRG step showed nothing at all in the agent's own live log until
    it finished, indistinguishable from a hang (found the hard way on this
    pipeline's very first real-time-watched run). This streams stdout line by
    line as it arrives, exactly like solange_hpc.py's own _run_with_progress
    does for the agent's top-level subprocess -- the same visibility this
    project has relied on all along to catch failures early, just missing one
    level down inside this orchestrator's own nested subprocess calls."""
    print(f"  $ {' '.join(cmd)}", flush=True)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    lines = []
    for line in iter(proc.stdout.readline, ""):
        print(line, end="", flush=True)
        lines.append(line)
    proc.wait()
    out = "".join(lines)
    if proc.returncode != 0:
        sys.exit(f"[{label}] failed with exit code {proc.returncode}")
    return out


def _api_get(api, path):
    try:
        with urllib.request.urlopen(api.rstrip("/") + path, timeout=30) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"  cluster lookup failed (treating as cache miss): {e}", file=sys.stderr)
        return {}


def _api_post(api, path, body):
    try:
        req = urllib.request.Request(
            api.rstrip("/") + path, data=json.dumps(body, default=str).encode(),
            method="POST", headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"  cluster save failed (non-fatal — this run's cluster just won't "
              f"be reusable next time): {e}", file=sys.stderr)
        return {}


def verify_residue(pdb_path, chain, resi, expect_resname):
    """Read ATOM records with fixed-width PDB columns, not whitespace split --
    whitespace splitting misidentified columns once already in this project's
    own history (the same check that caught 6W8D being a mutant structure)."""
    found = None
    neighbours = []
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        rec_chain = line[21].strip()
        try:
            rec_resi = int(line[22:26])
        except ValueError:
            continue
        rec_resname = line[17:20].strip()
        if rec_chain != chain:
            continue
        if rec_resi == resi:
            found = rec_resname
        if abs(rec_resi - resi) <= 3:
            neighbours.append((rec_resi, rec_resname))
    neighbours = sorted(set(neighbours))
    if found is None:
        sys.exit(f"RESIDUE CHECK FAILED: chain {chain} has nothing at position {resi}.\n"
                 f"  nearby residues: {neighbours}")
    if found.upper() != expect_resname.upper():
        sys.exit(f"RESIDUE CHECK FAILED: chain {chain} resi {resi} is {found}, "
                 f"expected {expect_resname}.\n"
                 f"  This structure may be a mutant, or the numbering may not match "
                 f"what was assumed.\n  nearby residues: {neighbours}")
    print(f"  residue check OK: chain {chain} resi {resi} = {found}")


def carve_and_probe(pdb, chain, resi, radius, spin, basis, avas, threshold, out_prefix):
    out_xyz = f"{out_prefix}_r{radius}.xyz"
    out = run(["python3", str(HERE / "build_qm_cluster.py"), "--pdb", pdb, "--chain", chain,
               "--resi", str(resi), "--radius", str(radius), "--out", out_xyz],
              "build_qm_cluster.py")
    print(out)
    if "HYDROGENS MISSING" in out:
        sys.exit("build_qm_cluster.py reports missing hydrogens on a structure this "
                 "job already protonated — something is wrong upstream. Stopping "
                 "rather than guessing.")
    m_avas = re.search(r'suggested --avas: "([^"]+)"', out)
    m_charge = re.search(r"starting point:\s*--charge\s*(-?\d+)", out)
    suggested_avas = m_avas.group(1) if m_avas else avas
    suggested_charge = int(m_charge.group(1)) if m_charge else 0
    use_avas = avas or suggested_avas

    probe_out = run(["python3", str(HERE / "avas_probe.py"), "--geometry", out_xyz,
                      "--charge", str(suggested_charge), "--spin", str(spin),
                      "--basis", basis, "--avas", use_avas, "--threshold", str(threshold)],
                     "avas_probe.py")
    print(probe_out)
    m = re.search(r"PROBE_RESULT ncas=(\d+) nelec=(\d+) occ=(\d+) virt=(\d+)", probe_out)
    if not m:
        sys.exit("avas_probe.py did not report a parseable result -- see output above.")
    ncas, nelec, occ, virt = (int(x) for x in m.groups())
    return dict(xyz=out_xyz, ncas=ncas, nelec=nelec, occ=occ, virt=virt,
                charge=suggested_charge, avas=use_avas, radius=radius)


def run_id_from(stdout):
    for tok in stdout.split():
        if tok.startswith("run_id="):
            return tok[len("run_id="):]
    return None


def class_from(stdout):
    m = re.search(r"^CLASS (\w)", stdout, re.MULTILINE)
    return m.group(1) if m else None


def agreement_from(stdout):
    m = re.search(r"agreement=(True|False)", stdout)
    return m.group(1) == "True" if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-id", help="4-char RCSB ID to fetch and protonate")
    ap.add_argument("--pdb", help="already-protonated local .pdb -- skips the protonate step")
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--resi", type=int, required=True)
    ap.add_argument("--expect-resname", required=True,
                    help="e.g. ARG -- what the residue at --chain/--resi must actually be")
    ap.add_argument("--radii", default="5.0,4.0,3.5,3.0",
                    help="tried in order, largest first, until AVAS active-space size "
                         "is acceptable")
    ap.add_argument("--max-orbitals", type=int, default=45)
    ap.add_argument("--avas", default=None,
                    help="omit to let each radius use build_qm_cluster.py's own "
                         "per-radius suggestion instead of one fixed criterion")
    ap.add_argument("--threshold", type=float, default=0.2)
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--key", required=True)
    ap.add_argument("--bond-dims", default="250,500,1000")
    ap.add_argument("--casci", action="store_true", default=True,
                    help="fixed AVAS orbitals, no CASSCF optimisation (default ON -- "
                         "this is also what solange_shci.py assumes, so the SHCI step "
                         "compares the SAME orbitals)")
    ap.add_argument("--skip-shci", action="store_true",
                    help="skip the SHCI leg (local testing without a Dice build only -- "
                         "the classifier's normal outcome needs BOTH methods, not one)")
    ap.add_argument("--dice-scripts", default=None,
                    help="defaults to <repo>/Dice/scripts next to this script")
    ap.add_argument("--sweep-eps", default="1e-2,1e-3,5e-4,1e-4")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--submit", required=True, help="SOLANGE API base URL")
    args = ap.parse_args()

    if not args.pdb and not args.pdb_id:
        sys.exit("need --pdb (already protonated) or --pdb-id (fetch + protonate)")
    dice_scripts = args.dice_scripts or str(HERE.parent.parent / "Dice" / "scripts")

    print("=" * 74)
    print(f"CLASSIFIER -- {args.out_prefix}")
    print("=" * 74)

    # ── PART 1 — cluster acquisition: reuse if this exact site was ever built
    # before, otherwise protonate/carve/probe from scratch and save the result.
    accepted = None
    if args.pdb_id:
        cached = _api_get(args.submit, f"/api/simulate/cluster/lookup?pdb_id={args.pdb_id}"
                                        f"&chain={args.chain}&resi={args.resi}").get("cluster")
        if cached and cached.get("expect_resname", "").upper() == args.expect_resname.upper():
            print(f"PART 1: cluster cache HIT for {args.pdb_id}/{args.chain}/{args.resi} "
                  f"(built {cached.get('created_at', '?')}, radius {cached.get('radius')}) "
                  f"-- skipping protonate/carve/probe entirely.")
            out_xyz = f"{args.out_prefix}_cached.xyz"
            Path(out_xyz).write_text(cached["geometry"])
            accepted = dict(xyz=out_xyz, ncas=cached["ncas"], nelec=cached["nelec"],
                            charge=cached["charge"], avas=cached["avas"],
                            radius=cached.get("radius"))
        elif cached:
            print(f"PART 1: cluster cache entry exists for {args.pdb_id}/{args.chain}/{args.resi} "
                  f"but its expect_resname ({cached.get('expect_resname')}) does not match this "
                  f"request ({args.expect_resname}) -- treating as a cache miss, not trusting it.")

    if not accepted:
        print("PART 1: no usable cached cluster -- building from scratch.")
        if args.pdb:
            protonated = args.pdb
            print(f"using existing protonated structure: {protonated}")
        else:
            protonated = f"{args.pdb_id}_protonated.pdb"
            if Path(protonated).exists():
                print(f"{protonated} already exists -- not re-running protonate.py "
                      f"(delete it first to force a redo)")
            else:
                print(f"step 1a: protonating {args.pdb_id} at pH {args.ph}")
                print(run(["python3", str(HERE / "protonate.py"), "--pdb-id", args.pdb_id,
                           "--out", protonated, "--ph", str(args.ph)], "protonate.py"))

        print(f"step 1b: verifying chain {args.chain} resi {args.resi} is {args.expect_resname}")
        verify_residue(protonated, args.chain, args.resi, args.expect_resname)

        radii = [float(r) for r in args.radii.split(",")]
        print(f"step 1c: carving + probing active-space size, radii tried in order: {radii}")
        for radius in radii:
            print(f"\n-- radius {radius} A --")
            result = carve_and_probe(protonated, args.chain, args.resi, radius, args.spin,
                                      args.basis, args.avas, args.threshold, args.out_prefix)
            print(f"  -> CAS({result['nelec']},{result['ncas']}), "
                  f"{result['occ']} occ / {result['virt']} virt")
            if result["ncas"] <= args.max_orbitals:
                accepted = result
                print(f"  ACCEPTED: {result['ncas']} orbitals <= --max-orbitals {args.max_orbitals}")
                break
            print(f"  REJECTED: {result['ncas']} orbitals > --max-orbitals {args.max_orbitals}")

        if not accepted:
            sys.exit("NO RADIUS in the tried list gave an acceptable active-space size. "
                     "This generic AVAS criterion may simply be wrong for this site -- "
                     "not a radius problem.")

        if args.pdb_id:
            save_resp = _api_post(args.submit, "/api/simulate/cluster/save", {
                "pdb_id": args.pdb_id, "chain": args.chain, "resi": args.resi,
                "expect_resname": args.expect_resname, "key": args.key,
                "geometry": Path(accepted["xyz"]).read_text(),
                "avas": accepted["avas"], "charge": accepted["charge"], "spin": args.spin,
                "ncas": accepted["ncas"], "nelec": accepted["nelec"], "radius": accepted["radius"],
            })
            print(f"  cluster cache: {'saved' if save_resp.get('saved') else 'NOT saved (' + str(save_resp.get('error')) + ')'} "
                  f"-- {'a later run of this exact site will reuse it' if save_resp.get('saved') else 'next run will rebuild from scratch'}")

    # ── PART 2 — classify with both methods, wait for their conclusion ──────
    print("\n" + "=" * 74)
    print(f"PART 2: classifying -- CAS({accepted['nelec']},{accepted['ncas']})")
    dmrg_cmd = ["bash", str(HERE / "run_dmrg.sh"),
                "--geometry", accepted["xyz"], "--charge", str(accepted["charge"]),
                "--spin", str(args.spin), "--basis", args.basis, "--avas", accepted["avas"],
                "--key", args.key, "--bond-dims", args.bond_dims, "--submit", args.submit,
                # --scratch keyed by a per-invocation random tag: solange_dmrg.py's own
                # default ("./tmp_dmrg") is one FIXED shared directory across every run
                # in this working directory, and its "[resume] loaded existing MPS"
                # logic then tries to continue from whatever a PRIOR, possibly
                # incompatible run left there -- which segfaults block2's C++ backend
                # rather than failing cleanly. Found the hard way on this pipeline's
                # first real run.
                "--scratch", f"./tmp_dmrg_{args.out_prefix}_{_RUN_TAG}"]
    if args.casci:
        dmrg_cmd.append("--casci")
    dmrg_out = run_streaming(dmrg_cmd, "run_dmrg.sh")  # streamed live already, not re-printed
    dmrg_run_id = run_id_from(dmrg_out)
    dmrg_class = class_from(dmrg_out)
    if not dmrg_run_id:
        sys.exit("DMRG step ran but no run_id was found in its output -- it may not "
                 "have been stored. See the output above.")
    print(f"DMRG stored: run_id={dmrg_run_id}  class={dmrg_class}")

    shci_class = shci_run_id = agreement = None
    if not args.skip_shci:
        print("\n" + "=" * 74)
        print("classifying with the second method -- SHCI, cross-validated against the DMRG record above")
        shci_cmd = ["python3", str(HERE / "solange_shci.py"),
                    "--geometry", accepted["xyz"], "--charge", str(accepted["charge"]),
                    "--spin", str(args.spin), "--basis", args.basis, "--avas", accepted["avas"],
                    "--key", args.key, "--dmrg-classification-id", dmrg_run_id,
                    "--dice-scripts", dice_scripts, "--sweep-eps", args.sweep_eps,
                    "--submit", args.submit]
        shci_out = run_streaming(shci_cmd, "solange_shci.py")  # streamed live already, not re-printed
        shci_run_id = run_id_from(shci_out)
        shci_class = class_from(shci_out)
        agreement = agreement_from(shci_out)
        print(f"SHCI stored: run_id={shci_run_id}  class={shci_class}" if shci_run_id
              else "SHCI ran but no run_id was found in its output.")

    print("\n" + "=" * 74)
    print("CLASSIFIER OUTCOME")
    print(f"  DMRG : class={dmrg_class}  run_id={dmrg_run_id}")
    if args.skip_shci:
        print("  SHCI : skipped (--skip-shci) -- this is a single-method result, not the "
              "classifier's normal two-method outcome")
    else:
        print(f"  SHCI : class={shci_class}  run_id={shci_run_id}")
        if agreement is True:
            print("  >>> METHODS AGREE <<<")
        elif agreement is False:
            print("  >>> ⚠ METHODS DISAGREE -- do not treat this classification as settled. "
                  "See the SHCI record's own delta_mha against the 1.6 mHa chemical-accuracy "
                  "bar for the size of the disagreement. <<<")
        else:
            print("  agreement: unknown (SHCI's own energy comparison was not found in its output)")
    print("=" * 74)
    print(f"DONE. dmrg_run_id={dmrg_run_id}")


if __name__ == "__main__":
    main()
