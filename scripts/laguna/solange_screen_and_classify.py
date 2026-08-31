#!/usr/bin/env python3
"""
solange_screen_and_classify.py — one job, PDB to notarized classification,
NO terminal interaction from the platform's end user.

This is what the "New Target from PDB" form in the Orchestration tab
dispatches (job_type='screen_classify'). It chains four steps that were,
until now, four separate manual commands run by hand on Laguna over the
course of one long night: protonate.py, build_qm_cluster.py + avas_probe.py
(the same shrink-the-radius-until-it-fits loop screen_target.py already
encoded), solange_dmrg.py --submit, and optionally solange_shci.py --submit.

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
      --radii 5.0,4.0,3.5,3.0 --max-orbitals 45 \\
      --bond-dims 250,500,1000 --casci \\
      --run-shci --dice-scripts ~/lili2014/Dice/scripts \\
      --submit https://qcaihpc-simulation-api.onrender.com
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def run(cmd, label):
    print(f"  $ {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr, file=sys.stderr)
        sys.exit(f"[{label}] failed with exit code {r.returncode}")
    return r.stdout


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
                         "this is also what solange_shci.py assumes, so a later SHCI "
                         "step compares the SAME orbitals)")
    ap.add_argument("--run-shci", action="store_true",
                    help="also run SHCI cross-validation immediately after DMRG, "
                         "against the DMRG record this job just created")
    ap.add_argument("--dice-scripts", default=None, help="required if --run-shci")
    ap.add_argument("--sweep-eps", default="1e-2,1e-3,5e-4,1e-4")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--submit", required=True, help="SOLANGE API base URL")
    args = ap.parse_args()

    if not args.pdb and not args.pdb_id:
        sys.exit("need --pdb (already protonated) or --pdb-id (fetch + protonate)")
    if args.run_shci and not args.dice_scripts:
        sys.exit("--run-shci requires --dice-scripts")

    print("=" * 74)
    print(f"screen_and_classify -- {args.out_prefix}")
    print("=" * 74)

    if args.pdb:
        protonated = args.pdb
        print(f"using existing protonated structure: {protonated}")
    else:
        protonated = f"{args.pdb_id}_protonated.pdb"
        if Path(protonated).exists():
            print(f"{protonated} already exists -- not re-running protonate.py "
                  f"(delete it first to force a redo)")
        else:
            print(f"step 1: protonating {args.pdb_id} at pH {args.ph}")
            print(run(["python3", str(HERE / "protonate.py"), "--pdb-id", args.pdb_id,
                       "--out", protonated, "--ph", str(args.ph)], "protonate.py"))

    print(f"step 2: verifying chain {args.chain} resi {args.resi} is {args.expect_resname}")
    verify_residue(protonated, args.chain, args.resi, args.expect_resname)

    radii = [float(r) for r in args.radii.split(",")]
    print(f"step 3: carving + probing active-space size, radii tried in order: {radii}")
    accepted = None
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

    print("\n" + "=" * 74)
    print(f"step 4: DMRG classification -- CAS({accepted['nelec']},{accepted['ncas']})")
    # Via run_dmrg.sh, NOT solange_dmrg.py directly -- it wraps with with_block2.sh
    # (the LD_PRELOAD fix for block2/MKL), which this job needs regardless of what
    # environment the agent process itself happened to start in.
    dmrg_cmd = ["bash", str(HERE / "run_dmrg.sh"),
                "--geometry", accepted["xyz"], "--charge", str(accepted["charge"]),
                "--spin", str(args.spin), "--basis", args.basis, "--avas", accepted["avas"],
                "--key", args.key, "--bond-dims", args.bond_dims, "--submit", args.submit]
    if args.casci:
        dmrg_cmd.append("--casci")
    dmrg_out = run(dmrg_cmd, "run_dmrg.sh")
    print(dmrg_out)
    dmrg_run_id = run_id_from(dmrg_out)
    if not dmrg_run_id:
        sys.exit("DMRG step ran but no run_id was found in its output -- it may not "
                 "have been stored. See the output above.")
    print(f"DMRG stored: run_id={dmrg_run_id}")

    if args.run_shci:
        print("\n" + "=" * 74)
        print("step 5: SHCI cross-validation against the DMRG record just created")
        shci_cmd = ["python3", str(HERE / "solange_shci.py"),
                    "--geometry", accepted["xyz"], "--charge", str(accepted["charge"]),
                    "--spin", str(args.spin), "--basis", args.basis, "--avas", accepted["avas"],
                    "--key", args.key, "--dmrg-classification-id", dmrg_run_id,
                    "--dice-scripts", args.dice_scripts, "--sweep-eps", args.sweep_eps,
                    "--submit", args.submit]
        shci_out = run(shci_cmd, "solange_shci.py")
        print(shci_out)
        shci_run_id = run_id_from(shci_out)
        print(f"SHCI stored: run_id={shci_run_id}" if shci_run_id
              else "SHCI ran but no run_id was found in its output.")

    print("=" * 74)
    print(f"DONE. dmrg_run_id={dmrg_run_id}")


if __name__ == "__main__":
    main()
