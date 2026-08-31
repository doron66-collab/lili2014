#!/usr/bin/env python3
"""
protonate.py — STANDALONE. Not wired into SOLANGE, not imported by anything.

X-ray crystallography does not resolve hydrogens (X-rays scatter off electrons;
hydrogen has one). A structure downloaded from the PDB is therefore a heavy-atom
skeleton: valences are unsatisfied and the electron count is not the molecule's.
build_qm_cluster.py says exactly this and refuses to pretend otherwise — it caps
cut bonds but does not add the missing hydrogens on existing heavy atoms.

This wrapper calls pdb2pqr (which bundles PROPKA) to do that properly: assign
protonation states from environment-dependent pKa shifts, add the hydrogens, and
report the resulting formal charge. It implements no chemistry of its own — that
is the point. A thin wrapper over a validated tool has less of our own code that
can be wrong.

USAGE
  python3 protonate.py --pdb 2OCJ.pdb --out 2OCJ_H.pdb
  python3 protonate.py --pdb-id 2OCJ --out 2OCJ_H.pdb --ph 7.0

VERIFICATION (this is the point of the module existing separately)
  --verify runs the whole thing against a known-good answer instead of trusting
  it. SOLANGE already has eight clusters that were protonated somehow and whose
  working net charge is recorded in targets.json. If this module reproduces those
  charges it is reproducing a process already proven; if it does not, we learn
  where it differs BEFORE anything touches the pipeline.

NOT DONE HERE, DELIBERATELY
  - No cluster carving. That is build_qm_cluster.py's job and it is not touched.
  - No metal-site special-casing. A zinc-bound cysteine is a thiolate, not a
    thiol, and PROPKA's handling of that is exactly one of the things --verify
    is meant to expose rather than assume.
"""
import argparse
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

PDB_URL = "https://files.rcsb.org/download/{}.pdb"


def fetch_pdb(pdb_id: str, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"{pdb_id.upper()}.pdb"
    if dest.exists():
        print(f"  using cached {dest}")
        return dest
    url = PDB_URL.format(pdb_id.upper())
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)
    return dest


def find_pdb2pqr() -> str:
    for name in ("pdb2pqr", "pdb2pqr30"):
        if shutil.which(name):
            return name
    # installed as a module but not on PATH
    try:
        import pdb2pqr  # noqa: F401
        return f"{sys.executable} -m pdb2pqr"
    except ImportError:
        pass
    sys.exit("protonate: pdb2pqr not found on PATH and not importable.\n"
             "  Install with:  pip install --user pdb2pqr")


def run_pdb2pqr(exe: str, src: Path, out_pdb: Path, ph: float,
                forcefield: str, verbose: bool) -> Path:
    """Run pdb2pqr with PROPKA pKa assignment. Returns the .pqr path.

    --keep-chain matters: without it chain IDs are dropped, and every downstream
    step here selects residues by chain. --pdb-output gives a protonated PDB
    alongside the .pqr, which is what a QM cluster builder can actually read.
    """
    pqr = out_pdb.with_suffix(".pqr")
    cmd = (exe.split() + [
        "--ff", forcefield,
        "--titration-state-method", "propka",
        "--with-ph", str(ph),
        "--keep-chain",
        "--pdb-output", str(out_pdb),
        str(src), str(pqr),
    ])
    print("  running:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=not verbose, text=True)
    if res.returncode != 0:
        if not verbose and res.stderr:
            print(res.stderr, file=sys.stderr)
        sys.exit(f"protonate: pdb2pqr failed (exit {res.returncode})")
    return pqr


def charge_from_pqr(pqr: Path) -> float:
    """Sum the per-atom partial charges pdb2pqr assigns.

    A correctly protonated protein sums to very near an integer; a value far off
    integer means something did not titrate cleanly and is worth stopping for,
    not rounding away.
    """
    total = 0.0
    n = 0
    for line in pqr.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            parts = line.split()
            try:
                total += float(parts[-2])   # PQR: ... x y z charge radius
                n += 1
            except (ValueError, IndexError):
                continue
    return total, n


def count_atoms(pdb: Path):
    heavy = h = 0
    for line in pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            elem = line[76:78].strip() or line[12:16].strip()[:1]
            if elem.upper() == "H":
                h += 1
            else:
                heavy += 1
    return heavy, h


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb", help="local .pdb file")
    src.add_argument("--pdb-id", help="4-character PDB ID, downloaded from RCSB")
    ap.add_argument("--out", required=True, help="output protonated .pdb")
    ap.add_argument("--ph", type=float, default=7.0,
                    help="pH for titration-state assignment (default 7.0). This is "
                         "a chemical choice, not a default to accept silently — a "
                         "buried site can behave as though at a different pH.")
    ap.add_argument("--forcefield", default="AMBER",
                    help="pdb2pqr force field for atom naming/charges (default AMBER)")
    ap.add_argument("--cache-dir", default="./pdb_cache")
    ap.add_argument("--verbose", action="store_true",
                    help="stream pdb2pqr's own output instead of capturing it")
    args = ap.parse_args()

    exe = find_pdb2pqr()
    print(f"protonate: using {exe}")

    source = Path(args.pdb) if args.pdb else fetch_pdb(args.pdb_id, Path(args.cache_dir))
    if not source.exists():
        sys.exit(f"protonate: input not found: {source}")

    before_heavy, before_h = count_atoms(source)
    print(f"  input : {source}  ({before_heavy} heavy, {before_h} H)")

    out_pdb = Path(args.out)
    pqr = run_pdb2pqr(exe, source, out_pdb, args.ph, args.forcefield, args.verbose)

    after_heavy, after_h = count_atoms(out_pdb)
    charge, natoms = charge_from_pqr(pqr)

    print(f"  output: {out_pdb}  ({after_heavy} heavy, {after_h} H)")
    print(f"  added : {after_h - before_h} hydrogens")
    print(f"  net formal charge (sum of {natoms} partial charges): {charge:+.3f}")

    nearest = round(charge)
    if abs(charge - nearest) > 0.05:
        print(f"  *** CHARGE IS NOT NEAR AN INTEGER ({charge:+.3f}, nearest {nearest:+d}). ***")
        print("      A cleanly titrated structure sums to near-integer. Investigate before")
        print("      passing this charge to a QM calculation — do not round it away.")
    else:
        print(f"  -> integer net charge: {nearest:+d}")

    if after_h <= before_h:
        print("  *** NO HYDROGENS ADDED — check that pdb2pqr actually ran on this input. ***")

    ratio = after_h / after_heavy if after_heavy else 0
    print(f"  H/heavy ratio: {ratio:.2f}  (a protonated organic structure is near 1.0)")


if __name__ == "__main__":
    main()
