#!/usr/bin/env python3
"""avas_probe.py -- STANDALONE. RHF + AVAS only, no DMRG.

Answers one question fast, in seconds rather than an sbatch job: how big is
the active space this cluster and this AVAS criterion actually produce?

Built because DNMT3A R882 at radius 5.0 gave CAS(144,81) -- discovered only
after paying for a full sbatch submission. This checks the same thing before
any DMRG time is spent.

USAGE
  python3 avas_probe.py --geometry cluster.xyz --charge 1 --basis sto-3g \\
      --avas "N 2p, O 2p" --threshold 0.2
"""
import argparse
import sys
from pathlib import Path


def read_xyz(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    if lines and lines[0].strip().isdigit():
        lines = lines[2:]
    return "\n".join(lines)


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
    args = ap.parse_args()

    from pyscf import gto, scf
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
        sys.exit("RHF did not converge -- fix charge/spin/protonation before probing AVAS.")
    print(f"RHF E = {mf.e_tot:.8f} Ha")

    ncas, nelec, mo = avas.avas(mf, aolabels, threshold=args.threshold)
    occ = nelec // 2
    virt = ncas - occ
    print(f"AVAS(threshold={args.threshold}) -> CAS({nelec},{ncas})  "
          f"{occ} occupied, {virt} virtual ({100*virt/ncas:.1f}%)")
    # Machine-readable line for the orchestrator to parse.
    print(f"PROBE_RESULT ncas={ncas} nelec={nelec} occ={occ} virt={virt}")


if __name__ == "__main__":
    main()
