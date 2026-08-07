#!/bin/bash
#SBATCH --job-name=solange_dmrg
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=dmrg_%j.out
#SBATCH --error=dmrg_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=doron.cohen@cgu.edu

# ── SOLANGE DMRG classifier job for USC Laguna (CARC) ───────────────────────
# Runs solange_dmrg.py's --geometry path (real active-site cluster, not a demo
# compound) as an INDEPENDENT batch job — unlike `salloc`, this does NOT die
# when your browser/SSH session disconnects. Submit it, close everything, come
# back later and check with `squeue --me` / `tail -f dmrg_<jobid>.out`.
#
# Defaults below match the TP53_C275F validation run (5A cluster around
# Cys275, sto-3g). Density fitting is ON by default with an explicit modest
# auxbasis (def2-universal-jkfit) for both RHF and CASSCF — an earlier
# --no-density-fit default here fixed a memory crash in RHF but made CASSCF
# itself take hours for even a tiny active space (it still transforms every
# basis function each macro-iteration without DF). The real fix for the
# original OOM was pinning the auxbasis, not dropping DF — see
# solange_dmrg.py's --df-auxbasis. Set NO_DENSITY_FIT_FLAG=--no-density-fit
# only if an exact, non-DF comparison is specifically needed.
# Override any of them via environment variables, e.g.:
#   sbatch --export=RADIUS_XYZ=tp53_c275_r7.0_h.xyz,BASIS=6-31g run_dmrg_geometry.sbatch
#
# DMRG_SCF=1 (set below) uses DMRG (block2) as CASSCF's own orbital-optimization
# solver instead of FCI — see scripts/laguna/dmrgscf_block2.py. This is what
# actually raises the ncas ceiling past ~16 orbitals; the FCI path this script
# used before cannot reach a real site's active-space size at all. Unset it
# (--export=DMRG_SCF=) to fall back to plain FCI-based CASSCF.
#
# MAX_MINUTES (below) matters as much as --time does. solange_dmrg.py's own
# run_dmrg() already knows how to stop gracefully and save a PROVISIONAL
# classification once a wall-clock budget is spent — but only if TOLD that
# budget via --max-minutes. Without it (the state this template shipped in
# until now), the only limit is SBATCH's --time, which SLURM enforces with a
# hard SIGKILL: no result written, nothing submitted, hours of RHF+AVAS+CASSCF
# gone with it, exactly what happened on the first real-geometry run (TP53
# Cys275 cluster, CAS(18,12)) — it never got past validate() in the surviving
# log because the job was killed partway through the bond-dimension sweep,
# not because anything was wrong scientifically. Default here (215 min) is
# --time's 240 with a 25-minute margin for the final write + --submit to
# complete before SLURM's own clock runs out; adjust both together if you
# override --time.
#
# --scratch's default was tied to $SLURM_JOB_ID until now, which meant a job
# stopped by --max-minutes (or, before that fix, killed outright by --time)
# could never actually resume on resubmission — a fresh sbatch call gets a
# fresh job ID, hence a fresh, empty scratch dir, hence bond_dims[0] again
# every time, exactly the outcome run_dmrg()'s own resume-by-reusing-scratch
# design was meant to avoid. Now keyed on the target instead, so re-running
# this same template for the same KEY picks the MPS up where it left off.

set -euo pipefail
echo "Job $SLURM_JOB_ID on $(hostname)  $(date)"
free -h

# --- environment ------------------------------------------------------------
module load conda 2>/dev/null || true
source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${SOLANGE_ENV:-base}" 2>/dev/null || true
python -c "import block2, pyscf" || { echo "block2/pyscf not importable — check conda env"; exit 1; }

# Fewer threads = lower peak memory (traded for speed) — this is what got the
# --no-density-fit run stable yesterday. Raise if a future run has memory
# headroom to spare and you want it faster instead.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-2}"

REPO="${SOLANGE_REPO:-$HOME/lili2014}"
cd "$REPO"

# --dmrg-scf and its two sub-flags are appended as an array, not inlined into
# the fixed argument list above them, so DMRG_SCF unset means the call below
# is byte-identical to before this option existed (no stray empty strings
# reaching argparse, which an inline "${X:-}" pattern would risk).
DMRG_SCF_ARGS=()
if [[ -n "${DMRG_SCF:-}" ]]; then
  DMRG_SCF_ARGS=(--dmrg-scf
                 --dmrg-scf-maxm "${DMRG_SCF_MAXM:-250}"
                 --dmrg-scf-scratch "${DMRG_SCF_SCRATCH:-./tmp_dmrgscf_orb_$SLURM_JOB_ID}")
fi

# ── The run ──────────────────────────────────────────────────────────────────
bash scripts/laguna/run_dmrg.sh \
    --key       "${KEY:-TP53_C275F}" \
    --side      "${SIDE:-native}" \
    --geometry  "${GEOMETRY:-tp53_c275_r5.0_h.xyz}" \
    --avas      "${AVAS:-S 3p}" \
    --basis     "${BASIS:-sto-3g}" \
    --charge    "${CHARGE:-0}" \
    --verbose   "${VERBOSE:-3}" \
    --max-memory "${MAX_MEMORY:-16000}" \
    --df-auxbasis "${DF_AUXBASIS:-def2-universal-jkfit}" \
    ${NO_DENSITY_FIT_FLAG:-} \
    --scratch   "${SCRATCH:-./tmp_dmrg_${KEY:-TP53_C275F}}" \
    --bond-dims "${BOND_DIMS:-250,500,1000,2000}" \
    --max-minutes "${MAX_MINUTES:-215}" \
    --stack-mem-gb "${STACK_MEM_GB:-4}" \
    "${DMRG_SCF_ARGS[@]}" \
    --submit

echo "Done $(date)."
