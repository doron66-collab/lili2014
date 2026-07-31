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

# ── SOLANGE DMRG classifier job for USC Laguna (CARC) ───────────────────────
# Runs solange_dmrg.py's --geometry path (real active-site cluster, not a demo
# compound) as an INDEPENDENT batch job — unlike `salloc`, this does NOT die
# when your browser/SSH session disconnects. Submit it, close everything, come
# back later and check with `squeue --me` / `tail -f dmrg_<jobid>.out`.
#
# Defaults below match the TP53_C275F validation run from 2026-07-30 (5A
# cluster around Cys275, sto-3g, no density-fit — that combination is what
# actually got RHF to converge without being OOM-killed on this cluster).
# Override any of them via environment variables, e.g.:
#   sbatch --export=RADIUS_XYZ=tp53_c275_r7.0_h.xyz,BASIS=6-31g run_dmrg_geometry.sbatch

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
    ${NO_DENSITY_FIT_FLAG:---no-density-fit} \
    --scratch   "${SCRATCH:-./tmp_dmrg_$SLURM_JOB_ID}" \
    --bond-dims "${BOND_DIMS:-250,500,1000,2000}" \
    --submit

echo "Done $(date)."
