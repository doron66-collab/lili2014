#!/bin/bash
#SBATCH --job-name=solange_hpc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=solange_%j.out
#SBATCH --error=solange_%j.err

# ── SOLANGE HPC job for the Laguna cluster ──────────────────────────────────
# Adjust --partition / --gres to match Laguna's actual GPU partition name and
# GPU type once `sinfo -s` / `nvidia-smi` are known. `gpu` matches the status
# page you sent (118 nodes). For the full-site DMRG runs of CDKN2A/ARID2 use
# --partition=largemem (and drop --gres) instead.

set -euo pipefail
echo "Job $SLURM_JOB_ID on $(hostname)  $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# --- environment: adapt to whatever Laguna provides -------------------------
# Try common module names; fall back to an existing conda env. Edit as needed.
module load cuda 2>/dev/null || true
module load anaconda3 2>/dev/null || module load conda 2>/dev/null || true
# source activate solange 2>/dev/null || true      # <- if you make a conda env

# One-time install into a user env if the packages are missing:
#   python -m pip install --user pyscf openfermion scipy numpy
#   python -m pip install --user pennylane pennylane-lightning[gpu]   # for --vqe on GPU

REPO="${SOLANGE_REPO:-$HOME/lili2014}"        # path to the cloned repo on Laguna
cd "$REPO"

# ── The run ─────────────────────────────────────────────────────────────────
# STEP 0 (do this first!): a dry-run that MUST reproduce the known 4q result,
# proving the toolchain is correct before you trust any scaled-up number.
#   python scripts/laguna/solange_hpc.py --compound acetamide --basis sto-3g \
#       --ncas 2 --nelecas 2 --key ARID2_LOF --side native \
#       --residue "Gln1118 dry-run" --out "$REPO/out"

# Scaled-up production run (edit compound/basis/ncas/nelecas per the plan table):
python scripts/laguna/solange_hpc.py \
    --compound "${COMPOUND:-acetamide}" \
    --basis    "${BASIS:-6-31g}" \
    --ncas     "${NCAS:-8}" \
    --nelecas  "${NELECAS:-8}" \
    --key      "${KEY:-ARID2_LOF}" \
    --side     "${SIDE:-native}" \
    --residue  "${RESIDUE:-Gln1118 amide}" \
    --vqe --vqe-steps 80 \
    --out "$REPO/out"

echo "Done $(date). Artifacts in $REPO/out — send both JSON files back to sync."
