#!/bin/bash
#SBATCH --job-name=sdhb_s3
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=sdhb_s3_%j.out
#SBATCH --error=sdhb_s3_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=doron.cohen@cgu.edu

# ── SDHB S3 [3Fe-4S] Gate 2 run as an independent batch job ─────────────────
# Why this exists: the previous run was launched with `nohup` inside an Open
# OnDemand Jupyter session. That session is itself a Slurm job with an 8-hour
# limit (job 2350784: TIMEOUT at 08:00:00), and when it ended Slurm killed
# everything inside it -- nohup only protects against a closed terminal, not
# against the enclosing job ending. A batch job submitted with sbatch lives
# on its own: submit it, close the browser, come back later.
#
#   cd ~/lili2014/scripts/laguna
#   sbatch run_gate2_sdhb_s3.sbatch
#   squeue --me
#   grep -E "\[scf\]|REFUSING|embedding-hs|DMRG M=" sdhb_s3_<jobid>.out
#
# Resubmitting the same command resumes, it does not restart: the SCF
# restarts from sdhb_s3_scf.chk, and each DMRG leg resumes its own MPS from
# its own scratch directory (tmp_dmrg_s3_hs / tmp_dmrg_s3).
#
# Budget: --time is a hard SIGKILL with nothing written. DMRG_MAX_MINUTES is
# a PER-LEG soft budget (run_dmrg() stops requesting larger bond dimensions
# once it is spent and still writes its result), so keep
#   SCF (~2-3 h) + 2 x DMRG_MAX_MINUTES  <  --time
# Default ladder is the protocol minimum (250,500,1000); extend to 2000 later
# by resubmitting with BOND_DIMS=250,500,1000,2000 -- the saved MPS resumes.
# If the partition rejects --time=24:00:00, sbatch says so immediately;
# lower --time and DMRG_MAX_MINUTES together.

set -uo pipefail
echo "Job ${SLURM_JOB_ID:-?} on $(hostname)  $(date)"
free -h

THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

REPO="${SOLANGE_REPO:-$HOME/lili2014}"
cd "$REPO/scripts/laguna" || { echo "repo not found: $REPO"; exit 1; }

# with_block2.sh sets up the block2/MKL environment (see that file) and then
# cd's to the repo root, so the run itself re-enters scripts/laguna: this
# script reads struct/, sifts/, mmcif_verify.py and its checkpoint files from
# there.
bash "$REPO/scripts/laguna/with_block2.sh" bash -c '
  cd scripts/laguna && exec python -u run_gate2_sdhb_s3.py "$@"' _ \
    --bond-dims "${BOND_DIMS:-250,500,1000}" \
    --dmrg-threads "$THREADS" \
    --dmrg-max-minutes "${DMRG_MAX_MINUTES:-540}" \
    ${EXTRA_ARGS:-}

echo "Done $(date)."
