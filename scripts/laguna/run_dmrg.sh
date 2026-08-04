#!/usr/bin/env bash
# run_dmrg.sh — run the SOLANGE DMRG classifier with the block2 / MKL environment
# set up AUTOMATICALLY, so the "cannot open libmkl_def.so.1" and "undefined symbol:
# omp_get_num_procs" errors (the packaging battle) never come back.
#
# That environment setup now lives in with_block2.sh, because it is needed by more
# than this script (solange_hpc.py --dmrg-scf and dmrgscf_block2.py's validation
# import block2 too, and fail identically without it). This stays as the familiar
# entry point for the classifier; every argument is passed straight through to
# solange_dmrg.py.
#
#   bash scripts/laguna/run_dmrg.sh --key CDKN2A_LOF --side native \
#       --ncas 10 --nelecas 10 --bond-dims 250,500,1000,2000 --submit
#
# Set SOLANGE_ENV / SOLANGE_REPO to match Laguna (see with_block2.sh).
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -u = unbuffered: the per-bond-dimension "DMRG M=" progress lines flush live to
# the agent's stdout stream instead of block-buffering until exit (which hid them
# from the live view). PYTHONUNBUFFERED can get reset across the conda activation
# inside with_block2.sh, so pass -u explicitly here.
exec bash "$HERE/with_block2.sh" python -u scripts/laguna/solange_dmrg.py "$@"
