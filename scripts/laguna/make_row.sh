#!/usr/bin/env bash
# make_row.sh — one-command Laguna run that creates a row in SOLANGE.
#
# What it does: runs a real CASSCF + (optional) statevector VQE on the cluster,
# seals a P1–P9 provenance record, and --submits it to the SOLANGE backend, where
# LEON re-verifies the P8 seal before the row is stored. The row then appears in
# Orchestration → "Ladder · Rung 2 — HPC" (Phase 3A-HPC).
#
# Adjust the two lines marked EDIT to match your Laguna setup, then:
#   bash scripts/laguna/make_row.sh          # real run → creates the row
#   bash scripts/laguna/make_row.sh --dry     # 5-sec sanity check, does NOT touch SOLANGE
#
# The only proof the row is really stored is the line "db=stored" in the output.
set -euo pipefail

# ── EDIT these two to your environment ───────────────────────────────────────
CONDA_ENV="${SOLANGE_ENV:-solange}"     # EDIT: conda env that has pyscf + pennylane + numpy
REPO_DIR="${SOLANGE_REPO:-$HOME/lili2014}"   # EDIT: where you cloned the repo on Laguna
# ─────────────────────────────────────────────────────────────────────────────

# The chemistry (a validated, known-good target — ARID2_LOF native, acetamide/Gln1118).
KEY="ARID2_LOF"; SIDE="native"; COMPOUND="acetamide"
BASIS="6-31g"; NCAS=6; NELECAS=6        # CAS(6,6) → 12 qubits: fast, and VQE fits statevector
RESIDUE="Gln1118 amide"

# module load conda is the step most often forgotten (→ "No module named numpy").
module load conda 2>/dev/null || true
# shellcheck disable=SC1091
source activate "$CONDA_ENV" 2>/dev/null || conda activate "$CONDA_ENV"
cd "$REPO_DIR"

if [[ "${1:-}" == "--dry" ]]; then
  echo "── SANITY DRY-RUN (CAS(2,2)/sto-3g, no --submit) ───────────────────────"
  python scripts/laguna/solange_hpc.py \
    --compound "$COMPOUND" --basis sto-3g --ncas 2 --nelecas 2 \
    --key "$KEY" --side "$SIDE" --residue "$RESIDUE dry-run" --out ./out
  echo "Dry-run OK — the environment works. Re-run without --dry to create the row."
  exit 0
fi

echo "── REAL RUN → will create a Phase 3A-HPC row in SOLANGE ────────────────"
python scripts/laguna/solange_hpc.py \
  --compound "$COMPOUND" --basis "$BASIS" --ncas "$NCAS" --nelecas "$NELECAS" \
  --key "$KEY" --side "$SIDE" --residue "$RESIDUE" \
  --vqe --out ./out --submit

echo
echo "Look for  db=stored  and  'matched → verified'  above."
echo "Then open SOLANGE → Orchestration → Rung 2 (Phase 3A-HPC). The live"
echo "auto-refresh should show the new row within a few seconds; if not, reload."
