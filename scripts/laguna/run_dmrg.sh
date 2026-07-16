#!/usr/bin/env bash
# run_dmrg.sh — run the SOLANGE DMRG classifier with the block2 / MKL environment
# set up AUTOMATICALLY, so the "cannot open libmkl_def.so.1" and "undefined symbol:
# omp_get_num_procs" errors (the packaging battle) never come back.
#
# It (1) activates your env + cd's to the repo, (2) finds block2's bundled MKL
# libraries, (3) makes sure libmkl_def.so.1 is present (the wheel omits it), and
# (4) builds the exact LD_PRELOAD that made block2 load, then runs the classifier.
# Every argument is passed straight through to solange_dmrg.py.
#
#   bash scripts/laguna/run_dmrg.sh --key CDKN2A_LOF --side native \
#       --ncas 10 --nelecas 10 --bond-dims 250,500,1000,2000 --submit
#
# Set SOLANGE_ENV / SOLANGE_REPO (or edit the two lines below) to match Laguna.
set -uo pipefail

CONDA_ENV="${SOLANGE_ENV:-solange}"           # env with block2/pyblock2 + pyscf
REPO_DIR="${SOLANGE_REPO:-$HOME/lili2014}"    # repo checkout on Laguna

module load conda 2>/dev/null || true
# shellcheck disable=SC1091
source activate "$CONDA_ENV" 2>/dev/null || conda activate "$CONDA_ENV" || {
  echo "run_dmrg: could not activate conda env '$CONDA_ENV' — set SOLANGE_ENV." >&2; exit 1; }
cd "$REPO_DIR" || { echo "run_dmrg: repo not found: $REPO_DIR — set SOLANGE_REPO." >&2; exit 1; }

# 1) locate block2's bundled library dir (block2.libs, created by the wheel) ────
BLK_LIBS="$(python - <<'PY'
import os, glob, importlib.util
spec = importlib.util.find_spec('block2')
if not spec or not spec.origin:
    print(''); raise SystemExit
site = os.path.dirname(os.path.dirname(spec.origin))          # .../site-packages
cands = glob.glob(os.path.join(site, 'block2.libs')) \
      + glob.glob(os.path.join(site, 'block2', '.libs'))
print(cands[0] if cands else '')
PY
)"
if [[ -z "$BLK_LIBS" || ! -d "$BLK_LIBS" ]]; then
  echo "run_dmrg: block2 not found in env '$CONDA_ENV'.  Install it:  pip install --user block2" >&2
  exit 1
fi
echo "run_dmrg: block2 libs → $BLK_LIBS"

# 2) ensure libmkl_def.so.1 exists in block2.libs (the wheel ships hashed MKL libs
#    but omits this one; without it block2 fails 'cannot open libmkl_def.so.1') ──
if ! ls "$BLK_LIBS"/libmkl_def.so.1 >/dev/null 2>&1; then
  SRC="$(python - <<'PY'
import glob, os, site
roots = []
p = os.environ.get('CONDA_PREFIX')
if p: roots.append(os.path.join(p, 'lib'))
try: roots += site.getsitepackages()
except Exception: pass
for r in roots:
    hits = glob.glob(os.path.join(r, '**', 'libmkl_def.so.1'), recursive=True)
    if hits: print(hits[0]); break
PY
)"
  if [[ -n "$SRC" ]]; then
    cp -n "$SRC" "$BLK_LIBS"/ && echo "run_dmrg: copied libmkl_def.so.1 ← $SRC"
  else
    echo "run_dmrg: WARNING libmkl_def.so.1 not found. If block2 errors on it, run:" >&2
    echo "         pip install --user mkl==2021.4.0   (then re-run — this copies it in)" >&2
  fi
fi

# 3) build LD_PRELOAD: block2's own intel_lp64 + gnu_thread + core, plus the GNU
#    OpenMP runtime (libgomp) — the latter fixes 'undefined symbol: omp_get_num_procs'
#    when only the MKL trio is preloaded. Order matters (lp64 → thread → core). ──
pick() { ls "$BLK_LIBS"/$1 2>/dev/null | head -n1; }
LP64="$(pick 'libmkl_intel_lp64*.so*')"
THRD="$(pick 'libmkl_gnu_thread*.so*')"
CORE="$(pick 'libmkl_core*.so*')"
GOMP="$(ls /usr/lib*/libgomp.so.1 /lib*/libgomp.so.1 2>/dev/null | head -n1)"
[[ -z "$GOMP" ]] && GOMP="$(find "${CONDA_PREFIX:-/nonexistent}" -name 'libgomp.so*' 2>/dev/null | head -n1)"

PRELOAD=""
for L in "$LP64" "$THRD" "$CORE" "$GOMP"; do
  [[ -n "$L" ]] && PRELOAD="${PRELOAD:+$PRELOAD:}$L"
done
if [[ -z "$LP64" || -z "$CORE" ]]; then
  echo "run_dmrg: WARNING could not find the MKL lp64/core libs in block2.libs;" >&2
  echo "         running without LD_PRELOAD — if it errors, that is the cause." >&2
fi
echo "run_dmrg: LD_PRELOAD → ${PRELOAD:-<none>}"
echo "run_dmrg: running → solange_dmrg.py $*"
echo "──────────────────────────────────────────────────────────────────────"

# 4) run, passing all arguments straight through ───────────────────────────────
LD_PRELOAD="${PRELOAD}" python scripts/laguna/solange_dmrg.py "$@"
