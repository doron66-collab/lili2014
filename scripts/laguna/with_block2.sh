#!/usr/bin/env bash
# with_block2.sh — run ANY command with the block2 / MKL environment set up.
#
# This is the packaging fix that used to live inside run_dmrg.sh, extracted so
# more than one entry point can use it. block2's wheel ships hashed MKL
# libraries but omits libmkl_def.so.1, and preloading the MKL trio without
# libgomp gives "undefined symbol: omp_get_num_procs" — so importing block2 in a
# plain `python` fails even though the package is installed correctly. Anything
# that imports block2 needs this environment, not just the DMRG classifier:
# solange_dmrg.py, solange_hpc.py --dmrg-scf, and dmrgscf_block2.py's own
# validation all hit the identical failure without it.
#
#   bash scripts/laguna/with_block2.sh python scripts/laguna/dmrgscf_block2.py
#   bash scripts/laguna/with_block2.sh python -u scripts/laguna/solange_hpc.py --dmrg-scf ...
#
# run_dmrg.sh delegates here, so the fix has exactly one home. Set SOLANGE_ENV /
# SOLANGE_REPO to match Laguna if the defaults are wrong.
set -uo pipefail

CONDA_ENV="${SOLANGE_ENV:-base}"              # env with block2/pyblock2 + pyscf
REPO_DIR="${SOLANGE_REPO:-$HOME/lili2014}"    # repo checkout on Laguna

module load conda 2>/dev/null || true
# If block2 already imports (you are already in the right env — e.g. the shell
# shows "(base)"), use the current environment as-is. Only activate CONDA_ENV
# otherwise, so this never needs a placeholder and never breaks a working session.
if ! python -c "import block2" 2>/dev/null; then
  # shellcheck disable=SC1091
  source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate "$CONDA_ENV" 2>/dev/null || source activate "$CONDA_ENV" 2>/dev/null || {
    echo "with_block2: block2 is not importable and env '$CONDA_ENV' could not be activated." >&2
    echo "  Activate your env first (e.g. 'conda activate base'), or set SOLANGE_ENV=<name>." >&2
    exit 1; }
fi
cd "$REPO_DIR" || { echo "with_block2: repo not found: $REPO_DIR — set SOLANGE_REPO." >&2; exit 1; }

# 1) locate block2's bundled library dir (block2.libs, created by the wheel) ────
BLK_LIBS="$(python - <<'PY'
import os, glob, importlib.util
spec = importlib.util.find_spec('block2')
if not spec or not spec.origin:
    print(''); raise SystemExit
# block2 may be a single .so (origin in site-packages) OR a package
# (origin in site-packages/block2/__init__...). Normalize to site-packages.
sp = os.path.dirname(spec.origin)
if os.path.basename(sp) == 'block2':
    sp = os.path.dirname(sp)
cands = glob.glob(os.path.join(sp, 'block2.libs')) \
      + glob.glob(os.path.join(sp, 'block2*.libs')) \
      + glob.glob(os.path.join(sp, 'block2', '.libs'))
print(cands[0] if cands else '')
PY
)"
if [[ -z "$BLK_LIBS" || ! -d "$BLK_LIBS" ]]; then
  echo "with_block2: block2 not found in env '$CONDA_ENV'.  Install it:  pip install --user block2" >&2
  exit 1
fi
echo "with_block2: block2 libs → $BLK_LIBS"

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
try: roots.append(site.getusersitepackages())
except Exception: pass
for r in roots:
    hits = glob.glob(os.path.join(r, '**', 'libmkl_def.so.1'), recursive=True)
    if hits: print(hits[0]); break
PY
)"
  if [[ -n "$SRC" ]]; then
    cp -n "$SRC" "$BLK_LIBS"/ && echo "with_block2: copied libmkl_def.so.1 ← $SRC"
  else
    echo "with_block2: WARNING libmkl_def.so.1 not found. If block2 errors on it, run:" >&2
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
  echo "with_block2: WARNING could not find the MKL lp64/core libs in block2.libs;" >&2
  echo "         running without LD_PRELOAD — if it errors, that is the cause." >&2
fi
echo "with_block2: LD_PRELOAD → ${PRELOAD:-<none>}"
echo "with_block2: running → $*"
echo "──────────────────────────────────────────────────────────────────────"

# 4) run the caller's command with that environment in place ───────────────────
exec env LD_PRELOAD="${PRELOAD}" "$@"
