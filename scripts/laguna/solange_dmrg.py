#!/usr/bin/env python3
"""
solange_dmrg.py — classical DMRG classifier for SOLANGE target A/B/C tiers.

Determines, from FIRST PRINCIPLES (not size), whether a target's active space is
classically solvable to chemical accuracy — the defensible basis for A/B/C:

  C  classical exact reaches chemical accuracy   (≤ ~18 active e⁻; CCSD(T)/FCI)
  B  DMRG (classical, approximate) reaches it     (converges at practical bond dim)
  A  even DMRG cannot                             (strong correlation / high
                                                    entanglement; quantum-necessary)

DMRG is a CLASSICAL algorithm (block2). It runs on CPU (largemem — big bond
dimensions need RAM) and can be GPU-accelerated; either way it is classical, NOT
quantum. The point of THIS script is exactly to find where classical (incl. DMRG)
stops reaching chemical accuracy — because that, not the 18e exact wall, is where
quantum is truly necessary.

Method: run DMRG at increasing bond dimension M. If the energy stops changing by
more than chemical accuracy (1 kcal/mol ≈ 1.6 mHa) at a practical M, DMRG has
delivered → the target is classically tractable (B). If it is still changing, or
the entanglement entropy is high (bond dim would blow up), DMRG has NOT delivered
→ quantum-necessary (A). The max bipartite entanglement entropy S_max is reported
as the physical reason.

USAGE (validate on a small case first):
  python solange_dmrg.py --compound acetamide --basis 6-31g --ncas 8 --nelecas 8 \
      --key ARID2_LOF --out ./out

HONEST SCOPE: this classifies whatever active space you give it. Classifying a
target's FULL functional site (e.g. ARID2 56-qubit site) rigorously requires
building that active space from the PDB structure first — a separate pipeline.
On the tiny model compounds it demonstrates the diagnostic, not a full-site verdict.
"""

import argparse
import hashlib
import json
import math
import socket
import subprocess
import sys
import time
import uuid
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))          # for solange_hpc
sys.path.insert(0, str(_HERE.parents[2]))      # repo root, for generate_expansion_jw

CHEM_ACC_MHA = 1.6          # 1 kcal/mol
# The bond dimension this pipeline is willing to spend. Not derived: a stated
# commitment, and the value that gives "practical" its meaning everywhere below.
# Cost is NOT a fixed M^3 law - measured empirically on this project's own runs
# at ~M^1.91, against Zhai et al. (2026)'s own measured ~D^3.05 on their FeMoco
# benchmark (2026-10 correction: this project's largest published-chemistry
# reference point, M=6000 on thousands of cores for weeks, CAS(113,76), traces
# to Zhai & Chan (2021, JCTC) and Zhai et al. (2026) - NOT to Reiher et al.
# (2017, PNAS), which is a quantum resource-estimation paper with no DMRG run
# of its own and is cited elsewhere in this project for that, unrelated, claim.
# 2000 is still a useful scale reference regardless of exponent: 512x the M=250
# the DMRG-SCF solver uses by default, and roughly what a single node reaches
# in hours. Calibration on an N2 stretch series at CAS(10,16) found the
# requirement in the 160-256 range, an order of magnitude under this ceiling -
# so it is generous for small systems and untested against large ones.
# Raising it does not make a target classically tractable; it changes what this
# pipeline is prepared to pay before it says so.
PRACTICAL_M  = 2000
# Exact diagonalisation builds every determinant in the active space, so its
# cost is combinatorial in the size of that space and carries no dependence on
# entanglement whatever. For CAS(n,n) at singlet spin the determinant count is
# C(n, n/2)^2:  16e -> 1.7e8,  18e -> 2.4e9,  20e -> 3.4e10,  22e -> 5.0e11.
# Eighteen is where routine exact diagonalisation stops being routine — a few
# billion determinants is reachable with effort, tens of billions is not.
#
# It is a practical convention, not a derived constant, and it is NOT "the
# classical wall": DMRG carries the classical frontier far past it (FeMoco at
# CAS(113,76)). What this boundary marks is where the classifier must START
# ASKING about entanglement. Below it the question is moot, because exact
# diagonalisation settles the space regardless of how entangled the state is.
# Above it that shortcut is gone and S_max decides.
EXACT_WALL_E = 18   # active electrons up to which exact diagonalisation is routine
# Max bipartite entanglement above which DMRG stops being practical AT
# PRACTICAL_M. Calibrated, not derived: N2 at equilibrium measures S_max = 0.42
# (weakly correlated) and stretched N2 measures 2.86 (strongly correlated), and
# 1.5 sits between them. A policy of this pipeline, not a constant of nature -
# it marks where this pipeline's committed bond dimension runs out.
#
# Deliberately NOT replaced by a value derived from M ~ C*e^S. Calibration
# (results/calibration/N2_CAS10-16_ccpvdz_2026-08-05.json) measured that
# relation directly and found a growth exponent of 0.31, not 1: the exponential
# form is an upper bound, following from a flat entanglement spectrum, and real
# spectra decay. Both its prefactor and its exponent depend on the chemistry and
# on the size of the active space, so a fit on N2 does not transfer to a
# 48-electron site.
S_HARD       = 1.5


# block2 pre-allocates its own memory pool and aborts with "exceeding allowed
# memory" when a sweep outgrows it. That pool is NOT the cluster's limit: the
# default is about 1 GB, roughly an eighth of the per-user cgroup ceiling here,
# and every run before this default was set was silently capped by it rather than
# by the machine. Sized to leave generous headroom under the cgroup ceiling.
DEFAULT_STACK_MEM_GB = 4.0


def detect_hardware(n_threads: int) -> str:
    """CPU/node identity for this run, added because dmrg_classifications carried
    no hardware column at all — unlike simulation_runs' p3_backend, which records
    the GPU actually used by the 3A-HPC CASSCF/VQE path (solange_hpc.py's
    detect_gpu(), via nvidia-smi). DMRGDriver here is constructed with only
    n_threads, no GPU/CUDA argument (see run_dmrg()) — DMRG in this pipeline is
    CPU-only. A GPU may still be physically present on the node (Laguna's HPC
    nodes carry NVIDIA L40S cards for the VQE path); this reports that presence
    for honesty, explicitly labelled unused, rather than letting it be read as
    the hardware this classification actually ran on."""
    host = socket.gethostname()
    cpu_model = None
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    cpu_model = line.split(":", 1)[1].strip()
                    break
    except Exception:
        pass
    cpu_desc = cpu_model or "CPU (model unavailable)"

    gpu_note = ""
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            stderr=subprocess.DEVNULL, timeout=10).decode().strip()
        if out:
            gpu_note = (f" · GPU present ({out.splitlines()[0]}) but NOT used — "
                        f"block2's DMRGDriver here is CPU-only")
    except Exception:
        pass

    return f"{cpu_desc} · {n_threads} threads · host={host}{gpu_note}"


class OrbitalTimeBudgetExceeded(Exception):
    """Raised from the CASSCF callback when the orbital phase outruns its budget.

    --max-minutes was written to protect the bond-dimension ladder, and it does.
    What it never covered was the phase BEFORE the ladder: DMRG-SCF orbital
    optimisation, which on a real protein cluster is the expensive part. On a
    163-atom site the DMRG solves themselves took under two seconds each while
    the macro-iterations around them took nearly four minutes, because CASSCF
    re-transforms integrals over 491 basis functions every time. Twenty
    macro-iterations of that overruns a four-hour allocation, and the scheduler
    then kills the job before the ladder starts - which is to say before the
    mechanism meant to save a partial result ever runs. Exactly the failure the
    budget was added to prevent, one stage earlier than it was looking.

    Stopping this way keeps whatever orbitals the optimisation had reached. They
    are usable and they are not converged, so every record built from them is
    marked accordingly rather than presented as a finished optimisation.

    Carries e_tot (the last completed macro-iteration's energy) alongside
    imacro. mc.e_tot is never assigned during mc1step.kernel() — 'e_tot' there
    is a plain local variable, only copied onto the CASSCF object by the
    wrapper AFTER kernel() returns normally (confirmed by reading mc1step.py
    directly: every 'e_tot' in the function body is a bare local, no
    'self.e_tot' assignment exists before the final return). Interrupting the
    loop early means that copy never happens, so mc.e_tot stays at its
    object-creation default (None) regardless of how many macro-iterations
    actually completed — float(mc.e_tot) then raises TypeError, turning a
    working guard into a crash (hit live on solange_hpc.py's identical
    --compound-path copy of this class, TP53_C275F CAS(14,14), macro-iteration
    9). The energy has to come from envs['e_tot'] — the callback's own
    locals(), captured at the moment of the raise — not from the mc object
    afterward.
    """
    def __init__(self, imacro, e_tot=None, mo=None):
        super().__init__(imacro)
        self.imacro = imacro
        self.e_tot = e_tot
        # Same story as e_tot, one level deeper: mc.mo_coeff is ALSO never
        # touched inside mc1step.kernel() — 'mo' is a bare local there too,
        # never written back to self.mo_coeff. An interrupted run leaves
        # mc.mo_coeff at its PRE-OPTIMIZATION value; get_h1eff()/get_h2eff()
        # called without an explicit mo_coeff= would silently read that stale
        # value. Caught live on solange_hpc.py's --compound-path copy of this
        # exact bug by that file's active-space-FCI consistency gate
        # (RuntimeError: active-space FCI != e_casscf) — the fix is to hand
        # the CURRENT mo through explicitly, here too.
        self.mo = mo


def run_dmrg(h1e, h2e, ecore, ncas, nelecas, bond_dims, scratch="./tmp_dmrg",
             n_threads=4, max_minutes=None, early_stop=True,
             stack_mem_gb=DEFAULT_STACK_MEM_GB, n_sweeps=10,
             noises=None, tol=1e-7, spin=0):
    """Run DMRG at increasing bond dimensions. Returns per-M energies + S_max.

    spin (added 2026-10-07, for SDHB S3's odd-electron CAS(53,38) ground
    state — na-nb=1, no spin=0 sector exists): the target 2S for block2's
    DMRGDriver.initialize_system(), same convention as PySCF's mol.spin.
    Defaults to 0, so every existing caller of this function (every prior
    run — C275F etc., all closed-shell active spaces) is completely
    unaffected; only a caller that explicitly passes spin=N changes behavior.

    HPC-ticket-aware: prints live per-M timing (so `tail -f` shows real progress,
    letting you judge whether to Ctrl+C before a fixed-walltime allocation ends),
    and honors max_minutes — once the cumulative wall-clock budget is exceeded, it
    stops requesting new (larger) bond dimensions and returns whatever it has
    rather than being killed mid-sweep with nothing recorded. Reusing the same
    `scratch` directory across separate job submissions lets block2 resume the MPS
    from where a prior run left off instead of restarting from bond_dims[0].

    noises/tol (Claude Science, 2026-10-02, found live on TP53_C275F_LADDER_32_256):
    the previous hardcoded noises=[1e-5, 1e-6, 0] against n_sweeps=10, with tol never
    passed (pyblock2 default 1e-8 Ha), combined into an UNDETECTED early-stop: block2
    pads a short noises list by repeating its LAST entry for every remaining sweep
    (confirmed against pyblock2/driver/core.py's own default, [1e-5]*5 + [0], which
    is itself only 6 entries against n_sweeps=10 -- the library's own default only
    makes sense if this is how it behaves), so every M ran 8 of its 10 sweeps at
    EXACTLY ZERO noise -- the standard DMRG mechanism for escaping a bad local
    minimum right after a bond-dimension increase, switched off almost immediately
    after every increase. Verified on real sweep_history.json data: M=16 stopped at
    9/10 sweeps (tol=1e-8 default + zero noise satisfied the convergence+zero-noise
    exit condition one sweep early). The new default below keeps noise on through
    sweep 14 (vs. sweep 2 before), so the same exit condition cannot fire before
    sweep 15, and passes tol explicitly instead of relying on pyblock2's default.
    """
    if noises is None:
        noises = [1e-5] * 8 + [1e-6] * 6 + [0] * 6
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=n_threads,
                     stack_mem=int(stack_mem_gb * (1 << 30)))
    drv.initialize_system(n_sites=ncas, n_elec=nelecas, spin=spin)
    # Free sector confirmation (Claude Science, 2026-10-02): this single call sets
    # the target (N, S) sector for BOTH the MPO built below and the random MPS
    # get_random_mps() creates right after -- so they are in the same symmetry
    # sector by construction, and a <ref|H_ladder|ref> expectation value between a
    # differently-sectored pair would have failed outright rather than returned a
    # number. Printed as a record, not as a test that can fail silently.
    print(f"  [sector] DMRGDriver.initialize_system(n_sites={ncas}, n_elec={nelecas}, "
          f"spin={spin}) -- target sector shared by this run's MPO and MPS", flush=True)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=h2e, ecore=ecore, iprint=0)
    # The module docstring above has claimed since this pipeline's early days that
    # reusing --scratch resumes a killed run instead of restarting the bond-dimension
    # ladder from scratch — that claim was never actually implemented: this used to be
    # an unconditional get_random_mps() call, so every re-submission silently redid the
    # full ladder from a fresh random MPS regardless of what was already on disk. Caught
    # 2026-08-14 by comparing wall-clock times across two R175H submissions that used the
    # same --scratch and took the same time at M=250 and M=500 — proof the "resume" was
    # never happening. Try loading a previously saved MPS by tag first; only fall back to
    # a fresh random one if none exists (first submission into this scratch dir) or the
    # load fails for any reason (block2 API/version mismatch, corrupted state, etc.) —
    # a failed load must never crash the run, only cost it the resume it would have given.
    #
    # A saved MPS's own (ncas, nelecas) is checked BEFORE ever calling load_mps() —
    # found live 2026-09-07: two different targets sharing one --scratch directory
    # (a real submission mistake, not this function's fault) produced a SEGFAULT
    # inside block2's C++ load path when the dimensions didn't match, not a catchable
    # Python exception — the try/except above is no protection against that, since a
    # segfault never reaches it. A small sidecar metadata file, written once
    # per (ncas, nelecas) and checked here, means a mismatch is refused in Python
    # before block2 ever sees the incompatible file, at the cost of one extra
    # resume being treated as "start fresh" if the metadata file itself is missing
    # (e.g. a scratch dir from before this check existed).
    meta_path = Path(scratch) / ".solange_dmrg_meta.json"
    dims_match = False
    if meta_path.exists():
        try:
            saved_dims = json.loads(meta_path.read_text())
            # spin added 2026-10-08 alongside run_dmrg(spin=...): two legs of one
            # target share (ncas, nelecas) but not spin, and an MPS from one spin
            # sector must never be "resumed" into another. Metadata written before
            # this existed has no "spin" key -> read as 0, which is exactly what
            # every pre-existing caller used, so their resume behaviour is unchanged.
            dims_match = (saved_dims.get("ncas") == ncas and saved_dims.get("nelecas") == nelecas
                          and saved_dims.get("spin", 0) == spin)
            if not dims_match:
                print(f"  [resume] {scratch} holds a saved MPS for a DIFFERENT active space "
                      f"(CAS({saved_dims.get('nelecas')},{saved_dims.get('ncas')}) on disk vs. "
                      f"CAS({nelecas},{ncas}) requested) — refusing to load it (this is what a "
                      f"mismatched load segfaults on inside block2) and starting fresh instead.",
                      flush=True)
        except Exception:
            dims_match = False  # unreadable/corrupt metadata — treat like no metadata at all
    if dims_match:
        try:
            ket = drv.load_mps(tag="KET", nroots=1)
            print(f"  [resume] loaded existing MPS from {scratch} — continuing, not restarting "
                  f"the bond-dimension ladder from scratch.", flush=True)
        except Exception as e:
            ket = drv.get_random_mps(tag="KET", bond_dim=min(bond_dims[0], 250), nroots=1)
            print(f"  [resume] no usable saved MPS in {scratch} ({type(e).__name__}) — "
                  f"starting from a fresh random MPS.", flush=True)
    else:
        ket = drv.get_random_mps(tag="KET", bond_dim=min(bond_dims[0], 250), nroots=1)
        if meta_path.exists():
            pass  # already explained above (dimension mismatch)
        else:
            print(f"  [resume] no active-space metadata recorded in {scratch} — "
                  f"starting from a fresh random MPS.", flush=True)
    Path(scratch).mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps({"ncas": ncas, "nelecas": nelecas, "spin": spin}))
    energies = []
    # discarded_weights: best-effort capture of block2's per-M truncation error,
    # for dmrg_extrapolate.py's E(w)->w=0 fit (Claude Science, 2026-09-30) — a
    # sturdier convergence test than raw ΔE between bond dims, which conflates
    # truncation error with an unconverged sweep (the R175H cold-start
    # artifact). The exact attribute (drv._dmrg.discarded_weights vs.
    # .sweep_discarded_weights) was found empirically on Laguna's pyblock2
    # build via a tiny probe, NOT from documentation — see that probe's
    # findings before trusting which one this reads. Wrapped in try/except,
    # matching s_max's own defensive style below: this is instrumentation, and
    # its absence must never fail a run that would otherwise report a valid
    # energy/S_max.
    discarded_weights = []
    # Full per-sweep history, one entry per bond dimension, each holding ALL
    # sweeps' energies/discarded-weights for that M (not just the last one).
    # Added 2026-10-02 (Claude Science): a negative extrapolate() slope has two
    # distinct causes needing OPPOSITE fixes -- too few sweeps at fixed M (more
    # sweeps, same M, helps) vs. the optimizer landing in different local minima
    # at different M (more sweeps do nothing; warm-start helps) -- and nothing
    # short of the per-sweep trace within each M can tell them apart. Keeping
    # only the last sweep's dw (as before) throws away exactly the data
    # dmrg_extrapolate.sweep_convergence() needs to distinguish the two cases
    # for free, after the fact, instead of guessing which re-run to try next.
    sweep_history = {}
    stop_reason = "completed"            # completed | converged | time_budget
    t_start = time.time()
    for M in bond_dims:
        elapsed_before = time.time() - t_start
        if max_minutes is not None and elapsed_before > max_minutes * 60 and energies:
            print(f"  [time budget] {elapsed_before/60:.1f}m elapsed > --max-minutes "
                  f"{max_minutes} — stopping before M={M}; returning {len(energies)} "
                  f"completed bond dim(s). Re-run with the same --scratch to resume.",
                  file=sys.stderr)
            stop_reason = "time_budget"
            break
        t0 = time.time()
        e = drv.dmrg(mpo, ket, n_sweeps=n_sweeps, bond_dims=[M],
                     noises=noises, thrds=[1e-9] * 3, tol=tol, iprint=0)
        dt = time.time() - t0
        energies.append((M, float(e)))
        try:
            # get_dmrg_results() is the DOCUMENTED accessor (pyblock2/driver/
            # core.py): returns (bond_dims, dws, energies), all three PER SWEEP,
            # with dws "the maximal discarded weight (sum of discarded
            # eigenvalues) for each sweep". Settled 2026-09-30 after two guessed
            # attributes (drv._dmrg.discarded_weights / .sweep_discarded_weights,
            # read directly with no accessor) both failed a real H8-vs-exact-FCI
            # validation -- discarded_weights turned out to BE the right field
            # (per-sweep, as guessed) but sweep_discarded_weights is untouched by
            # the Python driver entirely (a C++-internal field the wrapper never
            # populates) and should never have been read. dws[-1] is this
            # dmrg() call's LAST sweep (matches the energy recorded above, since
            # this call covers exactly one M) -- see dmrg_extrapolate.
            # last_sweep_per_bond_dim's docstring for why last-sweep, not max.
            _, dws, sweep_energies = drv.get_dmrg_results()
            dw = float(dws[-1])
            sweep_history[M] = {"dws": [float(x) for x in dws],
                                "energies": [float(x[0]) if hasattr(x, "__len__") else float(x)
                                            for x in sweep_energies]}
        except Exception:
            dw = None
        discarded_weights.append(dw)
        # stdout (not stderr) so the agent's live-progress streamer reliably captures
        # each "DMRG M=" line and surfaces it as a bond-dimension stage note.
        # n_done/n_sweeps printed explicitly (Claude Science, 2026-10-02): block2's
        # own tol+zero-noise early-exit condition can stop a given M short with no
        # other visible signal -- found live when M=16 (of this same ladder, before
        # this fix) silently ran only 9 of its 10 requested sweeps. Without this
        # count next to the energy, "converged" and "exited early" print identically.
        n_done = len(sweep_history[M]["dws"]) if M in sweep_history else "?"
        print(f"  DMRG M={M:5d}  E={float(e):.8f} Ha  [{n_done}/{n_sweeps} sweeps run, "
              f"{dt:.1f}s this M, {(time.time()-t_start)/60:.1f}m total]", flush=True)
        # Early stop: once the energy stops improving by more than chemical accuracy
        # between consecutive bond dims, larger M cannot change the verdict — the
        # answer has converged. This is exactly classify()'s own convergence test,
        # so stopping here loses no evidence (the last two points already agree),
        # and it saves the expensive high-M sweeps on easy/weakly-correlated cases.
        # A hard, strongly-correlated (Class-A) system keeps dropping, so it never
        # triggers and runs the full sweep — precisely where the evidence matters.
        if early_stop and len(energies) >= 2:
            dE = abs(energies[-1][1] - energies[-2][1]) * 1000.0     # mHa
            if dE < CHEM_ACC_MHA:
                print(f"  [early stop] ΔE(M={energies[-2][0]}→{M})={dE:.3f} mHa "
                      f"< {CHEM_ACC_MHA} (chemical accuracy) — converged; skipping "
                      f"larger bond dims ({len(bond_dims)-len(energies)} remaining).",
                      file=sys.stderr)
                stop_reason = "converged"
                break
    try:
        s_max = float(np.max(drv.get_bipartite_entanglement(ket)))
    except Exception:
        s_max = None
    if sweep_history:
        hist_path = Path(scratch) / "sweep_history.json"
        hist_path.parent.mkdir(parents=True, exist_ok=True)
        with open(hist_path, "w") as fh:
            json.dump(sweep_history, fh)
        print(f"  [sweep-history] wrote {hist_path} (per-sweep E/dw for every M on this "
              f"ladder -- feed to dmrg_extrapolate.sweep_convergence() to diagnose a "
              f"negative extrapolation slope without re-running)", flush=True)
    return energies, s_max, stop_reason, discarded_weights, sweep_history


# S_max is a bipartite entanglement entropy. CORRECTED 2026-09-29 (Claude
# Science, after reading the full dissertation): the naive ceiling S<=k*ln(4)
# (k = smaller side's orbital count) ignores particle-number conservation and
# is WRONG -- the real governing variable is FILLING (n_elec/(2*n_orb)), not
# orbital count alone. The exact combinatorial ceiling for a bipartition into
# a kL/kR orbital split is S <= ln(sum over (n_Lup,n_Ldown) sectors of
# min(d_L, d_R)), where d_L/d_R are the electron-sector Hilbert-space
# dimensions on each side (d = C(k,n_up)*C(k,n_down)) -- summed only over
# sectors consistent with the fixed total (nelec_alpha, nelec_beta). This
# REVERSES the previous naive-formula reading for the two cases that matter
# most: R175H (CAS(120,65), S_max=1.31) was 3.0% of the naive ceiling --
# below NEGCTRL_BORING's 4.1%, read as noise. Under the exact formula R175H
# is 8.2% -- ABOVE NEGCTRL_BORING's 4.2%: the first positive biological
# signal in the project, not noise (R175H has far FEWER orbitals than
# NEGCTRL_BORING but is far more densely FILLED -- 92.3% vs 60% -- which the
# naive k-only formula could not see). Uses log-space (lgamma) throughout so
# it never overflows float on the largest active spaces (CAS(338,235)).
# Reported as a diagnostic alongside S_max, not used to change any
# classification decision here -- see classify()'s own comment for why the
# trigger logic itself is being left alone for now (dissertation §06.i's
# "minimal path" fix: report both signals plainly, downgrade the CLAIM
# S_HARD supports, without re-running or re-deciding anything).
def _log_comb(n, k):
    if k < 0 or k > n:
        return None
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def entanglement_capacity_pct(s_max, ncas, nelecas=None, spin=0):
    if s_max is None or not ncas:
        return None
    if nelecas is None:
        # Fallback for callers that don't have electron counts on hand: the
        # old naive k*ln4 ceiling, clearly WRONG (see comment above) but
        # better than nothing and never silently claimed as the exact one.
        k = max(ncas // 2, 1)
        return 100.0 * s_max / (k * math.log(4))

    nelec_alpha = (nelecas + spin) // 2
    nelec_beta = (nelecas - spin) // 2
    kL = ncas // 2
    kR = ncas - kL
    log_terms = []
    for n_l_up in range(max(0, nelec_alpha - kR), min(kL, nelec_alpha) + 1):
        for n_l_dn in range(max(0, nelec_beta - kR), min(kL, nelec_beta) + 1):
            n_r_up = nelec_alpha - n_l_up
            n_r_dn = nelec_beta - n_l_dn
            if not (0 <= n_r_up <= kR and 0 <= n_r_dn <= kR):
                continue
            log_d_l = _log_comb(kL, n_l_up) + _log_comb(kL, n_l_dn)
            log_d_r = _log_comb(kR, n_r_up) + _log_comb(kR, n_r_dn)
            log_terms.append(min(log_d_l, log_d_r))
    if not log_terms:
        return None
    m = max(log_terms)
    log_capacity = m + math.log(sum(math.exp(t - m) for t in log_terms))
    return 100.0 * s_max / log_capacity


# METHODOLOGICAL LIMITATION, recorded rather than fixed here: this test is
# single-method. A high S_max means DMRG (an MPS/tensor-network method)
# specifically cannot represent the state at a practical bond dimension. It
# does NOT by itself rule out other classical families with different
# sensitivity to entanglement — selected CI, FCIQMC, non-MPS tensor-network
# geometries, neural quantum states. The dissertation's own prose already
# states the classification is "conditional on the current state of classical
# methods... and revisable as those methods advance" (§06.i) — this comment
# exists because that hedge was in the prose but not reflected anywhere in
# the code that actually decides a class. A Class A verdict here is accurate
# to what it claims to be: DMRG-specific evidence, not a cross-method
# exclusion. Treated as a genuinely open question, not one this platform
# currently answers.
# Fraction of PRACTICAL_M a ladder must reach before non-convergence alone is
# allowed to support Class A (Claude Science, 2026-10-03, found live on the
# corrected TP53_C275F ladder): S_max is an INTRINSIC property of the state —
# true at any M once measured — but "DMRG hasn't converged yet" is a PROCEDURAL
# fact about how far the ladder got, not about the molecule. The ladder that
# triggered this stopped at M=250, 6.2% of PRACTICAL_M=2000... er, of this
# project's own feasibility ceiling — "we haven't finished the classical
# computation" is not evidence FOR quantum necessity, it is an open question.
# OR-ing an intrinsic criterion with an incomplete procedural one let an
# unfinished run report a positive (Class A) verdict, which is the wrong
# direction to default an incomplete answer. This fraction is a judgment call,
# not a derived number: chosen conservatively (DMRG should be well into its
# practical range, not just past chemical-accuracy's own tiny multiple) and
# open to revision with more calibration runs.
NONCONVERGENCE_FEASIBILITY_FRACTION = 0.5


def classify(active_electrons, energies, s_max):
    """Map DMRG behaviour to an A/B/C/INCONCLUSIVE class with an explicit rationale.

    Two signals:
      • ΔE across the last two bond dims — direct evidence DMRG has (not) converged
        at a practical M. Only bites at large active spaces; small ones are exact.
      • S_max (max bipartite entanglement) — predicts the bond dimension the FULL
        site needs (M ~ e^S). This is the leading indicator, since a small model
        space is always DMRG-exact yet still reveals the correlation strength.

    INCONCLUSIVE (new, 2026-10-03): non-convergence alone, far short of this
    project's own feasibility ceiling (PRACTICAL_M), must not produce Class A —
    "quantum-necessary" is a claim about the molecule, and a ladder that simply
    stopped early has not earned it. S_max > S_HARD still produces Class A at any
    M, since that signal does not depend on how far the ladder climbed.
    """
    dE_final = abs(energies[-1][1] - energies[-2][1]) * 1000 if len(energies) >= 2 else None
    m_reached = energies[-1][0] if energies else None
    converged = (dE_final is not None) and (dE_final < CHEM_ACC_MHA) \
                and (m_reached is not None) and (m_reached <= PRACTICAL_M)
    strong = (s_max is not None) and (s_max > S_HARD)
    smx = "n/a" if s_max is None else f"{s_max:.2f}"
    near_feasibility_ceiling = (m_reached is not None) and \
        (m_reached >= NONCONVERGENCE_FEASIBILITY_FRACTION * PRACTICAL_M)

    if active_electrons <= EXACT_WALL_E:
        return "C", (f"{active_electrons}e ≤ {EXACT_WALL_E}e exact-classical wall — "
                     f"CCSD(T)/FCI reaches chemical accuracy; classical sufficient. "
                     f"(S_max={smx})")
    if converged and not strong:
        return "B", (f"DMRG reaches chemical accuracy at practical M={energies[-1][0]} "
                     f"(ΔE={dE_final:.2f} mHa) and entanglement is low (S_max={smx} < "
                     f"{S_HARD}) — classical (DMRG) delivers; quantum-advantaged, not necessary.")
    if strong:
        return "A", (f"provisional quantum-necessary — S_max={smx} > {S_HARD} (strong correlation; "
                     f"DMRG bond dim ~e^S blows up at the full {active_electrons}e site). Unresolved "
                     f"classically within this run's declared budget; not a proof of necessity, "
                     f"which no finite computation can establish.")
    if not converged:
        return "INCONCLUSIVE", (
            f"DMRG not at chemical accuracy by practical M={m_reached} (ΔE="
            f"{'n/a' if dE_final is None else round(dE_final, 2)} mHa), and S_max={smx} "
            f"is below {S_HARD} so entanglement alone does not support Class A. "
            + (f"M={m_reached} is only {100*m_reached/PRACTICAL_M:.0f}% of this project's "
               f"own feasibility ceiling (PRACTICAL_M={PRACTICAL_M}) — this is an unfinished "
               f"classical computation, not evidence of quantum necessity. Climb the ladder "
               f"before classifying."
               if not near_feasibility_ceiling else
               f"M={m_reached} is {100*m_reached/PRACTICAL_M:.0f}% of PRACTICAL_M="
               f"{PRACTICAL_M} — close to this project's feasibility ceiling without "
               f"converging; Class A is defensible but not yet declared automatically. "
               f"Review before classifying."))
    # Should be unreachable (converged-and-strong falls through to here only if
    # converged but strong, i.e. chemical accuracy reached yet S_max > S_HARD --
    # DMRG agrees with itself at a practical M while reporting high entanglement).
    return "A", (f"provisional quantum-necessary — S_max={smx} > {S_HARD} despite DMRG convergence at "
                 f"M={m_reached}. A direction within the declared budget, not a proof of necessity.")


def integrals_from_geometry(xyz_path, basis, avas_aos, charge=0, spin=0, verbose=0,
                             density_fit=True, max_memory=16000, df_auxbasis="def2-universal-jkfit",
                             max_cycle_macro=20, dmrg_scf=False, dmrg_scf_maxm=500,
                             dmrg_scf_scratch="./tmp_dmrgscf_orb", n_threads=4,
                             orbital_deadline=None, stack_mem_gb=None, casci=False,
                             avas_threshold=0.2, load_orbitals=None,
                             ncas_override=None, nelecas_override=None,
                             ncore_override=None,
                             skip_precondition=False):
    """Chemist-in-the-loop entry: given a QM-cluster geometry (xyz) and the target
    atomic orbitals, AVAS selects the active space automatically. Returns a dict
    shaped like run_casscf's output. The CLUSTER itself (which residues/atoms/metal,
    H-capping) is the chemist's input — that step is NOT auto-generated here.

    avas_threshold defaults to pyscf's own library default (0.2), never before
    exposed here — every prior --geometry run used it implicitly, silently,
    with no record that a choice was even made. Widening it (e.g. 0.1, 0.05)
    admits more virtual orbitals into the active space and is how three
    targets (TP53 R282W, KEAP1 G333C, TP53 C275F) were previously pushed from
    Class B to Class A — but that was done via an undocumented, unrecorded
    code path outside this script, not reproducibly through this CLI. Exposed
    now specifically so a widened-threshold run (and its negative-control
    counterpart on a chemically boring cluster of matched size) goes through
    the same auditable script and bond-dimension protocol as every other run,
    instead of an ad hoc one-off.

    charge/spin are ALSO the chemist's input, not guessed: pyscf defaults to a
    neutral (charge=0), closed-shell (spin=0) molecule, which silently assumes
    the atom count sums to an even number of electrons. A QM cluster cut from a
    protein (capped side chains, possibly missing hydrogens on the source
    structure) essentially never satisfies that by accident — get an odd
    electron count with the defaults and pyscf raises "Electron number N and
    spin 0 are not consistent" rather than guessing wrong silently.

    dmrg_scf: replace CASSCF's internal solver — ordinarily FCI, exact but
    combinatorial in ncas, impractical past ~16 active orbitals (see the
    size_note below) — with DMRG (block2) itself. This is DMRG-SCF, not the
    two-stage "CASSCF/FCI picks orbitals, DMRG only refines the final energy on
    whatever it was handed" this module otherwise does: here DMRG is inside the
    orbital-optimization loop, so the WHOLE procedure — orbital selection AND
    the active-space solve — scales polynomially in ncas, not exponentially.
    That is what actually lets --ncas grow past 16 toward a real ~48-155e site;
    plain CASSCF (dmrg_scf=False) cannot, no matter how the bond dims used
    downstream in run_dmrg() are set, because it never gets a converged active
    space at that size in the first place.

    dmrg_scf_maxm is deliberately a SEPARATE, usually much smaller bond
    dimension than --bond-dims: it only needs to be large enough for stable
    orbital gradients during optimization, not for the final converged energy
    (run_dmrg(), called afterward in main() on the resulting fixed integrals,
    is what sweeps up through the real --bond-dims list for the reported
    energy/S_max).

    The solver itself is dmrgscf_block2.Block2FCISolver — our own adapter over
    pyblock2, written because the usual pyscf-dmrgscf package shells out to a
    compiled `block2main` executable that this environment does not have (only
    the Python API is installed; confirmed empirically on Laguna). Its two
    verification checks are described in that module's docstring, and
    validate() is called before every optimization here: a 2-RDM convention
    error raises rather than converging silently to a wrong energy."""
    from pyscf import gto, scf, ao2mo, fci
    from pyscf.mcscf import avas
    geom = Path(xyz_path).read_text()
    # accept a raw xyz (skip the 2 header lines if present)
    lines = [l for l in geom.splitlines() if l.strip()]
    if lines and lines[0].strip().isdigit():
        lines = lines[2:]
    # pyscf's own default (4000 MB) is conservative for a desktop, not a Laguna
    # largemem node (1.5 TB) — a real protein cluster's CASSCF/FCI step routinely
    # wants more, and the alternative is a "needs N MB, over max_memory limit"
    # warning (or a hard failure) instead of just... using the RAM that's there.
    mol = gto.M(atom="\n".join(lines), basis=basis, charge=charge, spin=spin,
                verbose=verbose, max_memory=max_memory)
    # RHF requires closed-shell (spin=0); anything else needs ROHF instead. A QM
    # cluster this size (hundreds of atoms once real protonation is included,
    # easily 1000+ basis functions) makes conventional (non-density-fitted) SCF
    # AND CASSCF a real bottleneck — CASSCF's orbital optimization still touches
    # every basis function each macro-iteration, not just the active space, so
    # disabling DF doesn't just slow RHF, it makes CASSCF itself take hours even
    # for a tiny active space (learned the hard way: CAS(18,12) — trivial for
    # the FCI solver — sat for 2+ hours in a non-DF CASSCF on a 491-function
    # basis). Density fitting propagates from mf into mcscf.CASSCF(mf, ...)
    # automatically and keeps CASSCF fast too — the earlier OOM wasn't caused
    # by DF itself, it was pyscf silently falling back to an oversized
    # even-tempered auxiliary basis because sto-3g has no matching JKFIT
    # companion in the standard libraries. Naming a modest, general-purpose
    # auxbasis explicitly (default def2-universal-jkfit) fixes that mismatch
    # directly instead of giving up DF's speed for both RHF and CASSCF.
    # --no-density-fit (density_fit=False) is kept as an escape hatch only.
    print(f"  running RHF/ROHF{f'  (density-fitted, auxbasis={df_auxbasis})' if density_fit else ''} on "
          f"{mol.natm} atoms, {mol.nao} basis functions (this and the AVAS step "
          f"below print nothing further until they finish unless --verbose is "
          f"raised)...", flush=True)
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    if density_fit:
        mf = mf.density_fit(auxbasis=df_auxbasis)
    mf = mf.run()
    # Free provenance check (Claude Science, 2026-10-02): a 2.9 mHa discrepancy was
    # found between this run's own E_SCF and an "E_wt" figure reported elsewhere for
    # what was assumed to be the same cluster -- a gap 1.8x the project's own
    # chemical-accuracy threshold, so it blocks any final classification even though
    # it is 164x too small to explain the separate optimizer failure the embedding
    # gate catches. Printed here, for free, instead of guessed at: which SCF energy
    # this run actually is, and a cheap fingerprint of exactly what geometry/basis/
    # charge/spin produced it, so a mismatch against another reported number is a
    # hash comparison, not a re-run.
    geom_bytes = mol.atom_coords().tobytes() + repr(mol._basis).encode()
    print(f"  [provenance] E_SCF={mf.e_tot:.8f} Ha  converged={mf.converged}  "
          f"natm={mol.natm}  basis={basis}  charge={charge}  spin={spin}  "
          f"geom_sha256={hashlib.sha256(geom_bytes).hexdigest()[:16]}", flush=True)
    # ncas_override/nelecas_override (2026-09-30): when the orbitals being loaded
    # did NOT come from this function's own avas.avas(mf, ..., threshold=...) call
    # -- e.g. avas_by_count.py's explicit (n_occ, n_vir) selection, built because
    # AVAS's single threshold is a broken cut across the occupied/virtual sigma
    # distributions (see avas_by_count.py's docstring) -- calling avas.avas() here
    # anyway would silently compute a DIFFERENT (ncas, nelec) than the loaded
    # mo_coeff's actual shape, and CASCI would either crash on a shape mismatch or,
    # worse, run on the wrong active-space size against orbitals for a different
    # one. Skip AVAS entirely in that case; the caller is asserting these exact
    # numbers match the file.
    if ncas_override is not None and nelecas_override is not None:
        ncas, nelec = ncas_override, nelecas_override
        mo = None
        print(f"  --ncas/--nelecas given with --load-orbitals: trusting the caller's "
              f"CAS({nelec},{ncas}) instead of re-deriving it from AVAS's own threshold "
              f"sweep (this run's orbitals were selected by count, not by AVAS threshold "
              f"-- see avas_by_count.py).", flush=True)
    else:
        ncas, nelec, mo = avas.avas(mf, [s.strip() for s in avas_aos.split(",")], threshold=avas_threshold)
    # --load-orbitals: substitute a PREVIOUSLY-computed mo_coeff (saved by an
    # earlier run of this same script, --scratch/mo_coeff_final.npy) for AVAS's
    # own raw output, and skip the optimization loop entirely — added to let a
    # resubmission after an infrastructure failure (OOM, Slurm --time kill) reuse
    # orbitals a truncated run already spent real compute reaching, instead of
    # redoing the whole (multi-hour) orbital-optimization phase from scratch.
    # AVAS still runs above (cheap) so ncas/nelec are the SAME deterministic
    # values the original run used — only the orbitals themselves are replaced —
    # UNLESS ncas_override/nelecas_override were given (see above), in which case
    # AVAS was skipped entirely and ncas/nelec are the caller's asserted values.
    # This is NOT the "casci" flag's meaning below (AVAS's own untouched
    # output, trivially "converged" since nothing was attempted): these
    # orbitals came from a run that was ALREADY mid-optimization when it was
    # interrupted, so they must never be reported as converged.
    if load_orbitals:
        mo = np.load(load_orbitals)
        print(f"  --load-orbitals: loaded {load_orbitals} in place of AVAS's own output — "
              f"skipping orbital optimization entirely (one fixed diagonalization on these "
              f"orbitals, which are NOT re-verified as converged).", flush=True)
    # Report the active space BEFORE paying for CASSCF, not after: CASSCF cost is
    # roughly combinatorial in ncas, so a caller needs this number while they can
    # still Ctrl+C and narrow --avas, not only once the (possibly hours-long) run
    # has already finished or is still silently grinding.
    size_note = ""
    if ncas > 16:
        size_note = (f" *** LARGE for DMRG-SCF orbital optimization (maxM={dmrg_scf_maxm}) — "
                     f"expect it to be slow, not impossible ***" if dmrg_scf else
                     " *** LARGE for CASSCF/FCI (>16 active orbitals is often impractical "
                     "— consider narrowing --avas, or pass --dmrg-scf) ***")
    print(f"  AVAS selected active space: CAS({nelec},{ncas}){size_note}", flush=True)
    # Blocking precondition (2026-09-30): refuse to spend CASSCF/DMRG-SCF compute
    # on an active space whose verdict is already determined by its dimensions
    # alone. TP53_C275F cluster76 ran CAS(76,38) — occ=38, vir=0, exactly one
    # Slater determinant — to conclusion (~20 min DMRG-SCF each side) before
    # anyone checked that S_max=0.0 was mathematically guaranteed, not measured.
    # check_active_space.py's rules (R1 zero virtuals / R2 virtual fraction floor
    # / R3 entropy ceiling below S_HARD / R4 exactly-diagonalisable) are exactly
    # the check that data already on hand (ncas, nelecas) would have caught in
    # under a second. --skip-precondition exists only for a deliberately-C-class
    # active space (e.g. a demo/calibration run) where the "verdict is moot, not
    # missing" — do not use it to push past R1/R3 on a real classification run.
    if not skip_precondition:
        sys.path.insert(0, str(_HERE.parent))
        import check_active_space as cas_check
        _p = argparse.Namespace(s_hard=S_HARD, min_vir_frac=15.0, fci_floor=1e9, fci_ceiling=1e40)
        _fails = cas_check.report(nelec, ncas, _p)
        if _fails:
            sys.exit(f"\n*** REFUSING to run CASSCF/DMRG-SCF on CAS({nelec},{ncas}): "
                     f"{'; '.join(_fails)}. This run's classification is determined "
                     f"before it starts. Redesign the active space (see avas_by_count.py), "
                     f"or pass --skip-precondition if this is intentionally a boundary/"
                     f"calibration case, not a classification being reported. ***")
    from pyscf import mcscf
    # one_shot: structurally the SAME single-diagonalization path as --casci
    # (a CASCI object, no macro-iteration loop) — --load-orbitals forces it too,
    # since there is nothing left to optimize once the orbitals are already
    # fixed and loaded from disk. Kept as a separate variable from `casci`
    # itself because the two must NOT share casci's "trivially converged, since
    # AVAS's own untouched output was never asked to optimize" reporting below —
    # loaded orbitals came from a genuinely truncated prior run.
    one_shot = casci or bool(load_orbitals)
    if one_shot:
        # CASCI: freeze AVAS's orbitals as-is, do ONE diagonalization on them —
        # no macro-iteration loop at all. Added specifically because DMRG-SCF
        # (mcscf.CASSCF + Block2FCISolver) crashed block2 itself at CAS(48,28)
        # on the SECOND solver call (the first warm-started one, per
        # dmrgscf_block2.py's kernel() docstring) — a pattern this sidesteps
        # entirely by construction, since CASCI's kernel() calls the fcisolver
        # exactly once, cold, never warm-started. Trade-off: the orbitals are
        # AVAS's raw output, not iteratively improved against the active-space
        # energy — a real approximation, not a free win. Downstream (run_dmrg(),
        # classify()) is unaffected: it consumes the same h1e/h2e/ecore shape
        # regardless of which path produced them.
        if casci and not load_orbitals:
            print(f"  --casci: skipping orbital optimization, one fixed diagonalization "
                  f"on AVAS's orbitals as given.", flush=True)
        mc = mcscf.CASCI(mf, ncas, nelec)
    else:
        mc = mcscf.CASSCF(mf, ncas, nelec)
    if ncore_override is not None:
        # pyscf's default ncore = (mol.nelectron - nelecas) // 2 assumes the
        # active space sits immediately below the Fermi level (the standard
        # "top N occupied, bottom M virtual" convention) -- true for AVAS's
        # own raw output and for a prior solange_dmrg.py run's mo_coeff_final,
        # but FALSE for avas_mp2_select.py's output: its build_mo() places the
        # active block at an arbitrary column (active_start_col in its JSON),
        # after DROPPING every AVAS-pool candidate that correlation-ranking
        # did not select -- the saved file is not "all MOs, standard order."
        # Found live 2026-10-02 via a no-solver fingerprint check that showed
        # IDENTICAL integrals for genuinely different (canonical vs. split-PM
        # localized) orbital files on CAS(36,34): the default ncore=329 landed
        # entirely inside the UNTOUCHED virtual-beyond-AVAS-pool block (columns
        # 140:631), not the real active columns (106:140) -- so every prior
        # --load-orbitals run paired with avas_mp2_select.py's output on this
        # target computed on the WRONG 34 columns, not the chemist-selected,
        # correlation-ranked active space at all. Must be set explicitly
        # whenever the source is avas_mp2_select.py/avas_by_count.py, not left
        # to the default formula.
        print(f"  --ncore override: pyscf's default ncore would be "
              f"{(mf.mol.nelectron - nelec) // 2} (standard HOMO-adjacent convention); "
              f"using {ncore_override} instead (the real active_start_col from the "
              f"orbital-selection JSON) -- these differ whenever the source is "
              f"avas_mp2_select.py/avas_by_count.py, and using the wrong one silently "
              f"runs on the wrong 34 columns.", flush=True)
        mc.ncore = ncore_override
    if dmrg_scf:
        # Swap the orbital-optimization solver itself from FCI to DMRG (block2) —
        # see the dmrg_scf docstring above for why this, not just a bigger
        # --bond-dims list downstream, is what actually raises the ncas ceiling.
        # Uses our own pyblock2 adapter (dmrgscf_block2.py), not the
        # pyscf-dmrgscf package — that one needs a compiled `block2main` binary
        # this environment does not have. validate() cross-checks the adapter
        # against exact FCI on a small system and raises if it disagrees.
        from dmrgscf_block2 import Block2FCISolver, validate
        validate(scratch=dmrg_scf_scratch + "_validate", n_threads=n_threads,
                 verbose=True)
        # stack_mem_gb was previously NOT threaded through here, so this solver
        # silently used dmrgscf_block2's own DEFAULT_STACK_MEM_GB (4.0) no matter
        # what --stack-mem-gb was set to on the CLI -- that flag only ever reached
        # the FINAL run_dmrg() sweep below, not this orbital-optimization solver.
        # Invisible at small active spaces; at CAS(48,28) it crashed block2 itself
        # (std::length_error: cannot create std::vector larger than max_size()) on
        # the second solve, which is a memory-pool exhaustion, not a wiring bug.
        fci_kwargs = {"maxM": dmrg_scf_maxm, "scratch": dmrg_scf_scratch,
                      "n_threads": n_threads}
        if stack_mem_gb is not None:
            fci_kwargs["stack_mem_gb"] = stack_mem_gb
        mc.fcisolver = Block2FCISolver(**fci_kwargs)
        mc.internal_rotation = True  # DMRG-SCF needs this pyscf CASSCF option on
    if density_fit:
        mc = mc.density_fit(auxbasis=df_auxbasis)  # keep CASSCF's own integral transform DF-accelerated too
    mc.verbose = verbose
    mc.max_memory = max_memory  # explicit — don't rely on inheriting mol's setting
    # A hard ceiling, not a diagnosis: pyscf's own default (50) has no time bound.
    # This guarantees the run finishes within a predictable number of iterations
    # even if convergence is slow, at the cost of possibly stopping before full
    # convergence — mc.converged (checked below) reports whether that happened.
    mc.max_cycle_macro = max_cycle_macro
    if dmrg_scf:
        # Match run_casscf()'s DMRG-SCF settings — see the comment there. pyscf's
        # defaults are tighter than DMRG's own RDM truncation floor, so they can
        # never be met and the run just burns macro-iterations.
        mc.conv_tol      = 1e-6
        mc.conv_tol_grad = 1e-3
    if spin == 0 and not dmrg_scf and not one_shot:
        # fix_spin_ is a CASSCF orbital-optimization mixin method; CASCI has no
        # such loop to fix spin against, so this is skipped there rather than
        # risk an AttributeError on a class that may not define it.
        mc.fix_spin_(ss=0)  # only force the singlet when the input itself is closed-shell — FCI-solver option, not exposed by DMRGCI
    # Wall-clock guard on the orbital phase. Checked once per macro-iteration,
    # which is the only point where stopping leaves a coherent set of orbitals.
    truncated = False
    truncated_mo = None
    if orbital_deadline is not None:
        def _budget_callback(envs, _dl=orbital_deadline):
            if time.time() > _dl:
                raise OrbitalTimeBudgetExceeded(envs.get("imacro"), envs.get("e_tot"),
                                                envs.get("mo"))
        mc.callback = _budget_callback
    # casscf_converged tracks orbital-optimization convergence specifically -
    # separate from `truncated` (which only fires on the TIME-budget path).
    # Found live 2026-09-12: a run that hit max_cycle_macro without either
    # exception or mc.converged=True fell through as if nothing were wrong -
    # orbital_optimization_truncated stayed False, so main() never printed
    # the PROVISIONAL warning and the eventual classification was labeled
    # "FINAL result, not provisional" purely from DMRG's own bond-dimension
    # convergence, with no check at all on whether the ORBITALS underneath
    # it had actually settled. The same TP53_C275_NATIVE run's own log
    # printed "CASSCF did NOT converge" one line above that - the exact
    # oscillating-orbital-partition instability seen on NEGCTRL_BORING - and
    # nothing downstream acted on it. CASCI mode has no orbital optimization
    # to (not) converge, so it is always reported converged — EXCEPT
    # --load-orbitals, whose whole point is reusing orbitals a PRIOR run had
    # not yet finished optimizing when it was interrupted; reporting those as
    # "converged" via casci's normal convention would misreport exactly the
    # PROVISIONAL status this script's own orbital_optimization_converged flag
    # exists to catch (see the 2026-09-21/22 fixes elsewhere in this file).
    # external_orbitals: --load-orbitals used together with ncas_override/
    # nelecas_override (avas_by_count.py / avas_mp2_select.py's output) — a
    # DELIBERATELY built, fixed active space, not a truncated run's leftover
    # state. Conflating the two mislabels the record: Claude Science caught
    # this live 2026-09-30 reviewing the CAS(64,48) plan — "orbitals are
    # externally supplied" got reported as "NOT converged — loaded from a
    # prior truncated run", which is wrong provenance (the computation is
    # right, the record of WHY is not — exactly DP5's territory). Treated the
    # same way --casci's own orbitals are: nothing was asked to optimize, so
    # there is nothing to have failed to converge.
    external_orbitals = bool(load_orbitals) and ncas_override is not None and nelecas_override is not None
    try:
        e_casscf = mc.kernel(mo)[0]
        casscf_converged = (True if (casci or external_orbitals) else
                            False if load_orbitals else bool(mc.converged))
        conv_note = ("converged" if casscf_converged else
                     "NOT converged — loaded from a prior truncated run" if load_orbitals else
                     f"did NOT converge — hit max_cycle_macro={max_cycle_macro} first")
        if external_orbitals:
            conv_note = "orbitals externally supplied (fixed active-space selection) — optimization not applicable"
    except OrbitalTimeBudgetExceeded as exc:
        truncated = True
        casscf_converged = False
        # NOT mc.e_tot — see the class docstring: never populated on an
        # interrupted run, so float(mc.e_tot) raises TypeError on None instead
        # of the guard doing its job.
        e_casscf = float(exc.e_tot)
        truncated_mo = exc.mo   # likewise: mc.mo_coeff is stale (pre-optimization)
        conv_note = (f"STOPPED at the orbital time budget after macro-iteration {exc.imacro} "
                     f"— orbitals are usable but NOT converged")
    one_shot_label = ("CASCI (externally supplied, fixed active-space orbitals)" if external_orbitals else
                       "CASCI (externally-loaded orbitals)" if load_orbitals else
                       "CASCI (fixed AVAS orbitals)" if casci else "CASSCF")
    print(f"  {one_shot_label} {conv_note} · E={e_casscf:.8f}", flush=True)
    # Mandatory embedding gate (2026-10-02, Claude Science): a split (occ-among-
    # occ, vir-among-vir) active space always CONTAINS the HF determinant, so
    # E_CASCI/E_CASSCF must be <= E_SCF -- strictly. An energy ABOVE the SCF
    # reference means the active-space columns were sliced at the wrong offset,
    # not that the calculation is merely inaccurate. Caught live on this exact
    # target: E_CASCI = -4709.968 against E_SCF = -5333.125, +623 Ha above the
    # reference, traced to --active-start-col not matching pyscf's own
    # ncore=(mol.nelectron-nelecas)//2 arithmetic -- a comparison this cheap
    # would have ended that investigation in seconds instead of a day. Runs
    # unconditionally, not just for --load-orbitals, since the same failure
    # shape (active columns not where the energy bookkeeping assumes) is not
    # specific to that path.
    import basis_check as _bc
    _embed = _bc.check_embedding(e_casscf, float(mf.e_tot), ncore=mc.ncore,
                                 nelectron=mf.mol.nelectron, nelecas=nelec)
    for _n in _embed["notes"]:
        print(f"  [embedding] {_n}", flush=True)
    if not _embed["ok"]:
        for _f in _embed["failures"]:
            print(f"  [embedding] FAIL: {_f}", flush=True)
        sys.exit(f"\n*** REFUSING: embedding check failed on CAS({nelec},{ncas}) -- "
                 f"{_embed['verdict']}. Do not interpret this run; fix the active-space "
                 f"column offset / ncore before re-running. ***")
    # Explicit mo_coeff=truncated_mo when cut short — see OrbitalTimeBudgetExceeded's
    # docstring; None (the normal-completion case) is identical to omitting the
    # argument, since get_h1eff/get_h2eff both default to self.mo_coeff when None.
    h1e, ecore = mc.get_h1eff(mo_coeff=truncated_mo)
    h2e = ao2mo.restore(1, mc.get_h2eff(truncated_mo), ncas)
    # The actual orbitals h1e/h2e above were built from — truncated_mo when the
    # time budget cut CASSCF short, mc.mo_coeff (post-optimization) otherwise.
    # Returned so main() can persist it: until this, no code anywhere wrote the
    # rotated orbitals to disk, so a later SHCI cross-validation had no way to
    # reuse the SAME basis DMRG-SCF actually solved on — solange_shci.py always
    # built its own fixed, unrotated AVAS orbitals instead (see that script's
    # own module comment), making any DMRG-SCF-vs-SHCI energy delta partly an
    # artifact of comparing two different orbital bases, not a real solver
    # disagreement. Found live 2026-09-22 on TP53_R175_NATIVE.
    mo_used = truncated_mo if truncated else mc.mo_coeff
    if dmrg_scf:
        # Deliberately NOT re-verified against a full FCI diagonalization here —
        # that call (fci.direct_spin1.FCI(), below, for the plain-CASSCF path) is
        # exactly the combinatorial operation dmrg_scf exists to avoid; running
        # it anyway on a >16-orbital space would defeat the entire feature. The
        # consistency this loses is regained downstream: run_dmrg() in main()
        # sweeps the SAME (h1e, h2e, ecore) at the real --bond-dims and reports
        # its own S_max/convergence — that is this path's evidence, not e_fci.
        e_fci = None
    else:
        na = nelec // 2
        e_fci = fci.direct_spin1.FCI().kernel(h1e, h2e, ncas, (na, nelec - na), ecore=0.0)[0]
    return {"e_casscf": float(e_casscf), "ecore": float(ecore),
            "e_fci_active": None if e_fci is None else float(e_fci),
            "e_scf": float(mf.e_tot), "ncore": int(mc.ncore), "nelectron": int(mf.mol.nelectron),
            "h1e": h1e, "h2e": h2e, "ncas": int(ncas), "nelecas": int(nelec),
            "mo_coeff": mo_used,
            "orbital_optimization_truncated": bool(truncated),
            "orbital_optimization_converged": bool(casscf_converged),
            "orbital_optimization_method": (
                (f"CASCI, externally supplied fixed active-space orbitals from {load_orbitals} "
                 f"(deliberately built, e.g. avas_by_count.py/avas_mp2_select.py — optimization "
                 f"not applicable, not a truncated-run artifact) + " +
                 (f"DMRG solve (block2)" if dmrg_scf else "exact FCI solve") if external_orbitals else
                 f"CASCI, externally-loaded orbitals from {load_orbitals} (NOT re-optimized or "
                 f"re-verified as converged — inherited from a prior truncated run) + " +
                 (f"DMRG solve (block2)" if dmrg_scf else "exact FCI solve") if load_orbitals else
                 "CASCI, fixed AVAS orbitals (no optimization) + " +
                 (f"DMRG solve (block2)" if dmrg_scf else "exact FCI solve") if casci else
                 (f"DMRG-SCF (block2, maxM={dmrg_scf_maxm})" if dmrg_scf else "CASSCF (FCI solver)"))
                + (" — orbital optimisation STOPPED at its time budget, not converged"
                   if truncated else ""))}


# ── LEON seal (self-contained — mirrors backend/routes/leon.py's generic seal
# bit-for-bit, deliberately duplicated rather than imported: this script must run
# standalone on the Laguna cluster with no dependency on the SOLANGE backend
# package). LEON re-verifies this at ingestion; a mismatch is REJECTED, not stored.
def _seal_payload(record, exclude):
    return json.dumps({k: v for k, v in record.items() if k not in exclude},
                      sort_keys=True, default=str)


def _seal_hash(record, exclude):
    return hashlib.sha256(_seal_payload(record, exclude).encode()).hexdigest()


def _submit_dmrg(api, out):
    """POST the sealed DMRG classification to SOLANGE. Best-effort: prints the
    outcome but never raises — a submit failure must not discard the local JSON
    already written to --out."""
    import urllib.request
    try:
        url = api.rstrip("/") + "/api/simulate/hpc/dmrg/submit"
        body = json.dumps(out, default=str).encode()
        req = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            resp = json.loads(r.read().decode())
        print(f"SUBMITTED → {url}")
        print(f"  status={resp.get('status')}  seal_ok={resp.get('seal_ok')}  "
              f"db={resp.get('db_status')}  run_id={resp.get('run_id')}")
    except Exception as e:
        print(f"  SUBMIT FAILED (result is still safe locally in --out): {e}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="SOLANGE DMRG A/B/C classifier (Laguna largemem).")
    ap.add_argument("--compound", help="model compound (key in GEOM) — demo mode")
    ap.add_argument("--geometry", help="path to a QM-cluster .xyz (real functional site)")
    ap.add_argument("--avas", help="comma-separated AVAS target AOs, e.g. 'Zn 3d,S 3p'")
    ap.add_argument("--avas-threshold", type=float, default=0.2,
                    help="--geometry mode only: AVAS projection-score threshold (pyscf default "
                         "0.2, previously never exposed here — every prior run used it silently). "
                         "Widening this (0.1, 0.05) admits more virtual orbitals; see "
                         "integrals_from_geometry()'s docstring for why this needs to be an "
                         "explicit, recorded choice rather than an implicit library default.")
    ap.add_argument("--charge", type=int, default=0,
                    help="net molecular charge for --geometry mode (default 0/neutral). A QM "
                         "cluster cut from a protein rarely sums to a neutral closed shell by "
                         "accident — pyscf will error with 'Electron number N and spin S are "
                         "not consistent' if this is wrong; that error is the signal to set it.")
    ap.add_argument("--spin", type=int, default=0,
                    help="2S = Nalpha - Nbeta for --geometry mode (default 0, closed shell / "
                         "RHF). Non-zero switches to ROHF.")
    ap.add_argument("--no-density-fit", action="store_true",
                    help="--geometry mode only: disable density-fitted RHF/ROHF AND CASSCF "
                         "(both on by default — a QM cluster this size makes the conventional, "
                         "non-DF path slow for BOTH steps, not just RHF: CASSCF's orbital "
                         "optimization still touches every basis function each macro-iteration "
                         "even when the active space itself is tiny. Pass this only if an exact, "
                         "non-DF comparison is specifically needed — expect it to take hours.)")
    ap.add_argument("--df-auxbasis", default="def2-universal-jkfit",
                    help="--geometry mode only: auxiliary basis for density fitting (default "
                         "def2-universal-jkfit — a compact, general-purpose fitting basis that "
                         "works with any primary basis. Naming this explicitly matters most for "
                         "primary bases with no matching JKFIT companion in the standard "
                         "libraries, e.g. sto-3g: left to guess, pyscf falls back to an "
                         "oversized even-tempered auxiliary basis and DF can OOM instead of help)")
    ap.add_argument("--max-cycle-macro", type=int, default=20,
                    help="--geometry mode only: hard cap on CASSCF macro-iterations (default 20, "
                         "pyscf's own default is 50 with no time bound). This is a time ceiling, "
                         "not a convergence guarantee — the run prints whether it actually "
                         "converged or just hit this cap.")
    ap.add_argument("--max-memory", type=int, default=16000,
                    help="--geometry mode only: MB available to pyscf for SCF/CASSCF/FCI "
                         "(default 16000 — pyscf's own default of 4000 is sized for a "
                         "laptop, not a Laguna largemem node; raise this if you still see "
                         "a 'needs N MB, over max_memory limit' warning)")
    ap.add_argument("--casci", action="store_true",
                    help="--geometry mode only: freeze AVAS's orbitals as given and run ONE "
                         "diagonalization on them (CASCI) instead of iteratively optimizing "
                         "orbitals (CASSCF/DMRG-SCF). No macro-iteration loop at all, so the "
                         "warm-start solver pattern that crashed block2 at CAS(48,28) (see "
                         "integrals_from_geometry()'s casci branch) never triggers — the "
                         "fcisolver is called exactly once, cold. Trade-off: orbitals are AVAS's "
                         "raw output, not refined against the active-space energy. Combine with "
                         "--dmrg-scf to use the DMRG solver for that one solve instead of exact "
                         "FCI (still useful past ~16 orbitals, since CASCI's single solve is what "
                         "avoids the crash, not the choice of solver).")
    ap.add_argument("--load-orbitals", default=None,
                    help="--geometry mode only: path to a .npy mo_coeff matrix (as saved by a "
                         "prior run of THIS script to <its --scratch>/mo_coeff_final.npy) to use "
                         "in place of AVAS's own output, skipping orbital optimization entirely "
                         "(one fixed diagonalization, like --casci, but on these orbitals instead "
                         "of AVAS's raw ones). Reported as NOT converged regardless of whether the "
                         "source run's own message said so, since these orbitals were, by "
                         "construction, still mid-optimization when that run stopped. Exists to "
                         "let a resubmission after an infrastructure failure (Slurm --time kill, "
                         "OOM) reuse orbitals a truncated run already spent real compute reaching, "
                         "instead of redoing the whole (multi-hour) orbital-optimization phase — "
                         "found needed live 2026-09-22 after an OOM killed the DMRG ladder phase of "
                         "a run whose orbital phase alone had already taken ~14 hours.")
    ap.add_argument("--active-start-col", type=int, default=None,
                    help="--load-orbitals mode only, REQUIRED when the file came from "
                         "avas_mp2_select.py/avas_by_count.py (its 'active_start_col' field): sets "
                         "mc.ncore explicitly instead of pyscf's default ncore=(mol.nelectron-"
                         "nelecas)//2, which assumes the active block sits immediately below the "
                         "Fermi level. avas_mp2_select.py's build_mo() does NOT follow that "
                         "convention -- it places the active columns wherever AVAS's pool happened "
                         "to start and drops every non-selected pool candidate, so the default "
                         "formula silently reads the WRONG 34 columns. Omit only when --load-orbitals "
                         "points to a prior solange_dmrg.py run's own mo_coeff_final.npy, which IS in "
                         "standard order.")
    ap.add_argument("--dmrg-scf", action="store_true",
                    help="use DMRG (block2) as CASSCF's own orbital-optimization solver, "
                         "instead of FCI — this is what actually lets --ncas grow past the "
                         "~16-orbital FCI wall (see integrals_from_geometry()'s docstring); "
                         "uses the pyblock2 adapter in dmrgscf_block2.py, which cross-checks "
                         "itself against exact FCI before optimizing. Applies to both --geometry and "
                         "--compound/demo mode. The downstream DMRG sweep (run_dmrg(), at "
                         "--bond-dims) is unaffected either way — this flag only changes how "
                         "the active space's orbitals themselves get optimized.")
    ap.add_argument("--dmrg-scf-maxm", type=int, default=250,
                    help="bond dimension used DURING orbital optimization when --dmrg-scf is "
                         "set (default 250) — deliberately separate from --bond-dims, which is "
                         "the (usually larger) sweep run afterward on the resulting fixed "
                         "integrals for the reported energy/S_max. Raise this only if orbital "
                         "optimization itself fails to converge.")
    ap.add_argument("--dmrg-scf-scratch", default="./tmp_dmrgscf_orb",
                    help="block2 scratch dir for the --dmrg-scf orbital-optimization solver — "
                         "kept separate from --scratch (the downstream energy-sweep scratch) "
                         "so the two DMRG runs never collide over the same MPS files.")
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--ncas", type=int)
    ap.add_argument("--nelecas", type=int)
    ap.add_argument("--key", required=True, help="target key, e.g. ARID2_LOF")
    ap.add_argument("--side", default="native", choices=["native", "mutant"],
                    help="allele — used to auto-resolve the model compound from --key "
                         "when --compound is omitted (same mapping the HPC agent uses)")
    ap.add_argument("--no-early-stop", action="store_true",
                    help="run every requested bond dim even after the energy converges "
                         "(default: stop once ΔE between consecutive M < chemical accuracy, "
                         "since larger M cannot change the verdict)")
    ap.add_argument("--bond-dims", default="250,500,1000,2000",
                    help="comma-separated increasing bond dimensions")
    ap.add_argument("--n-sweeps", type=int, default=10,
                    help="DMRG sweeps per bond dimension (was hardcoded to 10; the project's "
                         "own calibration notes this default is itself uncalibrated, and that "
                         "n_sweeps=2 is a validated FAST TRIAGE screen, not a measurement -- see "
                         "the skill's classifier-parameters.md. Raise this (e.g. 20-30) for a "
                         "~48-orbital active space starting from a random MPS, where 10 sweeps "
                         "may not reach convergence and an unconverged run can look exactly like "
                         "the R175H cold-start artifact (large ΔE that is initialization noise, "
                         "not chemistry).")
    ap.add_argument("--out", default="./out")
    ap.add_argument("--verbose", type=int, default=0)
    ap.add_argument("--stack-mem-gb", type=float, default=DEFAULT_STACK_MEM_GB,
                    help=f"block2's own pre-allocated memory pool in GB (default "
                         f"{DEFAULT_STACK_MEM_GB}). This is internal to block2 and unrelated "
                         "to --max-memory, which bounds pyscf: exceeding it aborts the "
                         "process with 'exceeding allowed memory' and no dmesg trace. The "
                         "default was chosen to sit under a per-user cgroup ceiling that no "
                         "longer applies, so raise it for large sites - the reachable bond "
                         "dimension scales with it.")
    ap.add_argument("--threads", type=int, default=4,
                    help="CPU threads for block2 (was hard-coded to 4; raise this to use "
                         "more of a largemem node's cores and finish faster)")
    ap.add_argument("--scratch", default="./tmp_dmrg",
                    help="block2 scratch dir. Reuse the SAME path across separate HPC-ticket "
                         "submissions to resume an interrupted MPS instead of restarting from "
                         "bond_dims[0] each time.")
    ap.add_argument("--orbital-max-minutes", type=float, default=None,
                    help="wall-clock budget for the ORBITAL phase specifically, checked once "
                         "per CASSCF macro-iteration. Defaults to 60%% of --max-minutes. "
                         "--max-minutes alone never covered this phase: it bounds the "
                         "bond-dimension ladder, which runs afterwards, so an orbital "
                         "optimisation that overran simply got the whole job killed before "
                         "the ladder - and before the mechanism meant to save a partial "
                         "result - ever started.")
    ap.add_argument("--max-minutes", type=float, default=None,
                    help="stop requesting larger bond dimensions once this wall-clock budget "
                         "is exceeded, and return/save whatever completed so far — for "
                         "fixed-walltime HPC allocations (e.g. a 2-hour ticket) where getting "
                         "killed mid-sweep would lose everything.")
    ap.add_argument("--skip-precondition", action="store_true",
                    help="skip check_active_space.py's blocking precondition (R1 zero virtuals, "
                         "R2 virtual-fraction floor, R3 entropy ceiling below S_HARD, R4 exactly-"
                         "diagonalisable). Added 2026-09-30 after TP53_C275F cluster76 CAS(76,38) "
                         "ran to a S_max=0.0 conclusion that was mathematically guaranteed by its "
                         "occ/vir dimensions before the run started. Use this ONLY for a "
                         "deliberate boundary/calibration case where the verdict being moot is "
                         "the point, not to push past a failing check on a run meant to classify.")
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the sealed classification to SOLANGE so it's LEON-notarized and "
                         "stored immediately — safe even if this Laguna session later becomes "
                         "unreachable. Bare flag uses the default SOLANGE API URL; pass a value "
                         "to override.")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    bond_dims = [int(x) for x in args.bond_dims.split(",")]
    script_start = time.time()   # the ONE clock --max-minutes is a budget against — see below

    print("=" * 68)
    mo_path = None   # set below only for a --geometry run whose orbitals were saved
    if args.geometry:
        if not args.avas:
            ap.error("--geometry requires --avas (target AOs for active-space selection)")
        print(f"SOLANGE DMRG classifier · {args.key} · geometry={args.geometry} "
              f"· AVAS[{args.avas}]/{args.basis} · charge={args.charge} spin={args.spin}")
        # Default: 60% of --max-minutes to the orbital phase, the rest to the
        # bond-dimension ladder — a split, not a guess free of consequence, chosen
        # because the ladder's own early-stop makes it the cheaper phase to cut
        # short if the split is wrong, while a truncated orbital phase produces no
        # coherent handoff to the ladder at all.
        orbital_deadline = None
        if args.orbital_max_minutes is not None:
            orbital_deadline = time.time() + args.orbital_max_minutes * 60
        elif args.max_minutes is not None:
            orbital_deadline = time.time() + 0.6 * args.max_minutes * 60

        cas = integrals_from_geometry(args.geometry, args.basis, args.avas,
                                       charge=args.charge, spin=args.spin, verbose=args.verbose,
                                       density_fit=not args.no_density_fit, max_memory=args.max_memory,
                                       df_auxbasis=args.df_auxbasis, max_cycle_macro=args.max_cycle_macro,
                                       dmrg_scf=args.dmrg_scf, dmrg_scf_maxm=args.dmrg_scf_maxm,
                                       dmrg_scf_scratch=args.dmrg_scf_scratch, n_threads=args.threads,
                                       orbital_deadline=orbital_deadline, stack_mem_gb=args.stack_mem_gb,
                                       casci=args.casci, avas_threshold=args.avas_threshold,
                                       load_orbitals=args.load_orbitals,
                                       ncas_override=args.ncas if args.load_orbitals else None,
                                       nelecas_override=args.nelecas if args.load_orbitals else None,
                                       ncore_override=args.active_start_col,
                                       skip_precondition=args.skip_precondition)
        args.ncas, args.nelecas = cas["ncas"], cas["nelecas"]
        print(f"AVAS selected active space: CAS({args.nelecas},{args.ncas})")
        if cas.get("orbital_optimization_truncated"):
            print("  *** orbital optimisation hit its time budget — orbitals are usable but "
                  "NOT converged; the classification below is PROVISIONAL on that basis too ***")
        # Persist the ACTUAL rotated orbitals h1e/h2e were built from — in
        # --scratch (the stable, --key-named directory), not dmrg_scf_scratch
        # (block2's own internal MPS state, which does not hold this matrix at
        # all — confirmed live 2026-09-22 while investigating why a DMRG-SCF
        # run's orbitals could not be recovered after the fact for an SHCI
        # cross-validation). solange_shci.py's own --orbitals flag loads this
        # file to build h1e/h2e from the SAME basis DMRG-SCF solved on, instead
        # of its default fixed/unrotated AVAS orbitals.
        if cas.get("mo_coeff") is not None:
            Path(args.scratch).mkdir(parents=True, exist_ok=True)
            mo_path = Path(args.scratch) / "mo_coeff_final.npy"
            np.save(mo_path, cas["mo_coeff"])
            print(f"  saved rotated orbitals -> {mo_path} (for a matching-basis SHCI cross-validation)")
    else:
        # Auto-resolve the model compound from key/side (same mapping the HPC agent
        # uses) so a run needs only --key/--side/--ncas/--nelecas — the caller does
        # not have to know which GEOM compound models this gene.
        if not args.compound:
            try:
                from solange_hpc import _resolve_compound
                args.compound = _resolve_compound(args.key, args.side)
                print(f"--compound not given → resolved {args.key}/{args.side} → {args.compound}")
            except Exception as e:
                ap.error(f"could not resolve a model compound for {args.key}/{args.side} ({e}). "
                         f"Pass --compound <GEOM key> explicitly, or --geometry <xyz>+--avas.")
        if not (args.ncas and args.nelecas):
            ap.error("provide --ncas and --nelecas (or use --geometry+--avas for AVAS selection)")
        from solange_hpc import run_casscf
        print(f"SOLANGE DMRG classifier · {args.key} · {args.compound}/{args.basis} "
              f"· CAS({args.nelecas},{args.ncas})")
        if not args.skip_precondition:
            sys.path.insert(0, str(_HERE.parent))
            import check_active_space as cas_check
            _p = argparse.Namespace(s_hard=S_HARD, min_vir_frac=15.0, fci_floor=1e9, fci_ceiling=1e40)
            _fails = cas_check.report(args.nelecas, args.ncas, _p)
            if _fails:
                sys.exit(f"\n*** REFUSING to run CASSCF/DMRG-SCF on CAS({args.nelecas},{args.ncas}): "
                         f"{'; '.join(_fails)}. Pass --skip-precondition only if this is a "
                         f"deliberate boundary/calibration case, not a classification run. ***")
        # Same split as the --geometry branch above, and for the same reason: this
        # branch had NO orbital-phase time guard at all until a live --compound run
        # (ARID2_LOF, CAS(16,16)) ground past macro-iteration 170 with the energy
        # frozen — conv_tol was already satisfied, but conv_tol_grad, computed from
        # DMRG's own RDM noise floor, may never fall below its threshold, so the
        # loop had no way to stop short of the 200-macro-iteration cap.
        compound_orbital_deadline = None
        if args.orbital_max_minutes is not None:
            compound_orbital_deadline = time.time() + args.orbital_max_minutes * 60
        elif args.max_minutes is not None:
            compound_orbital_deadline = time.time() + 0.6 * args.max_minutes * 60
        cas = run_casscf(args.compound, args.basis, args.ncas, args.nelecas, args.verbose,
                         dmrg_scf=args.dmrg_scf, dmrg_scf_maxm=args.dmrg_scf_maxm,
                         dmrg_scf_scratch=args.dmrg_scf_scratch, n_threads=args.threads,
                         orbital_deadline=compound_orbital_deadline)
        if cas.get("orbital_optimization_truncated"):
            print("  *** orbital optimisation hit its time budget — orbitals are usable but "
                  "NOT converged; the classification below is PROVISIONAL on that basis too ***")
    # Two independent keys across the two code paths this branches into
    # (integrals_from_geometry's own "orbital_optimization_converged" vs.
    # solange_hpc.run_casscf's "converged") - normalized here once. A run
    # that hit max_cycle_macro without raising OrbitalTimeBudgetExceeded used
    # to fall through this entirely: orbital_optimization_truncated stayed
    # False (that flag only ever fires on the time-budget path), so neither
    # branch's own warning above ever printed, and the "FINAL result" message
    # below was reached with orbitals that were never actually checked.
    # exactly what happened live 2026-09-12 on TP53_C275_NATIVE, whose own
    # log said "CASSCF did NOT converge" one line above a run this code
    # would otherwise have called final.
    orbital_converged = cas.get("orbital_optimization_converged",
                                 cas.get("converged", True))
    if not orbital_converged and not cas.get("orbital_optimization_truncated"):
        print("  *** orbital optimisation hit max_cycle_macro without converging (not a time-"
              "budget stop) — orbitals are usable but NOT converged; the classification below "
              "is PROVISIONAL on that basis too ***")
    print(f"{cas['orbital_optimization_method']} E = {cas['e_casscf']:.8f} Ha")

    t0 = time.time()
    # run_dmrg()'s own max_minutes check starts counting from ITS OWN t_start (set
    # when it's called), not from this script's start — so passing args.max_minutes
    # through unchanged gave the ladder phase a FRESH full budget on top of whatever
    # the orbital-optimization phase (0.6*args.max_minutes, above) already spent,
    # instead of the REMAINING share of one shared budget the "60%/rest" design
    # comment above describes. A real bug, not a documentation gap: found live
    # 2026-09-22 when TP53_R175H's job ran orbitals for ~368m then was still mid-M=1000
    # (nowhere near a graceful stop) when Slurm's OWN --time=720m hard-killed it at
    # 720m with NOTHING saved — the internal graceful stop, at what was actually a
    # ~968m (368+600) effective budget, never had a chance to fire first. Computed
    # here from the REAL elapsed time (not the nominal 60/40 split), so it accounts
    # correctly regardless of whether the orbital phase finished early, late, or hit
    # its own time budget. max(0, ...) so an orbital phase that already overran the
    # WHOLE budget still lets the ladder attempt bond_dims[0] once (run_dmrg's own
    # "and energies" guard never stops before at least one M completes) rather than
    # passing a negative number.
    ladder_max_minutes = (
        max(0.0, args.max_minutes - (time.time() - script_start) / 60.0)
        if args.max_minutes is not None else None)
    energies, s_max, stop_reason, discarded_weights, sweep_history = run_dmrg(
                               cas["h1e"], cas["h2e"], cas["ecore"],
                               args.ncas, args.nelecas, bond_dims,
                               scratch=args.scratch, n_threads=args.threads,
                               max_minutes=ladder_max_minutes,
                               early_stop=not args.no_early_stop,
                               stack_mem_gb=args.stack_mem_gb, n_sweeps=args.n_sweeps)
    # (per-M timing is already printed live inside run_dmrg, as each M finishes —
    # so a `tail -f` on a background run shows real progress, not a single dump at exit.)
    print(f"max bipartite entanglement S_max = {s_max}")
    # Discarded-weight extrapolation (Claude Science, 2026-09-30): a sturdier
    # convergence read than raw ΔE, which conflates truncation error with an
    # unconverged sweep. Reported ALONGSIDE classify()'s own ΔE-based verdict,
    # not in place of it. The FIELD is now settled from pyblock2's own source
    # (drv.get_dmrg_results(), see above) -- discarded_weights is genuinely
    # per-sweep, sweep_discarded_weights was a red herring never touched by the
    # Python driver. What's NOT yet re-validated on a real (non-toy) active
    # space is whether THIS classification's bond-dimension ladder sits in the
    # asymptotic linear-in-w regime the fit assumes (the H8/CAS(8,8) validator
    # that first tried this field failed only because it had no such regime,
    # not because the field was wrong) -- window_stability() below is exactly
    # that check, run on this run's own ladder rather than assumed.
    extrap = None
    if all(w is not None for w in discarded_weights):
        try:
            import dmrg_extrapolate
            Ms = [m for m, _ in energies]
            Es = [e for _, e in energies]
            extrap = dmrg_extrapolate.extrapolate(Ms, Es, discarded_weights)
            if extrap.get("ok"):
                print(f"  [discarded-weight extrapolation] {extrap['verdict']} "
                      f"(residual {extrap['residual_mHa']:.3f} mHa vs chemical accuracy "
                      f"{extrap['chem_acc_Ha']*1000:.1f} mHa, R²={extrap['r2']:.4f}, "
                      f"w spans {extrap['w_decades']:.2f} decades)")
                for w in extrap.get("warnings", []):
                    print(f"  [discarded-weight extrapolation] WARNING: {w}")
            else:
                print(f"  [discarded-weight extrapolation] not run: {extrap['reason']}")
            if len(Ms) >= 4:
                stab = dmrg_extrapolate.window_stability(Ms, Es, discarded_weights)
                if stab.get("ok"):
                    print(f"  [window stability] {stab['verdict']} "
                          f"(intercept spread {stab['intercept_spread_mHa']:.3f} mHa "
                          f"across {len(stab['windows'])} fit windows)")
                    if extrap:
                        extrap["window_stability"] = stab
                else:
                    print(f"  [window stability] not run: {stab['reason']}")
        except Exception as e:
            print(f"  [discarded-weight extrapolation] skipped ({type(e).__name__}: {e})")
    else:
        print("  [discarded-weight extrapolation] skipped — discarded weight unavailable "
              "for one or more bond dimensions on this pyblock2 build")
    cap_pct = entanglement_capacity_pct(s_max, args.ncas, args.nelecas, args.spin)
    if cap_pct is not None:
        print(f"  ({cap_pct:.1f}% of this active space's own exact combinatorial "
              f"entanglement capacity (filling-aware, not orbital-count-only) — "
              f"size-normalised diagnostic, does not change the classification "
              f"below; see S_HARD's own comment)")
    # PROVISIONAL only when the sweep was cut short by the TIME budget — a convergence
    # early-stop is the opposite (the answer is final, we just skipped redundant high M).
    time_budget_hit = (stop_reason == "time_budget")
    if time_budget_hit:
        print(f"NOTE: stopped early at {len(energies)}/{len(bond_dims)} bond dimensions "
              f"(--max-minutes {args.max_minutes}). The class below is PROVISIONAL — "
              f"re-run with --scratch {args.scratch} to resume toward the full bond-dim list.")
    elif stop_reason == "converged" and orbital_converged:
        print(f"NOTE: converged early at M={energies[-1][0]} "
              f"({len(energies)}/{len(bond_dims)} bond dims run) — larger M cannot change "
              f"the verdict. This is a FINAL result, not provisional.")
    elif stop_reason == "converged" and not orbital_converged:
        print(f"NOTE: DMRG itself converged at M={energies[-1][0]}, but the ORBITALS underneath "
              f"it did not (see the orbital-optimisation warning above) — this is PROVISIONAL, "
              f"not final: a different orbital solution could give a different S_max.")

    orbital_provisional = not orbital_converged and not cas.get("orbital_optimization_truncated")
    # Mandatory preflight (Claude Science, 2026-10-03): the first run of this exact
    # ladder emitted CLASS A with "[discarded-weight extrapolation] skipped
    # (AttributeError: module 'dmrg_extrapolate' has no attribute 'extrapolate')"
    # printed one line above its own verdict -- a guard module had been silently
    # truncated to 0 bytes by a bad sync (ast.parse("") is valid syntax, so even a
    # syntax check passes). The exception was caught, "skipped" was printed
    # truthfully, and the verdict was emitted anyway with no convergence criterion
    # behind it. require_guards() makes that failure mode structural: if any guard
    # module cannot be imported, is zero bytes, or is missing a required symbol,
    # no verdict may be emitted at all. classify_preflight() separately checks the
    # LADDER itself against E_SCF -- the embedding gate above only ever validated
    # the reference solve, not the points classify() is about to read -- which is
    # exactly how seven ladder points sitting above E_SCF reached a verdict
    # unchecked on an earlier run of this same target.
    import guard_preflight
    guard_preflight.require_guards()
    if "e_scf" in cas and "ncore" in cas:
        _sweep_bd, _sweep_e = [], []
        for _m, _h in sweep_history.items():
            for _e in _h["energies"]:
                _sweep_bd.append(int(_m)); _sweep_e.append(_e)
        _pf = guard_preflight.classify_preflight(
            cas["e_casscf"], cas["e_scf"], label=args.key,
            ncore=cas["ncore"], nelectron=cas["nelectron"], nelecas=args.nelecas,
            bond_dims=[m for m, _ in energies], energies=[e for _, e in energies],
            discarded_weights=discarded_weights, s_max=s_max,
            sweep_bond_dims=_sweep_bd or None, sweep_energies=_sweep_e or None)
        for _b in _pf["blocking"]:
            print(f"  [preflight] BLOCKING: {_b}", flush=True)
        if not _pf["may_emit_verdict"]:
            sys.exit(f"\n*** REFUSING: classify_preflight found {len(_pf['blocking'])} blocking "
                     f"finding(s) on {args.key} -- {_pf['verdict']}. No class may be emitted "
                     f"from this run. ***")
    else:
        print("  [preflight] e_scf/ncore not recorded by this integrals path (--compound demo "
              "mode) — classify_preflight's ladder-vs-E_SCF check skipped; the embedding gate "
              "above still ran.", flush=True)
    cls, rationale = classify(args.nelecas, energies, s_max)
    print("-" * 68)
    # Class A is ALWAYS provisional (2026-10-09): classical failure "at any budget"
    # cannot be shown for a finite instance, and no fault-tolerant QPU exists to
    # show quantum success, so the strongest honest claim is "unresolved within
    # the declared budget". Class B, by contrast, can be certified.
    provisional_tag = (' (PROVISIONAL — time budget hit)' if time_budget_hit
                        else ' (PROVISIONAL — orbital optimization did not converge)' if orbital_provisional
                        else ' (PROVISIONAL — unresolved within declared budget, not a proof)' if cls == 'A'
                        else '')
    print(f"CLASS {cls}{provisional_tag}")
    print(f"  {rationale}")
    elapsed_s = round(time.time() - t0, 1)
    print(f"elapsed {elapsed_s}s")

    out = {
        "id": str(uuid.uuid4()),
        "key": args.key, "side": args.side, "compound": args.compound, "basis": args.basis,
        "ncas": args.ncas, "nelecas": args.nelecas,
        "e_casscf": cas["e_casscf"],
        "dmrg_energies": energies, "s_max": s_max,
        "discarded_weights": discarded_weights,
        "discarded_weight_extrapolation": extrap,
        "s_max_capacity_pct": entanglement_capacity_pct(s_max, args.ncas, args.nelecas, args.spin),
        "bqp_class": cls, "class_rationale": rationale,
        "time_budget_hit": time_budget_hit, "orbital_optimization_converged": orbital_converged,
        "bond_dims_requested": bond_dims,
        "elapsed_s": elapsed_s,
        "method": "DMRG (block2, classical) convergence + entanglement diagnostic",
        "orbital_optimization_method": cas["orbital_optimization_method"],
        # Was hardcoded "HPC/Laguna" regardless of where this actually ran --
        # harmless while every agent was this one Laguna account, but this
        # module has no functional Laguna dependency (confirmed 2026-10-03
        # while answering whether a chemist could point this at a different
        # HPC cluster via the hpc_dispatch pull-queue architecture): a
        # different agent/cluster would have silently misreported itself as
        # Laguna. Derived from the real hostname instead, same source
        # detect_hardware() already uses for the hardware column below.
        "provenance_source": f"HPC/{socket.gethostname()} (DMRG classifier)",
        "hardware": detect_hardware(args.threads),
    }
    # Recorded ONLY for a real --geometry run (not --compound demo mode): these
    # are exactly the reproducibility inputs an SHCI cross-validation needs to
    # rebuild the IDENTICAL active-space Hamiltonian later (solange_shci.py) --
    # without them, "Queue SHCI Cross-Validation" against this record would have
    # nothing to copy the geometry/AVAS/charge/spin from and would need them
    # re-typed by hand, exactly the friction this field exists to remove.
    if args.geometry:
        # Store the RAW .xyz content, not the local path — this record may be
        # read back much later (a different machine, a different run's scratch
        # dir) to re-dispatch the SAME active space to Rung 4 (QPU), and the
        # QPU agent writes this field verbatim as a new file's content
        # (solange_qpu.py's run_agent: geom_path.write_text(job["geometry"])).
        # A bare path string here would silently produce a garbage .xyz file.
        out["geometry"] = Path(args.geometry).read_text()
        out["avas"] = args.avas
        out["charge"] = args.charge
        out["spin"] = args.spin
        out["avas_threshold"] = args.avas_threshold
        # Unlike geometry/avas/charge/spin above, this is stored as a LOCAL PATH,
        # not portable content — an mo_coeff matrix is a binary numpy array, not
        # something that fits a JSON field the way the raw .xyz text does. This
        # is safe only because the whole SHCI cross-validation flow already runs
        # on this SAME Laguna filesystem (solange_hpc.py's agent dispatches
        # solange_shci.py as a local subprocess, never remotely) — the path is
        # meaningless off this cluster. None when orbitals weren't saved (no
        # --geometry run predating this field, or mo_coeff unexpectedly absent).
        out["orbitals_path"] = str(mo_path) if mo_path is not None else None
    # Seal at source (LEON re-verifies at ingestion — a mismatch is rejected, not
    # trusted). dmrg_seal_payload is stored verbatim so re-verification later is
    # exact-string, not float-reconstruction (the same robustness fix the P8 seal
    # got after floats/timestamps proved to reformat across a DB round-trip).
    out["dmrg_seal_payload"] = _seal_payload(out, exclude={"dmrg_hash", "dmrg_seal_payload"})
    out["dmrg_hash"] = hashlib.sha256(out["dmrg_seal_payload"].encode()).hexdigest()

    p = Path(args.out) / f"dmrg_class_{args.key}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"WROTE {p}")
    if args.submit:
        _submit_dmrg(args.submit, out)
    print("=" * 68)
    print("NOTE: DMRG is classical (CPU/largemem or GPU-accelerated), NOT quantum.")
    print("Full-site classification needs the target's active space built from PDB first.")


if __name__ == "__main__":
    main()
