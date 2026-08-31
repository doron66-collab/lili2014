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
# For scale, cost runs as M^3, so 2000 is 512x the M=250 the DMRG-SCF solver
# uses by default, and roughly what a single node reaches in hours; the largest
# published chemistry DMRG runs go to M=6000 on thousands of cores for weeks
# (FeMoco, CAS(113,76)). Calibration on an N2 stretch series at CAS(10,16)
# found the requirement in the 160-256 range, an order of magnitude under this
# ceiling - so it is generous for small systems and untested against large ones.
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
# 44-electron site.
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
             stack_mem_gb=DEFAULT_STACK_MEM_GB):
    """Run DMRG at increasing bond dimensions. Returns per-M energies + S_max.

    HPC-ticket-aware: prints live per-M timing (so `tail -f` shows real progress,
    letting you judge whether to Ctrl+C before a fixed-walltime allocation ends),
    and honors max_minutes — once the cumulative wall-clock budget is exceeded, it
    stops requesting new (larger) bond dimensions and returns whatever it has
    rather than being killed mid-sweep with nothing recorded. Reusing the same
    `scratch` directory across separate job submissions lets block2 resume the MPS
    from where a prior run left off instead of restarting from bond_dims[0].
    """
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=n_threads,
                     stack_mem=int(stack_mem_gb * (1 << 30)))
    drv.initialize_system(n_sites=ncas, n_elec=nelecas, spin=0)
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
    try:
        ket = drv.load_mps(tag="KET", nroots=1)
        print(f"  [resume] loaded existing MPS from {scratch} — continuing, not restarting "
              f"the bond-dimension ladder from scratch.", flush=True)
    except Exception as e:
        ket = drv.get_random_mps(tag="KET", bond_dim=min(bond_dims[0], 250), nroots=1)
        print(f"  [resume] no usable saved MPS in {scratch} ({type(e).__name__}) — "
              f"starting from a fresh random MPS.", flush=True)
    energies = []
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
        e = drv.dmrg(mpo, ket, n_sweeps=10, bond_dims=[M],
                     noises=[1e-5, 1e-6, 0], thrds=[1e-9] * 3, iprint=0)
        dt = time.time() - t0
        energies.append((M, float(e)))
        # stdout (not stderr) so the agent's live-progress streamer reliably captures
        # each "DMRG M=" line and surfaces it as a bond-dimension stage note.
        print(f"  DMRG M={M:5d}  E={float(e):.8f} Ha  [{dt:.1f}s this M, "
              f"{(time.time()-t_start)/60:.1f}m total]", flush=True)
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
    return energies, s_max, stop_reason


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
def classify(active_electrons, energies, s_max):
    """Map DMRG behaviour to an A/B/C class with an explicit rationale.

    Two signals:
      • ΔE across the last two bond dims — direct evidence DMRG has (not) converged
        at a practical M. Only bites at large active spaces; small ones are exact.
      • S_max (max bipartite entanglement) — predicts the bond dimension the FULL
        site needs (M ~ e^S). This is the leading indicator, since a small model
        space is always DMRG-exact yet still reveals the correlation strength.
    """
    dE_final = abs(energies[-1][1] - energies[-2][1]) * 1000 if len(energies) >= 2 else None
    converged = (dE_final is not None) and (dE_final < CHEM_ACC_MHA) \
                and (energies[-1][0] <= PRACTICAL_M)
    strong = (s_max is not None) and (s_max > S_HARD)
    smx = "n/a" if s_max is None else f"{s_max:.2f}"

    if active_electrons <= EXACT_WALL_E:
        return "C", (f"{active_electrons}e ≤ {EXACT_WALL_E}e exact-classical wall — "
                     f"CCSD(T)/FCI reaches chemical accuracy; classical sufficient. "
                     f"(S_max={smx})")
    if converged and not strong:
        return "B", (f"DMRG reaches chemical accuracy at practical M={energies[-1][0]} "
                     f"(ΔE={dE_final:.2f} mHa) and entanglement is low (S_max={smx} < "
                     f"{S_HARD}) — classical (DMRG) delivers; quantum-advantaged, not necessary.")
    reasons = []
    if strong:
        reasons.append(f"S_max={smx} > {S_HARD} (strong correlation; DMRG bond dim ~e^S "
                       f"blows up at the full {active_electrons}e site)")
    if not converged:
        reasons.append(f"DMRG not at chemical accuracy by practical M={energies[-1][0]} "
                       f"(ΔE={'n/a' if dE_final is None else round(dE_final,2)} mHa)")
    return "A", "quantum-necessary — " + "; ".join(reasons) + "."


def integrals_from_geometry(xyz_path, basis, avas_aos, charge=0, spin=0, verbose=0,
                             density_fit=True, max_memory=16000, df_auxbasis="def2-universal-jkfit",
                             max_cycle_macro=20, dmrg_scf=False, dmrg_scf_maxm=500,
                             dmrg_scf_scratch="./tmp_dmrgscf_orb", n_threads=4,
                             orbital_deadline=None, stack_mem_gb=None, casci=False):
    """Chemist-in-the-loop entry: given a QM-cluster geometry (xyz) and the target
    atomic orbitals, AVAS selects the active space automatically. Returns a dict
    shaped like run_casscf's output. The CLUSTER itself (which residues/atoms/metal,
    H-capping) is the chemist's input — that step is NOT auto-generated here.

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
    That is what actually lets --ncas grow past 16 toward a real ~44-155e site;
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
    ncas, nelec, mo = avas.avas(mf, [s.strip() for s in avas_aos.split(",")])
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
    from pyscf import mcscf
    if casci:
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
        print(f"  --casci: skipping orbital optimization, one fixed diagonalization "
              f"on AVAS's orbitals as given.", flush=True)
        mc = mcscf.CASCI(mf, ncas, nelec)
    else:
        mc = mcscf.CASSCF(mf, ncas, nelec)
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
    if spin == 0 and not dmrg_scf and not casci:
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
    try:
        e_casscf = mc.kernel(mo)[0]
        conv_note = ("converged" if mc.converged else
                     f"did NOT converge — hit max_cycle_macro={max_cycle_macro} first")
    except OrbitalTimeBudgetExceeded as exc:
        truncated = True
        # NOT mc.e_tot — see the class docstring: never populated on an
        # interrupted run, so float(mc.e_tot) raises TypeError on None instead
        # of the guard doing its job.
        e_casscf = float(exc.e_tot)
        truncated_mo = exc.mo   # likewise: mc.mo_coeff is stale (pre-optimization)
        conv_note = (f"STOPPED at the orbital time budget after macro-iteration {exc.imacro} "
                     f"— orbitals are usable but NOT converged")
    print(f"  {'CASCI (fixed AVAS orbitals)' if casci else 'CASSCF'} {conv_note} · E={e_casscf:.8f}", flush=True)
    # Explicit mo_coeff=truncated_mo when cut short — see OrbitalTimeBudgetExceeded's
    # docstring; None (the normal-completion case) is identical to omitting the
    # argument, since get_h1eff/get_h2eff both default to self.mo_coeff when None.
    h1e, ecore = mc.get_h1eff(mo_coeff=truncated_mo)
    h2e = ao2mo.restore(1, mc.get_h2eff(truncated_mo), ncas)
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
            "h1e": h1e, "h2e": h2e, "ncas": int(ncas), "nelecas": int(nelec),
            "orbital_optimization_truncated": bool(truncated),
            "orbital_optimization_method": (
                ("CASCI, fixed AVAS orbitals (no optimization) + " +
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
    ap.add_argument("--submit", nargs="?", const="https://qcaihpc-simulation-api.onrender.com",
                    help="POST the sealed classification to SOLANGE so it's LEON-notarized and "
                         "stored immediately — safe even if this Laguna session later becomes "
                         "unreachable. Bare flag uses the default SOLANGE API URL; pass a value "
                         "to override.")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    bond_dims = [int(x) for x in args.bond_dims.split(",")]

    print("=" * 68)
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
                                       casci=args.casci)
        args.ncas, args.nelecas = cas["ncas"], cas["nelecas"]
        print(f"AVAS selected active space: CAS({args.nelecas},{args.ncas})")
        if cas.get("orbital_optimization_truncated"):
            print("  *** orbital optimisation hit its time budget — orbitals are usable but "
                  "NOT converged; the classification below is PROVISIONAL on that basis too ***")
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
    print(f"{cas['orbital_optimization_method']} E = {cas['e_casscf']:.8f} Ha")

    t0 = time.time()
    energies, s_max, stop_reason = run_dmrg(cas["h1e"], cas["h2e"], cas["ecore"],
                               args.ncas, args.nelecas, bond_dims,
                               scratch=args.scratch, n_threads=args.threads,
                               max_minutes=args.max_minutes,
                               early_stop=not args.no_early_stop,
                               stack_mem_gb=args.stack_mem_gb)
    # (per-M timing is already printed live inside run_dmrg, as each M finishes —
    # so a `tail -f` on a background run shows real progress, not a single dump at exit.)
    print(f"max bipartite entanglement S_max = {s_max}")
    # PROVISIONAL only when the sweep was cut short by the TIME budget — a convergence
    # early-stop is the opposite (the answer is final, we just skipped redundant high M).
    time_budget_hit = (stop_reason == "time_budget")
    if time_budget_hit:
        print(f"NOTE: stopped early at {len(energies)}/{len(bond_dims)} bond dimensions "
              f"(--max-minutes {args.max_minutes}). The class below is PROVISIONAL — "
              f"re-run with --scratch {args.scratch} to resume toward the full bond-dim list.")
    elif stop_reason == "converged":
        print(f"NOTE: converged early at M={energies[-1][0]} "
              f"({len(energies)}/{len(bond_dims)} bond dims run) — larger M cannot change "
              f"the verdict. This is a FINAL result, not provisional.")

    cls, rationale = classify(args.nelecas, energies, s_max)
    print("-" * 68)
    print(f"CLASS {cls}{' (PROVISIONAL — time budget hit)' if time_budget_hit else ''}")
    print(f"  {rationale}")
    elapsed_s = round(time.time() - t0, 1)
    print(f"elapsed {elapsed_s}s")

    out = {
        "id": str(uuid.uuid4()),
        "key": args.key, "side": args.side, "compound": args.compound, "basis": args.basis,
        "ncas": args.ncas, "nelecas": args.nelecas,
        "e_casscf": cas["e_casscf"],
        "dmrg_energies": energies, "s_max": s_max,
        "bqp_class": cls, "class_rationale": rationale,
        "time_budget_hit": time_budget_hit, "bond_dims_requested": bond_dims,
        "elapsed_s": elapsed_s,
        "method": "DMRG (block2, classical) convergence + entanglement diagnostic",
        "orbital_optimization_method": cas["orbital_optimization_method"],
        "provenance_source": "HPC/Laguna (DMRG classifier)",
        "hardware": detect_hardware(args.threads),
    }
    # Recorded ONLY for a real --geometry run (not --compound demo mode): these
    # are exactly the reproducibility inputs an SHCI cross-validation needs to
    # rebuild the IDENTICAL active-space Hamiltonian later (solange_shci.py) --
    # without them, "Queue SHCI Cross-Validation" against this record would have
    # nothing to copy the geometry/AVAS/charge/spin from and would need them
    # re-typed by hand, exactly the friction this field exists to remove.
    if args.geometry:
        out["geometry"] = args.geometry
        out["avas"] = args.avas
        out["charge"] = args.charge
        out["spin"] = args.spin
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
