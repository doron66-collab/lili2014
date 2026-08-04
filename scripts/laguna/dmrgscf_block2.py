#!/usr/bin/env python3
"""
dmrgscf_block2.py — DMRG as PySCF's CASSCF active-space solver, via pyblock2.

WHY THIS EXISTS
CASSCF's default active-space solver is FCI: exact, but combinatorial in the
number of active orbitals, so it walls out around 16 and takes the whole
orbital-optimization loop with it. That wall — not the downstream bond-dimension
sweep — is what stops SOLANGE's DMRG rung from ever reaching a real ~44e site.
Replacing that inner solver with DMRG makes the WHOLE procedure (orbital
selection AND the active-space solve) polynomial rather than exponential.

The standard way to do that is the `pyscf-dmrgscf` package, which shells out to
a compiled `block2main` executable. That executable does not exist in this
project's environment: what is installed is `pyblock2`, a pure-Python API around
the same block2 library, with no CLI entry point for pyscf-dmrgscf to call. This
module is the adapter that closes that gap — it implements PySCF's fcisolver
interface directly on top of `pyblock2.driver.core.DMRGDriver`, the same driver
`solange_dmrg.run_dmrg()` already uses successfully.

THE DANGER, AND WHAT IS DONE ABOUT IT
Writing this by hand means asserting a convention for the two-particle reduced
density matrix (2-RDM) — specifically which index ordering block2 returns and
how it maps onto the one PySCF's CASSCF expects. Get that wrong and nothing
crashes: CASSCF converges smoothly to a wrong number. That failure mode is
exactly what this platform exists to prevent (DP1, "verify, don't trust"), and
it is why `run_casscf()`'s existing FCI consistency gate refuses to emit a
Hamiltonian whose active-space FCI disagrees with e_casscf — a gate that
necessarily disappears past 16 orbitals, precisely where this solver is needed.

So this module carries TWO independent checks, and neither is optional:

  1. ENERGY RECONSTRUCTION (every call, every size). The RDMs are contracted
     back against the very integrals they were computed from, and the resulting
     energy is compared to the energy block2 itself reported. A wrong index
     convention breaks this identity, so a silent convention error becomes a
     loud RuntimeError. This works at any active-space size — it needs no FCI.

  2. FCI CROSS-CHECK, ELEMENT-WISE (small sizes only, and this is the one that
     matters most). Check 1 has a real blind spot, and it is not theoretical:
     for real orbitals the electron-repulsion integrals carry full 8-fold
     permutational symmetry, and on this project's own block2 build EIGHT of
     the 24 possible axis orderings reconstruct the correct energy exactly.
     Any check that compares only energies — including comparing DMRG's energy
     to FCI's — waves all eight through. `validate()` therefore compares the
     density matrices themselves, entry by entry, against FCI's: a permuted
     tensor differs element-wise even when every contraction of it agrees.
     `--dmrg-scf` runs this before it will optimize anything.

Neither check makes this module chemically validated — that still needs a
computational chemist, and the dissertation says so (§09, Future Work). What
they do is make a convention error IMPOSSIBLE TO MISS rather than invisible,
which is the difference between a demonstrable engineering artifact and a
number nobody should trust.
"""
import numpy as np


# The 2-RDM index convention this adapter assumes block2 returns, stated once,
# explicitly, so the assumption is reviewable rather than buried in a reshape.
#
#   PySCF CASSCF expects (chemist's notation, matching (pq|rs) integral order):
#       dm2[p,q,r,s] = <p† r† s q>
#       E = Σ_pq h1[p,q]·dm1[p,q] + ½ Σ_pqrs eri[p,q,r,s]·dm2[p,q,r,s] + ecore
#
# The value below was NOT taken from documentation — the documented ordering
# (0,2,1,3) was tried first and failed check 1 outright on this build. It was
# measured by `--diagnose`, which scores all 24 orderings element-wise against
# FCI's own 2-RDM on this machine's actual pyblock2. Four orderings match, and
# --diagnose verifies they map the raw array onto the identical tensor: that is
# the 2-RDM's own symmetry (dm2[p,q,r,s] = dm2[r,s,p,q] = dm2[q,p,s,r] for a
# real wavefunction), not an ambiguity, so any of the four serves.
#
# Re-run `python dmrgscf_block2.py --diagnose` after any block2 upgrade. The
# checks below are what make that instruction safe to forget: a changed
# convention fails loudly rather than producing a plausible wrong number.
_BLOCK2_TO_PYSCF_2PDM_AXES = (0, 3, 1, 2)

# Reconstructed energy must match block2's own reported energy to better than
# this (Ha). Loose enough not to trip on DMRG's own truncation/rounding at
# modest bond dimension; far tighter than any convention error could survive.
ENERGY_RECONSTRUCTION_TOL = 1e-6

# Tolerance for the element-wise DMRG-vs-FCI RDM comparison in validate(). What
# sets this is DMRG's own RDM convergence, not the precision of the wiring: on
# this build a correctly-wired CAS(6,6) lands around 3e-7 (1-RDM) and 6e-7
# (2-RDM), while the energy agrees to ~1e-11. The placement is not delicate,
# because the two outcomes are orders of magnitude apart — a wrong axis order
# does not miss by 1e-5, it misses by O(1). Anything between ~1e-8 and ~1e-3
# separates them equally well.
FCI_AGREEMENT_TOL = 1e-6


def _as_alpha_beta(nelec):
    """PySCF passes nelec as either a total count or an (nalpha, nbeta) pair."""
    if isinstance(nelec, (int, np.integer)):
        nalpha = (int(nelec) + 1) // 2
        return nalpha, int(nelec) - nalpha
    return int(nelec[0]), int(nelec[1])


def _reconstruct_energy(h1e, eri_full, dm1, dm2, ecore):
    """Contract RDMs against their own integrals — see check 1 in the module
    docstring. Any index-convention error breaks this identity."""
    e1 = np.einsum("pq,pq->", h1e, dm1)
    e2 = 0.5 * np.einsum("pqrs,pqrs->", eri_full, dm2)
    return float(e1 + e2 + ecore)


class Block2FCISolver:
    """PySCF fcisolver interface backed by pyblock2's DMRGDriver.

    Plugged in as `mc.fcisolver = Block2FCISolver(...)`, this makes CASSCF's
    inner solve a DMRG sweep instead of an FCI diagonalization. PySCF calls
    kernel() once per macro-iteration and then asks for RDMs; `civec` in that
    interface is opaque to PySCF, so this class returns a token and keeps the
    real state (the MPS and its RDMs) internally.
    """

    def __init__(self, maxM=500, scratch="./tmp_dmrgscf_orb", n_threads=4,
                 n_sweeps=10, tol=1e-8, verbose=False):
        self.maxM = int(maxM)
        self.scratch = scratch
        self.n_threads = int(n_threads)
        self.n_sweeps = int(n_sweeps)
        self.tol = float(tol)
        self.verbose = verbose
        # PySCF reads these off the solver; keep the attributes it expects.
        self.nroots = 1
        self.spin = None
        # Populated by kernel(), consumed by make_rdm1/make_rdm12.
        self._dm1 = None
        self._dm2 = None
        self._last_energy = None

    # ── PySCF fcisolver interface ──────────────────────────────────────────

    def kernel(self, h1e, eri, norb, nelec, ci0=None, ecore=0, **kwargs):
        """Solve the active space by DMRG. Returns (energy, civec-token)."""
        from pyscf import ao2mo
        from pyblock2.driver.core import DMRGDriver, SymmetryTypes

        nalpha, nbeta = _as_alpha_beta(nelec)
        n_elec, spin = nalpha + nbeta, nalpha - nbeta
        # CASSCF hands the solver symmetry-packed integrals; block2 wants the
        # full 4-index tensor. restore(1, ...) is a no-op if already full.
        eri_full = ao2mo.restore(1, np.asarray(eri), norb)
        h1e = np.asarray(h1e)

        driver = DMRGDriver(scratch=self.scratch, symm_type=SymmetryTypes.SU2,
                            n_threads=self.n_threads)
        driver.initialize_system(n_sites=norb, n_elec=n_elec, spin=spin)
        mpo = driver.get_qc_mpo(h1e=h1e, g2e=eri_full, ecore=ecore, iprint=0)
        ket = driver.get_random_mps(tag="CASCI", bond_dim=min(self.maxM, 250), nroots=1)
        energy = float(driver.dmrg(mpo, ket, n_sweeps=self.n_sweeps,
                                   bond_dims=[self.maxM],
                                   noises=[1e-5, 1e-6, 0], thrds=[self.tol] * 3,
                                   iprint=1 if self.verbose else 0))

        dm1 = np.asarray(driver.get_1pdm(ket))
        dm2 = np.asarray(driver.get_2pdm(ket)).transpose(*_BLOCK2_TO_PYSCF_2PDM_AXES)

        # ── Check 1: energy reconstruction (see module docstring) ──────────
        # Done HERE, inside every macro-iteration, not once at the end: a
        # convention error must stop the optimization rather than let it
        # converge to a confident wrong answer.
        e_check = _reconstruct_energy(h1e, eri_full, dm1, dm2, ecore)
        if abs(e_check - energy) > ENERGY_RECONSTRUCTION_TOL:
            raise RuntimeError(
                f"Block2FCISolver: RDMs do not reproduce block2's own energy "
                f"({e_check:.10f} vs {energy:.10f} Ha, Δ={abs(e_check-energy):.2e}). "
                f"This means the 2-RDM index convention assumed by this adapter "
                f"(_BLOCK2_TO_PYSCF_2PDM_AXES={_BLOCK2_TO_PYSCF_2PDM_AXES}) does not "
                f"match what this build of pyblock2 returns. REFUSING to continue — "
                f"CASSCF would otherwise converge silently to a wrong energy. Fix the "
                f"axis order in dmrgscf_block2.py and re-run validate() before using it.")

        self._dm1, self._dm2, self._last_energy = dm1, dm2, energy
        return energy, "block2-mps"      # token; PySCF never inspects it

    def make_rdm1(self, civec, norb, nelec):
        if self._dm1 is None:
            raise RuntimeError("Block2FCISolver.make_rdm1 called before kernel()")
        return self._dm1

    def make_rdm12(self, civec, norb, nelec):
        if self._dm2 is None:
            raise RuntimeError("Block2FCISolver.make_rdm12 called before kernel()")
        return self._dm1, self._dm2

    def spin_square(self, civec, norb, nelec):
        """SU2 mode is spin-adapted by construction, so the state is a pure
        spin eigenstate; report it from the requested nelec rather than
        measuring. PySCF only uses this for reporting, never for the solve."""
        nalpha, nbeta = _as_alpha_beta(nelec)
        s = 0.5 * (nalpha - nbeta)
        return s * (s + 1), 2 * s + 1


# ── Check 2: the cross-check that actually earns confidence ────────────────

def _random_test_system(norb, nelec, seed):
    """A small random Hamiltonian for the checks below.

    Deliberately random rather than a real molecule: accidental degeneracies
    and spatial symmetry in a tidy molecular system are exactly what can mask
    an index error. The two-electron array is symmetrized to carry the same
    8-fold permutational symmetry real integrals have, so the tests exercise
    the realistic (harder) case rather than an artificially easy one.
    """
    rng = np.random.default_rng(seed)
    h1e = rng.normal(size=(norb, norb))
    h1e = 0.5 * (h1e + h1e.T)                  # h1 must be symmetric
    a = rng.normal(size=(norb, norb, norb, norb))
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    eri = a + a.transpose(2, 3, 0, 1)
    return h1e, eri, _as_alpha_beta(nelec)


def validate(norb=6, nelec=6, seed=0, maxM=500, scratch="./tmp_dmrgscf_validate",
             n_threads=4, verbose=True):
    """Compare this solver's RDMs, ELEMENT BY ELEMENT, against exact FCI's.

    This is the check that catches what energy reconstruction cannot, and it
    must compare the density matrices themselves — not energies. Comparing
    energies would be no test at all here: the DMRG energy comes from block2's
    own sweep, never touching the RDMs, so it agrees with FCI regardless of how
    the 2-RDM axes are ordered. (That was a real flaw in an earlier version of
    this function; the diagnostic scan found 8 axis orderings that all
    reproduce the correct energy, and an energy-only check would have waved
    every one of them through.)

    Element-wise comparison has no such blind spot: a wrong axis order permutes
    the tensor, and a permuted tensor differs from the reference in individual
    entries even when every contraction of it agrees.

    Returns (max_dm1_error, max_dm2_error). Raises RuntimeError on mismatch.
    """
    from pyscf import fci

    h1e, eri, (nalpha, nbeta) = _random_test_system(norb, nelec, seed)
    ecore = 0.0

    # More sweeps than a production run needs. The residual here is DMRG's own
    # RDM convergence, and at the default 10 sweeps it lands ~8e-7 — under the
    # tolerance, but close enough to it that a slightly worse-converging build
    # would fail this check while being perfectly correctly wired. A false alarm
    # that blocks a working setup is its own kind of failure, so the margin is
    # bought here rather than by loosening the tolerance.
    solver = Block2FCISolver(maxM=maxM, scratch=scratch, n_threads=n_threads,
                             n_sweeps=24)
    e_dmrg, civec = solver.kernel(h1e, eri, norb, (nalpha, nbeta), ecore=ecore)
    dm1, dm2 = solver.make_rdm12(civec, norb, (nalpha, nbeta))

    fci_solver = fci.direct_spin1.FCI()
    e_fci, fci_civec = fci_solver.kernel(h1e, eri, norb, (nalpha, nbeta), ecore=ecore)
    dm1_ref, dm2_ref = fci_solver.make_rdm12(fci_civec, norb, (nalpha, nbeta))

    d_e   = abs(e_dmrg - e_fci)
    d_dm1 = float(np.max(np.abs(dm1 - dm1_ref)))
    d_dm2 = float(np.max(np.abs(dm2 - dm2_ref)))
    if verbose:
        print(f"  [validate] CAS({nelec},{norb}) random Hamiltonian · "
              f"ΔE={d_e:.2e} · max|Δ 1-RDM|={d_dm1:.2e} · max|Δ 2-RDM|={d_dm2:.2e}")

    # Energy first: if the solver itself is wrong, no RDM comparison is meaningful.
    if d_e > FCI_AGREEMENT_TOL:
        raise RuntimeError(
            f"Block2FCISolver energy disagrees with FCI: {e_dmrg:.10f} vs {e_fci:.10f} Ha "
            f"(Δ={d_e:.2e} > {FCI_AGREEMENT_TOL:.0e}) on a CAS({nelec},{norb}) random "
            f"Hamiltonian. This is a solver/parameter problem (bond dimension, sweeps, "
            f"symmetry sector), not an index convention — no axis order can repair it.")
    if d_dm1 > FCI_AGREEMENT_TOL or d_dm2 > FCI_AGREEMENT_TOL:
        raise RuntimeError(
            f"Block2FCISolver's density matrices do not match FCI's, even though the "
            f"energies agree: max|Δ 1-RDM|={d_dm1:.2e}, max|Δ 2-RDM|={d_dm2:.2e} "
            f"(tolerance {FCI_AGREEMENT_TOL:.0e}) on a CAS({nelec},{norb}) random "
            f"Hamiltonian. That combination is the signature of a wrong 2-RDM axis "
            f"order (_BLOCK2_TO_PYSCF_2PDM_AXES={_BLOCK2_TO_PYSCF_2PDM_AXES}): the "
            f"energy is a contraction and survives the permutation, the tensor does "
            f"not. Run `python dmrgscf_block2.py --diagnose` to get the right order "
            f"from this build. Do NOT use --dmrg-scf until this passes — past ~16 "
            f"orbitals there is no FCI reference left to catch it.")
    return d_dm1, d_dm2


def diagnose(norb=6, nelec=6, seed=0, maxM=500, scratch="./tmp_dmrgscf_diagnose",
             n_threads=4):
    """Determine the correct 2-RDM axis order EMPIRICALLY, from this build.

    Run this when kernel()'s energy-reconstruction check fails, which means the
    convention block2 documents is not the one it returns here. Rather than
    guess again, this asks the installed library directly: solve one small
    system, then score all 24 axis permutations of the returned 2-RDM against
    the two independent tests.

    Scoring is ELEMENT-WISE against FCI's own 2-RDM, not by reconstructed
    energy. Energy is a contraction: with the integrals' 8-fold permutational
    symmetry, many wrong orderings contract to exactly the right number — an
    earlier version of this scan scored that way and reported 8 equally-good
    candidates. Comparing the tensors entry by entry has no such degeneracy and
    normally singles out one answer.
    """
    import itertools
    from pyscf import ao2mo, fci
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes

    h1e, eri, (nalpha, nbeta) = _random_test_system(norb, nelec, seed)
    ecore = 0.0

    eri_full = ao2mo.restore(1, np.asarray(eri), norb)
    driver = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=n_threads)
    driver.initialize_system(n_sites=norb, n_elec=nalpha + nbeta, spin=nalpha - nbeta)
    mpo = driver.get_qc_mpo(h1e=h1e, g2e=eri_full, ecore=ecore, iprint=0)
    ket = driver.get_random_mps(tag="DIAG", bond_dim=min(maxM, 250), nroots=1)
    e_dmrg = float(driver.dmrg(mpo, ket, n_sweeps=10, bond_dims=[maxM],
                               noises=[1e-5, 1e-6, 0], thrds=[1e-8] * 3, iprint=0))
    dm1 = np.asarray(driver.get_1pdm(ket))
    dm2_raw = np.asarray(driver.get_2pdm(ket))

    fci_solver = fci.direct_spin1.FCI()
    e_fci, fci_civec = fci_solver.kernel(h1e, eri_full, norb, (nalpha, nbeta), ecore=ecore)
    dm1_ref, dm2_ref = fci_solver.make_rdm12(fci_civec, norb, (nalpha, nbeta))

    print(f"  DMRG energy {e_dmrg:.10f} Ha · FCI energy {float(e_fci):.10f} Ha "
          f"(Δ={abs(e_dmrg-float(e_fci)):.2e})")
    if abs(e_dmrg - float(e_fci)) > FCI_AGREEMENT_TOL:
        print("  WARNING: DMRG and FCI already disagree on the ENERGY, before any RDM is "
              "involved. That is a solver/parameter problem (bond dimension, sweeps, "
              "symmetry sector), not an index convention — fix that first; no axis order "
              "can repair it.")

    d1 = float(np.max(np.abs(dm1 - dm1_ref)))
    print(f"  1-RDM agrees with FCI to max|Δ|={d1:.2e}"
          + ("" if d1 <= FCI_AGREEMENT_TOL else "   *** 1-RDM ITSELF IS WRONG — the "
             "problem is not just the 2-RDM axis order ***"))

    print(f"  scanning all 24 axis permutations of the returned 2-RDM "
          f"(shape {dm2_raw.shape}) against FCI's, element-wise …")
    scored = []
    for perm in itertools.permutations(range(4)):
        err = float(np.max(np.abs(dm2_raw.transpose(*perm) - dm2_ref)))
        scored.append((err, perm))
        if err <= FCI_AGREEMENT_TOL:
            print(f"    {perm}  max|Δ|={err:.2e}   <-- MATCHES")
    scored.sort()
    matches = [p for e, p in scored if e <= FCI_AGREEMENT_TOL]

    print("-" * 68)
    if not matches:
        print("  NO axis permutation reproduces FCI's 2-RDM. The problem is not the index")
        print(f"  order: get_2pdm here returns something other than the spin-traced 2-RDM")
        print(f"  this adapter assumes (a spin-resolved array, a different normalization,")
        print(f"  or a different particle convention). Closest candidates, for reference:")
        for err, perm in scored[:3]:
            print(f"    {perm}  max|Δ|={err:.2e}")
        print(f"  Inspect its shape ({dm2_raw.shape}) against pyblock2's docs for this")
        print("  build before going further — do not use --dmrg-scf.")
    elif len(matches) == 1:
        print(f"  UNIQUE match: set _BLOCK2_TO_PYSCF_2PDM_AXES = {matches[0]}")
        print("  Then run the validation (no --diagnose) to confirm end to end.")
    else:
        # Several matches is the EXPECTED result, not a failure, and does not mean
        # the answer is ambiguous. The 2-RDM has its own permutational symmetry —
        # dm2[p,q,r,s] = dm2[r,s,p,q] always, and = dm2[q,p,s,r] for a real
        # wavefunction — so a set of axis orders map the raw array onto the very
        # same tensor. Verify that here rather than assert it: if every match
        # produces an identical array, any one of them is correct.
        first = dm2_raw.transpose(*matches[0])
        spread = max(float(np.max(np.abs(dm2_raw.transpose(*p) - first)))
                     for p in matches[1:])
        print(f"  {len(matches)} permutations match element-wise: {matches}")
        if spread <= FCI_AGREEMENT_TOL:
            print(f"  They produce the SAME tensor (max pairwise |Δ|={spread:.2e}) — this is")
            print("  the 2-RDM's own permutational symmetry (dm2[p,q,r,s] = dm2[r,s,p,q] =")
            print("  dm2[q,p,s,r] for a real wavefunction), not an ambiguity. Any one is")
            print("  correct.")
            print(f"  USE: _BLOCK2_TO_PYSCF_2PDM_AXES = {matches[0]}")
        else:
            print(f"  They do NOT produce the same tensor (max pairwise |Δ|={spread:.2e}),")
            print("  so this is a real ambiguity: the test state is accidentally symmetric.")
            print("  Re-run with a different --seed; whichever matches at several seeds is")
            print("  the correct one.")
    print("-" * 68)
    return matches


if __name__ == "__main__":
    import sys
    if "--diagnose" in sys.argv:
        # Empirically determine the axis order from the installed block2.
        print("=" * 68)
        print("Block2FCISolver DIAGNOSTIC — determining the 2-RDM axis convention")
        print("=" * 68)
        diagnose()
    else:
        # Standalone: `python dmrgscf_block2.py` runs the cross-check at a couple of
        # sizes. All must pass before --dmrg-scf is trustworthy on this machine.
        print("=" * 68)
        print("Block2FCISolver validation — DMRG-SCF adapter vs exact FCI")
        print("=" * 68)
        for norb, nelec in [(4, 4), (6, 6), (8, 8)]:
            validate(norb=norb, nelec=nelec, scratch=f"./tmp_dmrgscf_validate_{norb}")
        print("-" * 68)
        print("All cross-checks PASSED — the 2-RDM convention in this adapter matches "
              "this build of pyblock2.")
        print("=" * 68)
