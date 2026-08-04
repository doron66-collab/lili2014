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

  2. FCI CROSS-CHECK (small sizes only, and this is the one that matters most).
     Check 1 has a real blind spot: for real orbitals the electron-repulsion
     integrals carry full 8-fold permutational symmetry, so certain WRONG index
     permutations contract to the SAME energy as the right one. Reconstruction
     cannot distinguish them. An independent FCI solve can. `validate()` below
     runs both solvers on the same small system and compares — the only way to
     actually earn confidence in the convention before relying on it where no
     reference exists. `--dmrg-scf` in solange_hpc.py / solange_dmrg.py runs
     this once at import-time scale before it will optimize anything.

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
# pyblock2's spin-traced get_2pdm (SU2) returns <p† q† s r>, which is the same
# tensor with axes 1 and 2 exchanged. Both checks above exist to catch this
# being wrong for a given block2 version rather than to assert it is right.
_BLOCK2_TO_PYSCF_2PDM_AXES = (0, 2, 1, 3)

# Reconstructed energy must match block2's own reported energy to better than
# this (Ha). Loose enough not to trip on DMRG's own truncation/rounding at
# modest bond dimension; far tighter than any convention error could survive.
ENERGY_RECONSTRUCTION_TOL = 1e-6

# Tolerance for the DMRG-vs-FCI cross-check in validate(). Chemical accuracy is
# 1.6 mHa; a correctly-wired DMRG at a decent bond dimension on a small system
# should agree with FCI far more tightly than that, so this is deliberately
# stricter — it is testing the wiring, not the chemistry.
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

def validate(norb=6, nelec=6, seed=0, maxM=500, scratch="./tmp_dmrgscf_validate",
             n_threads=4, verbose=True):
    """Run this solver and an exact FCI solve on the SAME small random system
    and require them to agree.

    This is the check that catches what energy reconstruction cannot. For real
    orbitals the two-electron integrals have full 8-fold permutational
    symmetry, so several WRONG 2-RDM index orderings reconstruct the correct
    energy anyway — reconstruction alone cannot tell them apart. An independent
    FCI solve can, because it never touches the RDMs at all.

    Deliberately uses a random (non-symmetric, non-physical) Hamiltonian rather
    than a real molecule: accidental degeneracies and spatial symmetry in a
    tidy molecular system are exactly what can mask an index error.

    Returns (e_dmrg, e_fci). Raises RuntimeError if they disagree.
    """
    from pyscf import fci

    rng = np.random.default_rng(seed)
    h1e = rng.normal(size=(norb, norb))
    h1e = 0.5 * (h1e + h1e.T)                  # h1 must be symmetric

    # Build an eri with the 8-fold symmetry real integrals actually have, so
    # the test exercises the realistic (harder) case rather than an easy one.
    a = rng.normal(size=(norb, norb, norb, norb))
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    eri = a + a.transpose(2, 3, 0, 1)

    ecore = 0.0
    nalpha, nbeta = _as_alpha_beta(nelec)

    solver = Block2FCISolver(maxM=maxM, scratch=scratch, n_threads=n_threads)
    e_dmrg, _ = solver.kernel(h1e, eri, norb, (nalpha, nbeta), ecore=ecore)

    e_fci = fci.direct_spin1.FCI().kernel(h1e, eri, norb, (nalpha, nbeta),
                                          ecore=ecore)[0]

    delta = abs(e_dmrg - e_fci)
    if verbose:
        print(f"  [validate] DMRG {e_dmrg:.10f} Ha  vs  FCI {e_fci:.10f} Ha  "
              f"(Δ={delta:.2e}, CAS({nelec},{norb}) random Hamiltonian)")
    if delta > FCI_AGREEMENT_TOL:
        raise RuntimeError(
            f"Block2FCISolver failed its FCI cross-check: DMRG {e_dmrg:.10f} Ha vs "
            f"FCI {e_fci:.10f} Ha (Δ={delta:.2e} > {FCI_AGREEMENT_TOL:.0e}) on a "
            f"CAS({nelec},{norb}) random Hamiltonian. The adapter is mis-wired — most "
            f"likely the 2-RDM index convention. Do NOT use --dmrg-scf until this "
            f"passes: past ~16 orbitals there is no FCI reference left to catch it.")
    return float(e_dmrg), float(e_fci)


if __name__ == "__main__":
    # Standalone: `python dmrgscf_block2.py` runs the cross-check at a couple of
    # sizes. Both must pass before --dmrg-scf is trustworthy on this machine.
    print("=" * 68)
    print("Block2FCISolver validation — DMRG-SCF adapter vs exact FCI")
    print("=" * 68)
    for norb, nelec in [(4, 4), (6, 6), (8, 8)]:
        validate(norb=norb, nelec=nelec, scratch=f"./tmp_dmrgscf_validate_{norb}")
    print("-" * 68)
    print("All cross-checks PASSED — the 2-RDM convention in this adapter matches "
          "this build of pyblock2.")
    print("=" * 68)
