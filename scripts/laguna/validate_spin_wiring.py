#!/usr/bin/env python3
"""
validate_spin_wiring.py -- confirm `spin=0` passed to DMRGDriver.initialize_system
under SymmetryTypes.SU2 is actually wired through and interpreted as the target
TOTAL spin (2S), not silently ignored or interpreted as something else.

Claude Science's a.1 review (2026-09-30, round 2) accepted the methodological
argument that SU2 symmetry makes a non-singlet structurally unrepresentable when
spin=0 is requested -- so a singlet/triplet degeneracy stress test doesn't add
information. What IS worth checking, cheaply, is the WIRING: that the value 0
actually reaches initialize_system and produces a different, sensible result
than a different target spin.

Test: a tiny 4-electron/4-orbital random Hamiltonian (same style as
dmrgscf_block2.py's own _random_test_system -- deliberately random, not a real
molecule, to avoid masking a wiring error behind spatial symmetry), solved twice
under SU2: once with spin=0 (singlet target) and once with spin=2 (triplet
target, 2S=2). If the wiring is correct these two solves are physically
DIFFERENT variational problems and should generally return different energies.
If they return the identical energy, or an error, that is the signature of the
spin argument not doing what it is assumed to do.

This does not (and cannot, from here) name the definitive pyblock2 API for
reading back the state's own quantum numbers -- if pyblock2 exposes one, print
whatever ket.info/ket attributes look relevant so a human can confirm this
directly rather than trusting the energies difference alone.

Usage: python validate_spin_wiring.py
"""
import numpy as np


def _random_test_system(norb, nelec_total, seed=0):
    rng = np.random.default_rng(seed)
    h1e = rng.normal(size=(norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    a = rng.normal(size=(norb, norb, norb, norb))
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    eri = a + a.transpose(2, 3, 0, 1)
    return h1e, eri


def solve(norb, nelec, spin, scratch):
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    h1e, eri = _random_test_system(norb, nelec)
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=2,
                     stack_mem=int(1.0 * (1 << 30)))
    drv.initialize_system(n_sites=norb, n_elec=nelec, spin=spin)
    mpo = drv.get_qc_mpo(h1e=h1e, g2e=eri, ecore=0.0, iprint=0)
    ket = drv.get_random_mps(tag=f"KET_s{spin}", bond_dim=100, nroots=1)
    e = float(drv.dmrg(mpo, ket, n_sweeps=20, bond_dims=[100],
                       noises=[1e-4, 1e-5, 0], thrds=[1e-10] * 3, iprint=0))
    interesting = [a for a in dir(ket) if 'info' in a.lower() or 'target' in a.lower()
                   or 'quanta' in a.lower() or 'spin' in a.lower()]
    print(f"  spin={spin}: E={e:.10f}   ket attrs of interest: {interesting}")
    for a in interesting:
        try:
            print(f"    ket.{a} = {getattr(ket, a)}")
        except Exception as exc:
            print(f"    ket.{a} -> could not read ({exc})")
    return e


def main():
    norb, nelec = 4, 4  # same electron count, same orbitals -- only spin differs
    print("Solving the SAME random Hamiltonian at spin=0 (singlet target) and "
          "spin=2 (triplet target, 2S=2). Under correct SU2 wiring these are "
          "different variational problems.")
    e0 = solve(norb, nelec, 0, "./tmp_spinwire_s0")
    e2 = solve(norb, nelec, 2, "./tmp_spinwire_s2")
    print("=" * 60)
    print(f"E(spin=0) = {e0:.10f}")
    print(f"E(spin=2) = {e2:.10f}")
    diff = abs(e0 - e2)
    print(f"|difference| = {diff:.3e} Ha")
    if diff < 1e-8:
        print("*** SUSPICIOUS: the two solves returned (numerically) the SAME "
              "energy. Either this random Hamiltonian is spin-degenerate by "
              "accident (re-run with a different seed to check) or the `spin` "
              "argument is not being wired through as expected -- do not treat "
              "solange_dmrg.py's SU2 singlet assumption as confirmed on this "
              "build until re-run with another seed rules out coincidence.")
    else:
        print("Energies differ as expected -- spin=0 vs spin=2 solve genuinely "
              "different problems on this pyblock2 build, consistent with `spin` "
              "being interpreted as the target 2S under SU2 symmetry.")


if __name__ == "__main__":
    main()
