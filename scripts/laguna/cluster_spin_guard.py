"""
cluster_spin_guard.py  --  preflight invariants for open-shell metal-cluster QM regions.

Written after the SDHB S3 [3Fe-4S] parity crash (2026-10-07).

The crash was NOT a code bug: the three ligating cysteines had been built as
neutral thiols (-SH) instead of thiolates (-S-), adding one spurious hydrogen
per residue.  With an ODD number of ligating cysteines that flips the electron
parity and PySCF refuses to build the molecule.  With an EVEN number it does
not -- the same bug on SDHB S1 or S2 (four cysteines each) passes every check
the pipeline currently has and runs to completion with FOUR extra electrons
silently injected into the cluster.

The existing charge check could not catch this, because it compared
EXPECTED_NET_CHARGE (computed from the formal-charge recipe: core + n*(-1))
against the charge passed to gto.M() (the same recipe).  It compared a constant
to itself.  Every function below derives its quantity from the ATOM LIST.
"""

Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "FE": 26}

# formal d-electron count and maximum spin alignment, by iron oxidation state
_D_ELECTRONS = {2: 6, 3: 5, 4: 4}
_S_TIMES_TWO = {2: 4, 3: 5, 4: 4}


def sum_nuclear_charge(symbols):
    """Total nuclear charge of the built region, from the atoms themselves."""
    return sum(Z[s.upper()] for s in symbols)


def assert_spin_parity(nelectron, spin, label=""):
    """
    N_alpha = (N + spin)/2 must be an integer, so N and spin share parity.

    This is the ONE check that catches any odd-count error in the build --
    a missing cap, a spurious hydrogen, a charge wrong by an odd amount --
    without knowing anything about the chemistry.  It is cheap and it is
    unconditional.  Run it on every open-shell target, not only the ones
    with an odd number of ligands.
    """
    if (nelectron - spin) % 2 != 0:
        raise AssertionError(
            "{}electron/spin parity conflict: N={} ({}), spin={} ({}). "
            "Something in the built region is off by an odd number -- check "
            "ligand protonation first, then cap hydrogens, then the charge."
            .format(label + ": " if label else "",
                    nelectron, "even" if nelectron % 2 == 0 else "odd",
                    spin, "even" if spin % 2 == 0 else "odd"))


def assert_ligands_deprotonated(coords, symbols, sg_indices, cutoff=1.6):
    """
    A cysteine that ligates an Fe-S cluster is a THIOLATE.  A protonated SG
    does not coordinate, and the extra proton corrupts both the electron count
    and the formal charge.  Standard protonation tools add HG to every CYS
    unless the residue is typed as metal-coordinating (CYM), so this is the
    default failure mode, not an exotic one.

    coords: (natm, 3) array in Angstrom; sg_indices: indices of the ligating S.
    """
    import numpy as np
    coords = np.asarray(coords, dtype=float)
    offenders = []
    for i in sg_indices:
        for j, s in enumerate(symbols):
            if s.upper() == "H" and np.linalg.norm(coords[i] - coords[j]) < cutoff:
                offenders.append((i, j))
    if offenders:
        raise AssertionError(
            "ligating sulfur(s) carry hydrogen: {} -- these cysteines were "
            "built as thiols, not thiolates. Remove the S-H hydrogens; the net "
            "charge does NOT need a compensating change, because the correct "
            "charge was already specified for the thiolate form."
            .format(offenders))


def charge_from_formal_structure(core_charge, n_thiolate):
    """Formal charge of [cluster]^core_charge capped by n thiolates."""
    return core_charge - n_thiolate


def assert_charge_consistent(symbols, nelectron, charge, label=""):
    """Non-circular: the atoms must actually account for the electrons."""
    expected = sum_nuclear_charge(symbols) - charge
    if expected != nelectron:
        raise AssertionError(
            "{}atom list and electron count disagree: sum(Z)={}, charge={:+d} "
            "=> N should be {}, but the molecule reports {}."
            .format(label + ": " if label else "", sum_nuclear_charge(symbols),
                    charge, expected, nelectron))


def ferromagnetic_spin(iron_oxidation_states):
    """
    PySCF `spin` (= 2S) for the fully spin-aligned reference, and the d-electron
    count, from the formal oxidation states.  Both are locked by the core charge
    -- they are not free parameters to reconcile a parity conflict with.
    """
    d = sum(_D_ELECTRONS[ox] for ox in iron_oxidation_states)
    two_s = sum(_S_TIMES_TWO[ox] for ox in iron_oxidation_states)
    core_charge = sum(iron_oxidation_states) - 2 * 0  # caller adds the sulfide term
    return {"d_electrons": d, "spin": two_s, "S": two_s / 2.0,
            "s_squared": (two_s / 2.0) * (two_s / 2.0 + 1.0),
            "sum_fe_charge": core_charge}


def preflight(symbols, coords, nelectron, charge, spin, sg_indices,
              iron_oxidation_states, core_charge, label=""):
    """Run every invariant. Call this before gto.M(), not after the SCF."""
    assert_charge_consistent(symbols, nelectron, charge, label)
    assert_ligands_deprotonated(coords, symbols, sg_indices)
    assert_spin_parity(nelectron, spin, label)
    fm = ferromagnetic_spin(iron_oxidation_states)
    n_thiolate = len(sg_indices)
    q = charge_from_formal_structure(core_charge, n_thiolate)
    if q != charge:
        raise AssertionError(
            "{}formal charge {:+d} (core {:+d} + {} thiolates) != charge passed "
            "to the molecule ({:+d})".format(label + ": " if label else "",
                                             q, core_charge, n_thiolate, charge))
    if spin not in (0, 1, fm["spin"]):
        raise AssertionError(
            "{}spin={} is neither the ground state (0 or 1) nor the "
            "ferromagnetic reference ({})".format(label + ": " if label else "",
                                                  spin, fm["spin"]))
    return {"n_electrons": nelectron, "charge": charge, "spin": spin,
            "expected_s_squared": (spin / 2.0) * (spin / 2.0 + 1.0), **fm}
