# Pre-registration — SDHB S3 [3Fe-4S], DMRG tractability measurement

**Frozen:** 2026-10-10, before any DMRG run on this target.
**Target:** SDHB (UniProt P21912) centre S3, the [3Fe-4S]¹⁺ cluster of PDB 8GS8
chain B, iron ligated by Cys196 / Cys243 / Cys249 (assigned geometrically,
Fe–S < 2.8 Å; the ligation is cross-wired relative to the sequence motifs).
**Clinical anchor:** C196Y, PPGL4, dbSNP rs876658367 — a pathogenic missense
directly at an S3 ligand. The wild-type site is what is measured; the variant is
what makes the site worth measuring.

---

## 1. Model

| Item | Value |
|---|---|
| Structure | 8GS8 chain B, cryo-EM 2.86 Å |
| Core geometry | used as deposited — Fe–Fe 2.699 Å vs 2.694 ± 0.035 Å across 24 matched-ligation X-ray reference sites ≤ 1.2 Å (+0.1 σ) |
| QM region | backbone-extended model (job 2350825): cluster + 3 ligating Cys + backbone (N, CA, C, O, CB) of SER198, TYR199, HIS244, THR245, ILE246, MET247, ASN248, THR250 |
| Why extended | the 43-atom model omits 10 of the 12 backbone N–H···S donors and sits in a pocket whose protein environment is net +2 within 10 Å while the cluster carries −2 |
| Charge / spin | −2 ; 2S = 15 (S = 15/2, ferromagnetic) and 2S = 1 (S = 1/2, ground) |
| Basis | def2-TZVP on Fe and S; def2-SVP on the added backbone shell |
| Solvation | ddCOSMO, ε = 4 |
| Orbitals | ROKS/BP86 on the **high-spin** state, **one orbital set used for both spin legs** |

**Orbital-source rationale.** At the HF level the ligand-radical state lies
288.56 mHa *below* the intended Fe(III)₃ state (gas-phase HF energy of the BP86
determinant −7308.67732798 vs the ROHF ligand-radical minimum −7308.96589), so
HF orders the two states backwards. BP86 generates orbitals only; DMRG
recomputes the energy inside the active space, so the functional does not need
to be accurate for spin-state energetics.

**Shared orbitals are load-bearing, not a convenience.** Per-leg optimised
orbitals absorb part of the antiferromagnetic state's difficulty into the
orbital rotation and shrink any spin-resolved ratio. They also require an SCF at
2S = 1, which is multireference and will land on broken symmetry.

---

## 2. State gate (must pass before AVAS)

Applied to meta-Löwdin spin populations. **Spin populations are not observables**
— they are partition-dependent — so the gate is written in scale-free terms
wherever possible.

| Gate | Bound | Measured (43-atom, BP86/ddCOSMO) |
|---|---|---|
| `max(ligand spin) / min(Fe spin)` | ≤ 0.35 | 0.245 |
| `Σ spin(Fe) / 2S` | 0.60 – 0.85 | 0.691 |
| Fe–Fe spin spread | ≤ 0.5 | 0.022 |
| Σ\|spin\| on non-Fe, non-S atoms | ≤ 0.3 | 0.253 |
| Sign | all Fe same sign | pass |
| Total spin population | == 2S exactly | 15.0000 |

**Reported, not gated:** all seven sulfurs, partitioned sulfide / thiolate.
Measured: sulfides 0.774 / 0.781 / 0.821 / 0.845, thiolates 0.373 / 0.374 /
0.408. The four sulfides are empirically **one class** (the μ₃ value 0.781 sits
inside the μ₂ range) — the μ₂/μ₃ split proposed on 2026-10-10 is withdrawn, see
§6. The thiolates are a separate class at roughly half the sulfide value.

**Scope note.** The Fe-spread gate is valid for S3 and S1, where all irons are
formally equivalent. It is **not** valid for S2 [4Fe-4S]²⁺, which is
mixed-valence Fe(III)₂Fe(II)₂ and where a spread is expected.

---

## 3. Active space — FROZEN

**Target AO set:** Fe 3d (15) + Fe 3d′ correlating shell (15) + bridging-sulfide
S 3p (12) = 42 target AOs. Thiolate S 3p are **excluded**; the measured spin
populations support this (thiolates carry roughly half what the sulfides carry).

**AVAS threshold: 0.2, pre-set, not tuned.** AVAS may return more orbitals than
target AOs — in a covalent Fe–S cluster both the bonding and the antibonding
Fe–S partner can clear the threshold. That is expected and is not a defect.
`openshell_option=3` (every SOMO kept active).

**The 3d′ shell is not optional, and this is why:**

| Space | 2S | nα | nβ | α virtuals | filling | exact Schmidt rank | M=4000 captures |
|---|---|---|---|---|---|---|---|
| CAS(49,32) | 15 | **32** | 17 | **0** | 76.6 % | 5.27e4 | 7.6 % |
| CAS(49,32) | 1 | 25 | 24 | 7 | 76.6 % | 9.03e6 | 0.044 % |
| CAS(49,47) | 15 | 32 | 17 | **15** | **52.1 %** | 3.64e12 | 1.1e-9 |
| CAS(49,47) | 1 | 25 | 24 | 22 | 52.1 % | 7.00e13 | 5.7e-11 |

In CAS(49,32) at 2S = 15, **nα = ncas**: every α spin-orbital is occupied and the
α string space has dimension exactly 1. The two legs then differ by **171× in
exact Schmidt rank** for a reason that has nothing to do with spin coupling.
Adding the 3d′ shell takes the filling to 52.1 % — matching TP53 C275F's 52.9 %,
the only prior space that was not starved — and reduces the leg asymmetry to
19.2×. `ncore` is unchanged at 130, because the added shell is empty.

**Implementation** (either route, recorded in provenance): augment the AVAS
output with the 15 lowest virtual orbitals ranked by Fe-d projection onto the
*full-basis* Fe d AOs (5 per iron, count pre-registered), **or** run AVAS against
a minimal reference carrying two d shells for Fe. Neither has been tested here.

**Recorded as an alternative, not the reported space:** CAS(57,39), Fe 3d + all
S 3p (36 target AOs), which is the space the earlier `"S 3p"` AVAS label
produced by also matching the three thiolates.

---

## 4. Metric — FROZEN

### 4.1 Primary

```
R = dE(SDHB S3, 2S=1, M -> 2M) / dE(control, same nelecas, same ncas, 2S=1, same M -> 2M)
```

This is the same R validated on TP53 C275F (7.20 wild type, 7.43 mutant), moved
to an odd-electron system. It is measured on the **S = 1/2 leg only** — the leg
that carries the antiferromagnetic coupling.

### 4.2 Why the spin-resolved ratio is demoted

`R_spin = dE(2S=1)/dE(2S=15)` was the design of 2026-10-07. It is demoted to a
reported companion because the two legs differ combinatorially at fixed
(nelecas, ncas): 171× in exact Schmidt rank at CAS(49,32), 19.2× at CAS(49,47).
The asymmetry is **irreducible** — a spin-aligned determinant genuinely has fewer
configurations — and it cannot be cancelled by an external control, because a
chemistry-free system with 15 unpaired electrons is 15 non-interacting radical
centres, for which dE → 0 in both legs and the ratio is 0/0.

`R_spin` is therefore **reported alongside its analytic combinatorial floor**
(the Schmidt-rank ratio at the space actually used), and no verdict is drawn
from it alone.

### 4.3 The S = 15/2 leg keeps three jobs

1. Orbital generation (§1).
2. The exchange coupling J, with an Fe–Fe sensitivity scan over the reference
   range 2.522–2.744 Å.
3. The state gate (§2).

---

## 5. Control — FROZEN

| Requirement | Value |
|---|---|
| System | saturated alkyl radical, C 2p target AOs — the same chemistry-free family as the n-C12H26 control used for TP53 C275F |
| Match | identical `(nelecas, ncas, 2S)` — **2S must match**, because the Hilbert-space structure depends on it |
| Size | forced to exactly the (nelecas, ncas) the target returns; at CAS(49,47) this is ≈ 16 carbons (C 2p on 16 C = 48 orbitals) |
| Basis | def2-TZVP, same family as the target |
| Protocol | identical M ladder, sweep count and noise schedule |
| Orbitals | RHF/ROHF is acceptable here — the control is weakly correlated, so HF and BP86 orbitals are near-identical; the difference is recorded, not corrected |
| Replicates | **3 independent controls** (distinct conformers / seeds), to measure the control's own noise band |

Matching size and filling alone is **not** sufficient. The control must match
2S, because nα and nβ — and therefore the number of α virtuals and the exact
Schmidt rank — are set by it.

---

## 6. Corrections carried into this document

- **The μ₂/μ₃ sulfide gate split is withdrawn.** It was justified on the grounds
  that a μ₃ sulfide receives covalent spin from three irons and a μ₂ from two,
  so the μ₃ should carry more. The measurement says otherwise: μ₃ = 0.781 sits
  inside the μ₂ range (0.774–0.845), and the Fe–S bond lengths are identical
  across all four sulfides (2.254–2.261 Å, spread 0.007 Å), so geometry does not
  explain it either. A sulfide has three 3p orbitals; a μ₂ can commit close to a
  full p orbital per bond while a μ₃ must share three across three bonds, so
  per-bond covalency falls as the bridging count rises and the totals come out
  comparable. The spin on a bridging sulfide is set by its p-orbital budget, not
  by how many metals it touches.
- **The absolute per-sulfur bound (0.2–0.8) is withdrawn** — wrong in form:
  absolute, multi-class, and partition-dependent. Replaced by §2.
- **CAS(53,38) has no provenance.** The specification issued on 2026-10-07 was
  CAS(39,27) (15 Fe 3d + 12 bridging S 3p). CAS(53,38) does not appear in the
  advisory record and should not be attributed to it.

---

## 7. Threshold — deliberately NOT frozen

**No R threshold is set, and none should be.**

R is reported descriptively, with an uncertainty band propagated from the three
control replicates. A threshold is a governance decision about what "Class A"
means — it belongs to the investigator, not to the chemistry advice — and any
threshold calibrated on the first metal system ever run would be calibrated on
a target rather than on controls. That ordering is what required five prior
retractions.

What *is* frozen is the condition under which R may be interpreted at all:

```
R is interpretable only if  |R - 1| > k * sigma_R
   sigma_R propagated from the spread of the 3 control replicates
   k calibrated from CONTROLS ONLY, never from targets, and set in a later
     document once more than one control family has been run
```

Outcome space: `HARDER THAN CONTROL` / `INDISTINGUISHABLE FROM CONTROL` /
`INCONCLUSIVE` (ladder did not meet its own convergence criteria).

---

## 8. Open prediction, already on record

For job 2350825 (backbone-extended model, same BP86 functional), fixed before
that job finished:

```
|delta spin(S3_mu3)| < |delta spin(S1)|, |delta spin(S2)|, |delta spin(S4)|
and     mean spin(Fe) rises above 3.457
```

The μ₃ sulfide is buried at the cuboid centre and gains **zero** new backbone
N–H donors, while S4 gains 2 (ASN248, MET247), S2 gains 2 (ILE246, MET247) and
S1 gains 1 (HIS244). This is an ordering prediction requiring no calibration.
If S3 moves like the others, the hydrogen-bond diagnosis is wrong and the cause
is the functional or the solvation model.

---

## 9. What this document does not establish

- No DMRG has been run on this target. R has not been measured.
- The 3d′ augmentation has not been tested; the two implementation routes in §3
  are untested recommendations.
- `M_FEASIBLE` for this system has not been measured from wall-clock.
- The claim that per-bond covalency falls with bridging count (§6) is an
  explanation offered after the fact for a measurement that contradicted the
  prior reasoning. It is not independently tested here.
