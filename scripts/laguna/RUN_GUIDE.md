# SOLANGE — Run Guide (HPC · DMRG · Quantum)

A clean, ordered runbook for launching simulations from SOLANGE onto **Laguna
(USC CARC)** and **real IBM quantum hardware**. Work top to bottom. Every command
is copy-and-run — no placeholders to edit (the one-time setup handles that).

The platform is organized as an **escalation ladder** (small → large → quantum):

| Rung | What | Where it runs | How you launch it |
|------|------|---------------|-------------------|
| 1 · Laptop | in-browser VQE | your browser | click **Run Live Simulation** in SOLANGE |
| 2 · HPC | classical CASSCF/VQE | Laguna | **queue** from UI → agent pulls |
| 3 · DMRG | A/B/C classifier | Laguna | **copy command** from UI → run |
| 4 · Quantum | real IBM Heron QPU | IBM cloud | **queue** from UI → QPU agent pulls |

Runs come back sealed and notarized by **LEON** into Rung 2/3/4 of the
**Orchestration** tab.

---

## 0. Every new terminal — one word

Open a terminal on Laguna (JupyterLab: **+** → Terminal), then:

```bash
solange
```

That loads conda, activates the `base` env, and `cd`s into `~/lili2014`. Verify
once if you like:

```bash
python -c "import block2, pyscf; print('env OK')"
```

Keep the repo current when a guide step says so:

```bash
git pull origin claude/code-access-clarification-ab1W8
```

---

## 1. HPC classical runs — Rung 2 (queue + agent)

The classical Laguna agent pulls jobs you queue in the UI. It runs **only** while
alive in your session; it claims **only** what you queued.

### 1a. Start the agent (once — survives closing the terminal)

```bash
bash scripts/laguna/agent_keepalive.sh start
```

Watch it come alive; the SOLANGE agent dot turns green within ~15 s:

```bash
tail -f ~/.solange/agent.log      # Ctrl+C exits the tail — it does NOT stop the agent
```

Manage it any time:

```bash
bash scripts/laguna/agent_keepalive.sh status
bash scripts/laguna/agent_keepalive.sh stop
```

### 1b. Queue a job from SOLANGE

1. Load an NGS report so the **Per-Mutation Routing** list appears.
2. Tick a mutation.
3. (optional) tick **run VQE (≤20q)** — only for CAS(8,8)/CAS(10,10).
4. Click **▶ Dispatch selected to HPC**.

The agent pulls it within ~15 s → row appears in **Rung 2**.

### 1c. Run directly, without the queue (alternative)

```bash
bash scripts/laguna/make_row.sh          # real run → creates a Rung 2 row
bash scripts/laguna/make_row.sh --dry     # 5-sec sanity check, does NOT touch SOLANGE
```

---

## 2. DMRG runs — Rung 3 (copy command)

DMRG is a classical, largemem Laguna run. It is launched by a copy-command
(not a queue).

### 2a. From SOLANGE

1. Tick a mutation → click **⬇ DMRG command** → **Copy**.
2. Paste into your terminal → Enter.

### 2b. The command (what the button gives you)

```bash
bash scripts/laguna/run_dmrg.sh --key TP53_C275F --side native --ncas 8 --nelecas 8 --bond-dims 250,500,1000,2000 --submit
```

`run_dmrg.sh` auto-detects your env and sets up the block2/MKL libraries — no
setup needed. On success: `db=stored` → row in **Rung 3**, notarized by LEON.

- Model compound is auto-resolved from the key.
- The sweep **early-stops** once the energy converges (nothing lost — the verdict
  can't change). Use `--no-early-stop` to force the full sweep.
- **Real functional site** (past the classical wall → possible Class A):
  ```bash
  bash scripts/laguna/run_dmrg.sh --key <KEY> --geometry site.xyz --avas "Zn 3d,S 3p" --submit
  ```

### 2c. SHCI — a second, independent classical classifier (Rung 3)

SHCI reaches its own A/B/C verdict (convergence-only signal — see
`solange_shci.py`'s module docstring for what it does and does not check yet).
It does NOT require a DMRG record to run against — standalone classification
and DMRG cross-validation are two independent uses of the same script.

One-time: create the table (Supabase SQL editor) before the first submit —
run both blocks even on a table that already exists from before this was
added (the `alter table` lines are additive and safe to re-run):

```sql
create table if not exists public.shci_crossvalidations (
  id uuid primary key,
  created_at timestamptz not null default now(),
  key text, dmrg_classification_id uuid,
  ncas int, nelec int, e_shci numeric,
  sweep_eps text, method text, elapsed_s numeric,
  provenance_source text, hardware text,
  e_dmrg_ref numeric, delta_mha numeric, agreement boolean,
  shci_seal_payload text, shci_hash text
);
alter table public.shci_crossvalidations add column if not exists bqp_class text;
alter table public.shci_crossvalidations add column if not exists class_rationale text;
alter table public.shci_crossvalidations add column if not exists shci_energies jsonb;
```

Needs a Dice build (see the SHCI/Dice build notes — Boost/Eigen/HDF5 modules,
`make -j4 Dice EIGEN=... HDF5=... BOOST=...`).

**Standalone classification** — no DMRG record needed or referenced:

```bash
python3 scripts/laguna/solange_shci.py --geometry cluster.xyz --charge <C> --spin <S> \
  --basis sto-3g --avas "<AVAS criterion>" --key <KEY> \
  --dice-scripts ~/lili2014/Dice/scripts --sweep-eps 1e-2,1e-3,5e-4,1e-4 --submit
```

**Cross-validated against an existing DMRG classification** — add
`--dmrg-classification-id`, the DMRG record's own `id` from its
`/hpc/dmrg/submit` response. The backend computes `delta_mha`/`agreement`
itself from that record's own stored energy at ingestion — never from this
script's own claim (DP1, verify-don't-trust) — and refuses (409) if the active
space doesn't match the DMRG record's exactly, or if that record doesn't
exist (404). Either way, SHCI's own `bqp_class` verdict is recorded and shown
regardless of whether a DMRG record was named:

```bash
python3 scripts/laguna/solange_shci.py --geometry <SAME xyz the DMRG run used> \
  --charge <SAME> --spin <SAME> --basis <SAME> --avas "<SAME AVAS as the DMRG run>" \
  --key <KEY> --dmrg-classification-id <uuid-from-dmrg-submit> \
  --dice-scripts ~/lili2014/Dice/scripts --sweep-eps 1e-3,5e-4,1e-4 --submit
```

**IMPORTANT — orbital basis must match, or the comparison is not apples-to-apples.**
The `--ncas`/`--nelecas` check above only confirms the active-space *size* matches;
it says nothing about whether the two methods solved on the *same orbitals*. The
command above builds SHCI's own fixed, unrotated AVAS orbitals — a fair
comparison only if the referenced DMRG record was run with `--casci` (no
orbital optimization of its own). If it was run with `--dmrg-scf` (required
above ~16 active orbitals — i.e. every real CAS(100+,60+) run so far), its
orbitals were CASSCF-rotated, a genuinely different basis, and a raw energy
delta between the two is partly an artifact of that mismatch, not a real
solver disagreement (found live 2026-09-22 on TP53_R175_NATIVE — Δ=42.276 mHa,
fully explainable by the basis mismatch alone, with no way from the numbers
alone to tell how much, if any, was real). solange_dmrg.py's `--geometry` path
now saves the exact rotated orbitals it solved on to `<its --scratch>/mo_coeff_final.npy`
— pass that file via `--orbitals` (with `--ncas`/`--nelecas` set explicitly,
since AVAS does not run in this mode) to make SHCI solve on that SAME basis:

```bash
python3 scripts/laguna/solange_shci.py --geometry <SAME xyz the DMRG run used> \
  --charge <SAME> --spin <SAME> --basis <SAME> --avas "<SAME AVAS as the DMRG run>" \
  --orbitals <that DMRG run's --scratch>/mo_coeff_final.npy \
  --ncas <SAME as the DMRG record> --nelecas <SAME as the DMRG record> \
  --key <KEY> --dmrg-classification-id <uuid-from-dmrg-submit> \
  --dice-scripts ~/lili2014/Dice/scripts --sweep-eps 1e-3,5e-4,1e-4 --submit
```

### 2d. Queue SHCI from the browser — "▶ Queue SHCI" button (Rung 3)

A DMRG record with stored geometry/AVAS AND its own saved orbitals
(`orbitals_path` — every real `--geometry` run submitted since 2026-09-22)
shows a **▶ Queue SHCI** button on its own row in the Rung 3 table. Clicking
it copies that record's own geometry/AVAS/charge/spin/orbitals — nothing
re-typed — and queues an SHCI job for the SAME classical agent that already
runs HPC/DMRG jobs (`solange_hpc.py --agent`, started once via
`agent_keepalive.sh start`). No terminal interaction needed per job once that
agent is running, and no need to think about which orbital basis to match —
`orbitals_path` being present is exactly what makes the SAME basis get used
automatically (`--orbitals`, under the hood; see §2c's cross-validation note
below for why this matters). A record with geometry/AVAS but no
`orbitals_path` (predates 2026-09-22) shows **"SHCI needs re-run"** instead —
its rotated orbitals, if any, were never saved and cannot be recovered
after the fact; re-run DMRG to get a matchable record.

One-time migrations (Supabase SQL editor) — additive, safe to re-run:

```sql
-- dmrg_classifications: store the reproducibility inputs a later SHCI
-- cross-validation needs to copy (only populated for a real --geometry run,
-- not --compound demo mode).
alter table public.dmrg_classifications add column if not exists geometry text;
alter table public.dmrg_classifications add column if not exists avas text;
alter table public.dmrg_classifications add column if not exists charge int;
alter table public.dmrg_classifications add column if not exists spin int;
-- avas_threshold: the AVAS projection-score threshold actually used (pyscf
-- default 0.2 unless explicitly widened via solange_dmrg.py's --avas-threshold).
alter table public.dmrg_classifications add column if not exists avas_threshold numeric;

-- orbital_optimization_converged: whether CASSCF/DMRG-SCF's own orbital
-- optimization actually converged (mc.converged), independent of whether
-- the subsequent DMRG bond-dimension sweep converged. Found live 2026-09-12:
-- a run whose orbitals hit max_cycle_macro without converging (not a time-
-- budget stop) was previously indistinguishable from a fully-converged run
-- in the stored record, even though the console log itself said "CASSCF did
-- NOT converge" - the classification could not honestly be called FINAL.
alter table public.dmrg_classifications add column if not exists orbital_optimization_converged boolean;

-- orbital_optimization_method: e.g. "DMRG-SCF (block2, maxM=250)" vs "CASCI,
-- fixed AVAS orbitals (no optimization) + ...". solange_dmrg.py always
-- submitted this; it was never in the backend's DB whitelist either (found
-- live 2026-09-22 auditing two historical SHCI cross-validations and finding
-- the field simply absent from both stored records).
alter table public.dmrg_classifications add column if not exists orbital_optimization_method text;

-- orbitals_path: a LOCAL Laguna filesystem path to the .npy mo_coeff matrix
-- h1e/h2e were actually built from (solange_dmrg.py's --geometry path saves
-- this since 2026-09-22). Lets "Queue SHCI" pass --orbitals automatically so
-- SHCI solves on the SAME basis a --dmrg-scf DMRG record actually used,
-- instead of building its own unrotated AVAS orbitals — the root cause of a
-- spurious Class-A "disagreement" found live on TP53_R175_NATIVE (Δ=42.276
-- mHa, fully explainable by the basis mismatch alone). NOT portable content
-- like geometry/avas/charge/spin below — meaningless off this cluster, but
-- that's fine since the whole SHCI cross-validation flow already runs here.
alter table public.dmrg_classifications add column if not exists orbitals_path text;

-- hpc_dispatch: carries an SHCI job's parameters from the "Queue SHCI" button
-- to the agent that picks it up.
alter table public.hpc_dispatch add column if not exists geometry text;
alter table public.hpc_dispatch add column if not exists avas text;
alter table public.hpc_dispatch add column if not exists charge int;
alter table public.hpc_dispatch add column if not exists spin int;
alter table public.hpc_dispatch add column if not exists sweep_eps text;
alter table public.hpc_dispatch add column if not exists dmrg_classification_id uuid;
-- orbitals_path: mirrors dmrg_classifications' own column above — carries the
-- referenced DMRG record's saved mo_coeff_final.npy path through to the agent.
alter table public.hpc_dispatch add column if not exists orbitals_path text;

-- hpc_dispatch: opt-in Zero-Noise Extrapolation for a Rung 4 (QPU) job — set
-- per-job, since it genuinely multiplies QPU time/cost (unlike the always-on
-- TREX readout-error mitigation, which only needs extra calibration circuits).
alter table public.hpc_dispatch add column if not exists zne boolean;
```

The agent needs `boost/1.85.0` loadable via `module load` for Dice's shared
library — the dispatch path wraps the SHCI subprocess call in `bash -lc` so
the `module` shell function is available even if the agent itself was started
without it (`RUN_GUIDE.md`'s own earlier troubleshooting note for this exact
`libboost_mpi.so.1.85.0` error applies here too, just handled automatically
instead of by hand).

### 2e. "The Classifier" — PDB to a two-method decision, zero terminal (Rung 3)

The block at the top of the Rung 3 card (PDB ID, chain, residue #, expected
residue, target key, optional AVAS, radii, max orbitals, spin) queues
`job_type='screen_classify'`. The agent runs `solange_screen_and_classify.py`,
which does two parts:

**Part 1 — cluster acquisition.** Checks `GET /api/simulate/cluster/lookup`
for this exact (pdb_id, chain, resi) first — if this site was ever built
before, its cached geometry/AVAS/charge/spin/active-space are reused
directly, no re-protonation, no re-carving, no re-probing. Otherwise:
`protonate.py` (fetch + add missing hydrogens) → carve at each radius with
`build_qm_cluster.py` and probe with `avas_probe.py`, shrinking the radius
until the active space fits `--max-orbitals` → the accepted cluster is saved
via `POST /api/simulate/cluster/save` so the next run of this same site
skips straight to Part 2.

**Part 2 — classify with both methods, mandatory.** `run_dmrg.sh --submit`,
then `solange_shci.py --submit` cross-validated against the DMRG record Part
2 just created. SHCI is NOT optional here — the classifier's decision is the
two-method outcome, not one method's opinion (pass `--skip-shci` directly to
the script only for local testing without a Dice build). The job's final
output states the outcome explicitly: both classes, and whether their
energies agree.

Nobody needs to open a terminal for this — provided the classical agent is
already running (§1a).

**What it still does not decide for you** (same limits already stated in
`screen_target.py` and `gate2_requirements.py`, now inherited rather than
re-litigated): the AVAS chemical criterion (a generic default is used unless
you supply one — and that default has already been shown, in this project's
own investigation, to pick up unintended atoms on a protein-heavy cluster),
protonation pH (fixed at 7.0), and a metal-containing site's spin/oxidation
state (spin is a plain dropdown here, not a chemist's sign-off — Gate 2,
tracked separately above, does not gate this pipeline yet).

One-time migrations (Supabase SQL editor), additive:

```sql
alter table public.hpc_dispatch add column if not exists pdb_id text;
alter table public.hpc_dispatch add column if not exists chain text;
alter table public.hpc_dispatch add column if not exists resi int;
alter table public.hpc_dispatch add column if not exists expect_resname text;
alter table public.hpc_dispatch add column if not exists radii text;
alter table public.hpc_dispatch add column if not exists max_orbitals int;
alter table public.hpc_dispatch add column if not exists skip_shci boolean;

create table if not exists public.qm_clusters (
  id uuid primary key,
  created_at timestamptz not null default now(),
  pdb_id text, chain text, resi int, expect_resname text,
  key text, geometry text, avas text, charge int, spin int,
  ncas int, nelec int, radius numeric,
  unique(pdb_id, chain, resi)
);
```

(`avas`, `spin`, `sweep_eps` were already added in §2d's migration and are
reused here.)

A failure at any step (protonation, no radius fitting `max_orbitals`, DMRG,
or SHCI) stops the whole job and reports why in its status note — it does not
retry with different parameters or guess a fix.

### 2f. Gate 2 — mechanism category & chemist sign-off (tracking only)

A new card at the top of the Orchestration tab lets you record, per target,
which mechanism category applies (covalent reactive cysteine, metal redox
center, protein-interface disruption, structural-stabilizer local comparison,
catalytic loss-of-function) and the sign-offs that category requires. **This
does not block anything above it** — Gate 3 (DMRG/SHCI) answers a different
question (is this classically tractable) that routine screening depends on
and this does not gate. Gate 2 only becomes load-bearing once a mechanism-
specific energy calculator (e.g. a covalent ΔG‡ pipeline) exists and checks
it — that tool is not built yet, so today this is documentation, not a lock.

One-time migration (Supabase SQL editor):

```sql
create table if not exists public.gate2_records (
  target text primary key,
  category text,
  metal_present boolean,
  spin_assigned_by text,
  oxidation_state_assigned_by text,
  protonation_assigned_by text,
  scope_signoff_by text,
  ts_search_configured boolean,
  source text,
  updated_by text,
  updated_at timestamptz
);
```

No agent involvement — this is a plain save/read against the backend
(`/api/gate2/record`, `/api/gate2/list`), same as any other form on the page.

### 2h. Gate 2 as a real pre-dispatch gate — Custom Model Compound

`POST /api/gate2/dispatch_custom_compound` is the one Gate 2 endpoint that
actually blocks a dispatch (not just tracks one): it upserts the target's
category + sign-offs, computes `missing_requirements`, and refuses to queue
an `hpc_dispatch` row (`job_type='dmrg'`, geometry mode) if anything is
missing. For a target that is neither a real PDB site (Gate 1) nor already
in the fixed compound library — e.g. a metal-cluster proxy compound.

Its geometry/AVAS/charge/spin/basis are also saved to a PERMANENT store,
separate from the disposable `hpc_dispatch` queue row — a "Clear Queue"
click deletes queued/running/failed rows outright, and with them the only
copy of a hand-pasted `.xyz` geometry, found live 2026-09-02 re-testing the
same compound a second time. The "↻ reload" button on each Gate 2 table row
reads from this table (`GET /api/gate2/custom_compound/{target}`), not from
the queue.

One-time migration (Supabase SQL editor):

```sql
create table if not exists public.custom_compounds (
  target text primary key,
  geometry text,
  avas text,
  basis text,
  charge int,
  spin int,
  updated_by text,
  updated_at timestamptz
);
```

### 2i. Gate 2's first real "oven" — structural_stabilizer_local_comparison

`POST /api/gate2/stabilizer/compare` is the first mechanism-category
calculator to actually consume a Gate 2 sign-off (the other four categories
still have no downstream tool — see the module docstring). Given a target,
a wild-type DMRG classification id, and a mutant DMRG classification id, it:

1. Requires a Gate 2 record for that target in the
   `structural_stabilizer_local_comparison` category with every required
   sign-off on file (`missing_requirements` empty) — refuses (409) otherwise.
2. Fetches both DMRG records and requires them to match on every field that
   defines the computational setup — `basis`, `avas`, `avas_threshold`,
   `charge`, `spin`, `ncas`, `nelecas` — refusing (409, with the specific
   mismatched fields named) if they don't. Found live 2026-09-12: two real
   DMRG records for the SAME key had silently drifted on `avas_threshold`
   (one recorded it, one predated that column) — this check is what a
   genuinely apples-to-oranges comparison looks like before it reaches
   anyone as a number.
3. Computes ΔE = E(mutant) − E(WT) from each record's own converged
   `dmrg_energies`, seals the result (LEON's generic seal — tamper-evident,
   not the full P1-P9 physics check that schema needs a JW circuit for),
   and stores it.

`GET /api/gate2/stabilizer/list` reads every saved comparison back, newest
first. **The 2026-09-12 "first real result" (TP53_C275F vs. its own
wild-type cluster, ΔE = 2.61 mHa) is RETRACTED, found 2026-09-23**: both
sides' PDB files recorded residue 275 as CYS — no Cys→Phe substitution had
actually been applied to the structure the "mutant" run used. The apples-
to-apples check this endpoint performs (same basis, AVAS criterion,
active-space size) verifies computational SETUP, not that the two sites
being compared are actually chemically different at the mutation site —
this gap let an unmutated pair through undetected. A real C275F mutant was
subsequently built (PyMOL mutagenesis, verified atom-by-atom against the
source PDB) and re-run; that exposed a second, more general problem —
raw ΔE between clusters of DIFFERENT elemental composition (any real
substitution) is dominated by atomic-composition energy (measured ≈+166 Ha
here, matching a simple atomic-Hartree-Fock estimate almost exactly) and
carries no local-stability signal above that noise floor. **This
calculator has no valid output for a real point mutation as currently
implemented.** See the isodesmic-correction note and the Y220C cross-check
below for the (unvalidated) direction being explored instead.

One-time migration (Supabase SQL editor):

```sql
create table if not exists public.stabilizer_comparisons (
  id uuid primary key,
  created_at timestamptz not null default now(),
  target text,
  wt_dmrg_id uuid,
  mutant_dmrg_id uuid,
  e_wt_ha numeric,
  e_mutant_ha numeric,
  delta_e_ha numeric,
  delta_e_mha numeric,
  shared_setup jsonb,
  requested_by text,
  note text,
  seal text
);
```

**Isodesmic-correction direction (2026-09-23, manual, not built into the
pipeline):** balancing a reference reaction so atom types/counts cancel on
both sides (native cluster + [mutant's model compound] → mutant cluster +
[native's model compound], reusing the small-molecule RHF energies already
in `all_mutations_casscf.json`) reduces the raw atomic-composition artifact
above to a chemically plausible range: for C275F, +166 Ha → −82 mHa. Cross-
checked the same evening against TP53 Y220C, which has a real published
value (ΔTm = −8.4°C, ΔΔG ≈ +4 kcal/mol destabilizing; Joerger & Fersht) —
built from the actual native (2OCJ) and mutant (2VUK, a real crystal
structure, no in-silico mutagenesis needed) clusters, RHF-only (bypassing
`--dmrg-scf` entirely, see the instability note above — RHF is unaffected
by it and is the correct level anyway, since it matches the pre-existing
reference-compound energies). Result: **destabilizing in sign, consistent
with experiment, but ~50 kcal/mol vs. the real ~4 kcal/mol — off by roughly
an order of magnitude.** Working explanation: the reference compounds are
at their own relaxed/optimized geometries while the cluster is frozen at
its crystallographic geometry, a strain mismatch the correction does not
currently account for. **Established:** the correction points the right
direction on an independent test case. **Not established:** its magnitude
is trustworthy. **Not claimed:** any number from this approach is
citable — it needs geometry-consistent reference compounds and chemist
review before it is.

### 2g. Gate 1 — structural resolvability, with an end-of-day promotion step

A card at the very top of the Orchestration tab (before Rung 1) lets you look
up or record whether a mutation has a resolvable structure to anchor an
active space to. `targets.json` is authoritative for every mutation it
already covers (its `structure_caveat` field) — this card's own database,
`gate1_checks`, is only a **working cache** for mutations `targets.json` has
never seen (an NGS report can surface any gene). It is not a second
permanent record.

**End of day**, run `scripts/laguna/promote_gate1_checks.py` — it reads
everything accumulated in `gate1_checks`, merges each verdict into
`targets.json`, empties the cache, and tells you to review the diff before
committing:

```
python3 scripts/laguna/promote_gate1_checks.py --email you@x.com --password ...
```

(or `SOLANGE_EMAIL`/`SOLANGE_PASSWORD`, same convention as `solange_hpc.py`'s
`--agent` mode). Add `--dry-run` to see what would change without writing
anything. This script never runs automatically — run it by hand, then
`git diff targets.json`, then `python scripts/laguna/verify_consistency.py`,
then commit and push as usual.

One-time migration (Supabase SQL editor):

```sql
create table if not exists public.gate1_checks (
  target text primary key,
  pdb_id text,
  chain text,
  resi int,
  resolved boolean,
  reason text,
  source text,
  checked_by text,
  updated_at timestamptz
);
```

**What this does not do yet:** the actual PDB lookup (fetch a structure, check
whether it covers a given residue) is not automated — `POST /api/gate1/check`
only records a verdict a human already worked out by hand, the same way the
five entries in the frontend's old `STRUCTURALLY_UNRESOLVED` list were
originally decided. Automating that lookup is future work.

---

## 3. Quantum runs — Rung 4 (real IBM hardware)

Real quantum time is scarce (~10 min/month on the Open plan). Order matters:
**credentials → dry-run (free) → hardware.**

### 3a. One-time — IBM credentials

You need two things from **quantum.cloud.ibm.com**:
- an **IBM Cloud API key** (API keys page → **Create** → copy the value *once*),
- your instance **CRN** (Resource list → your Quantum instance → copy CRN;
  starts `crn:v1:bluemix:...`).

Save them (stored in `~/.qiskit`, never in the repo):

```bash
python -c "from qiskit_ibm_runtime import QiskitRuntimeService; QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token='<API_KEY>', instance='<CRN>', overwrite=True)"
```

Verify (lists backends — **free**, no quantum time):

```bash
python -c "from qiskit_ibm_runtime import QiskitRuntimeService; print([b.name for b in QiskitRuntimeService().backends()])"
```

Expect: `['ibm_fez', 'ibm_marrakesh', 'ibm_kingston']`.

> ⚠ Never paste your API key/CRN into chat or commit it. Keys are shown once —
> if a value is lost or wrong, create a new API key.

### Which backend? (any of the three)

Your Open instance has **three** QPUs, all IBM **Heron r2**, 156 qubits:
**`ibm_fez`**, **`ibm_marrakesh`**, **`ibm_kingston`**. Any online one works —
`--backend` is your choice, and the examples below use `ibm_kingston` only as a
stand-in. Pick whichever is online and least busy (a machine can go into
maintenance):

```bash
# list online backends + queue depth, pick the shortest — free, no quantum time
python -c "from qiskit_ibm_runtime import QiskitRuntimeService; s=QiskitRuntimeService(); [print(b.name, 'online' if b.status().operational else 'DOWN', b.status().pending_jobs, 'queued') for b in s.backends()]"
```

Swap `--backend <name>` in any command below (or `--agent --backend <name>`) to
whichever you chose.

### 3b. One-time — the queue migration (for the QPU agent)

Run once in the **Supabase SQL editor**:

```sql
alter table public.hpc_dispatch add column if not exists job_type text not null default 'hpc';
update public.hpc_dispatch set job_type = 'hpc' where job_type is null;
```

### 3c. Always check the pipeline first — free

```bash
python scripts/laguna/solange_qpu.py --check-credentials
python scripts/laguna/solange_qpu.py --key TP53_C275F --side native --dry-run --submit
```

The dry-run runs the whole pipeline on a local simulator (creates a
`3B-QPU-dryrun` row). `Δ(measured−HF)=0.00 mHa` is expected — no hardware noise.

### 3d. Option A — queue + QPU agent (like HPC)

Start the QPU agent (spends real quantum time on each claimed job). It runs
**every** claimed job on the `--backend` you give it, so pick an online one:

```bash
python scripts/laguna/solange_qpu.py --agent --backend ibm_kingston   # or ibm_fez / ibm_marrakesh
```

Then in SOLANGE: tick a mutation → **▶ Queue for QPU** → confirm the
real-quantum-time warning. The QPU agent pulls it → row in **Rung 4**.

### 3e. Option B — run one job directly (no agent)

```bash
python scripts/laguna/solange_qpu.py --key TP53_C275F --side native --hardware --backend ibm_kingston --submit
```

On real hardware `Δ(measured−HF)` shows **actual hardware noise** (a few mHa) —
that is the point, not an error. If `ibm_kingston` is busy, swap
`--backend ibm_marrakesh` or `ibm_fez`.

### 3f. Recover a completed job without spending quantum time

If a run hung after IBM finished it:

```bash
python scripts/laguna/solange_qpu.py --retrieve <JOB_ID> --key TP53_C275F --side native --submit
```

---

## Quick reference

```bash
# every terminal
solange

# HPC (Rung 2)
bash scripts/laguna/agent_keepalive.sh start        # start agent, then queue from UI
bash scripts/laguna/make_row.sh                     # or run directly

# DMRG (Rung 3)
bash scripts/laguna/run_dmrg.sh --key <KEY> --side native --ncas 8 --nelecas 8 --bond-dims 250,500,1000,2000 --submit

# SHCI (Rung 3, second independent classifier) — standalone, own A/B/C verdict
python3 scripts/laguna/solange_shci.py --geometry site.xyz --charge <C> --spin <S> \
  --basis sto-3g --avas "<AVAS criterion>" --key <KEY> \
  --dice-scripts ~/lili2014/Dice/scripts --sweep-eps 1e-2,1e-3,5e-4,1e-4 --submit

# SHCI cross-validated against an existing DMRG record — add its own id
python3 scripts/laguna/solange_shci.py --geometry site.xyz --charge <SAME> --spin <SAME> \
  --basis <SAME> --avas "<SAME as the DMRG run>" --key <KEY> \
  --dmrg-classification-id <uuid-from-dmrg-submit> --dice-scripts ~/lili2014/Dice/scripts --submit

# QPU (Rung 4) — --backend is any of: ibm_fez | ibm_marrakesh | ibm_kingston
python scripts/laguna/solange_qpu.py --check-credentials
python scripts/laguna/solange_qpu.py --key <KEY> --side native --dry-run --submit      # free
python scripts/laguna/solange_qpu.py --agent --backend <BACKEND>                        # agent, then queue from UI
python scripts/laguna/solange_qpu.py --key <KEY> --side native --hardware --backend <BACKEND> --submit   # direct
```

---

## Troubleshooting (things we actually hit)

| Symptom | Cause | Fix |
|---------|-------|-----|
| `CondaError: Run 'conda init'` | conda hook not loaded | it's in the `solange` function now — just type `solange` |
| `command not found: solange` in a new terminal | login shell didn't source `~/.bashrc` | already patched into your login file; if it recurs, `cd ~/lili2014` works (base is active) |
| `cannot open libmkl_def.so.1` / `undefined symbol` | block2 packaging | never run a block2 script with a bare `python` — prefix it: `bash scripts/laguna/with_block2.sh python <script> …` (this is what `run_dmrg.sh` does for you) |
| agent dot red, `1 queued` | no agent running | start it: `agent_keepalive.sh start` (HPC) / `solange_qpu.py --agent` (QPU) |
| `Provided API key could not be found` | wrong/expired API key | create a **new** IBM Cloud API key, copy once, `save_account` again |
| `not a valid instance name` | hidden char in the CRN, or wrong value | re-copy the CRN with the copy button; `save_account` again |
| `REFUSING hardware without credentials` | no IBM token/account | do §3a (save_account), then retry |
| `--dmrg-scf` crashes at CAS(~30,30) or larger (`Intel MKL ERROR: Parameter N incorrect on entry to DGEMM`, `std::length_error`, `terminate called without an active exception`, `malloc(): unsorted double linked list corrupted`, plain `Segmentation fault` — a *different* one of these each attempt, on the *identical* command) | genuine, non-deterministic instability in this Laguna pyblock2/block2 build once the active space grows past roughly CAS(30,30) — not root-caused (see note below) | if the run only needs a total energy at one fixed level of theory (e.g. an isodesmic-style comparison), bypass `--dmrg-scf`/block2 entirely: a standalone plain-PySCF RHF calculation on the same geometry/charge/basis is unaffected and completes in seconds. Otherwise: no reliable workaround yet — re-test before trusting `--dmrg-scf` above ~30 active orbitals |
| two concurrent `--dmrg-scf` runs both silently corrupt each other's resume check (`[resume] ./tmp_dmrg holds a saved MPS for a DIFFERENT active space...` followed by a crash, even though neither run's own `--scratch` flag was misspelled) | both jobs defaulted to (or were pasted with a stale copy of) the same `--scratch ./tmp_dmrg` and ran at the same time | always pass a distinct `--scratch <path>` per concurrent job; confirm it actually reached the process by checking `with_block2.sh`'s own printed `running →` line, not just the command you meant to paste |

All runs are re-verified and notarized by **LEON** on ingestion, independent of
Laguna/IBM connectivity afterward.

### `--dmrg-scf` instability at large active spaces (found 2026-09-23)

While validating an isodesmic-correction approach on TP53 Y220C (native Tyr220
vs. the real 2VUK Cys220 mutant structure), `--dmrg-scf` crashed repeatedly and
non-deterministically at active spaces around CAS(54,31)–CAS(60,68), on both
native and mutant clusters. The `[validate]` CAS(6,6) self-check inside the
same script (the element-wise DMRG-vs-FCI RDM comparison — see
`dmrgscf_block2.py`'s module docstring) always passed first, ruling out the
2-RDM axis-convention wiring as the cause; the crash is inside block2's own
`driver.dmrg(...)` C++ call, on the real active space, not in this project's
adapter code.

**Ruled out, one at a time, with direct evidence — not by assumption:**
- Thread count (`--threads 2` vs the default 4, with `OMP_NUM_THREADS`/
  `MKL_NUM_THREADS` set to match) — crashed both ways.
- block2's own internal memory pool (`--stack-mem-gb` raised from the 4.0
  default to 32) — crashed both ways.
- The shell's stack `ulimit -s` (8192 kB default vs. `unlimited`) — crashed
  both ways; a run that had survived once at the small limit later crashed
  at the unlimited one, ruling out a simple stack-overflow explanation too.
- AVAS criterion including sulfur (`S 3p`) — crashed with S included in
  four separate attempts, but a same-size run (CAS(54,31)) *without* S in
  the criterion also crashed at least once — so sulfur is not the sole
  trigger either, even though it looked that way after the first few tries.

**A separate, real bug found and fixed along the way:** two concurrent
`solange_dmrg.py --dmrg-scf` processes sharing the default `--scratch
./tmp_dmrg` corrupt each other's resume-metadata check
(`.solange_dmrg_meta.json` / MPS files) — always pass a distinct `--scratch`
per concurrent job, and verify it actually reached the process via
`with_block2.sh`'s own `running →` line (a pasted stale command silently
missing the flag is what caused this the first two times it was "tested"
tonight).

**Working conclusion:** this Laguna pyblock2/block2 build has a real,
non-deterministic instability once the active space grows past roughly
CAS(30,30)-ish, independent of thread count, memory pool size, shell stack
limit, and active-space atomic composition. Not yet root-caused — candidates
include an internal block2 buffer-sizing bug at scale, or an MKL/ABI mismatch
specific to this environment's `with_block2.sh` LD_PRELOAD setup. **Do not
treat `--dmrg-scf` as reliable above ~30 active orbitals on Laguna without
re-testing.**

**Practical workaround used tonight:** for the Y220C check, only a total
energy at one fixed level of theory (RHF, matching the pre-existing
p-cresol/methanethiol reference-compound energies) was actually needed —
not a converged DMRG classification. A ~20-line standalone PySCF script
(no block2, no AVAS, no CASCI) gave both cluster RHF energies reliably in
seconds. Reach for this whenever the question is "what does a fixed-method
total energy say," not "what is this active space's true entanglement" —
the latter still needs `--dmrg-scf` (or a smaller, FCI-tractable
`--casci` space) to answer at all.

**Code-level attempt, same night, partial result — `dmrgscf_block2.py`'s
`Block2FCISolver.kernel()`:** the cold-solve ramp was originally handed to
block2 as ONE `driver.dmrg()` call carrying the whole 50→maxM bond-dimension
list internally. `run_dmrg()` (the downstream ladder, never observed to
crash all night on the same active spaces) instead calls `driver.dmrg()`
once per bond dimension, each with a single-element list. Splitting the
cold-solve ramp into one call per rung, mirroring that pattern, raised the
success rate on the previously-100%-crashing CAS(60,35)+sulfur case to 2 of
3 attempts — a real improvement, not nothing — but did **not** eliminate the
crash (a third, identical attempt crashed again). The cold-solve step itself
shows **no timing regression**: 871.6s / 852.7s across two successful native
runs, essentially identical. (An earlier version of this note claimed a
5-6x slowdown, based on watching cumulative CPU time on a still-running
process and wrongly assuming it was stuck at the cold-solve step it had
already passed, when it had actually moved on to the downstream
bond-dimension ladder — corrected once the actual per-step log timings were
read directly, rather than inferred from an in-progress `ps` snapshot.)
**Status: improves the odds, does not fix the bug. Do not report this as
resolved.**

**If the instability needs to actually go away, not just improve:**
`CheMPS2` is a more realistic path than continuing to patch this adapter —
unlike pyblock2, it has an **official, community-maintained** PySCF CASSCF
integration (`pyscf.dmrgscf`) built and tested by others for exactly this
"DMRG as CASSCF's inner solver" role, at active-space sizes matching where
this build of block2 struggles (20-40 orbitals). Switching would retire
`dmrgscf_block2.py`'s from-scratch 2-RDM convention validation entirely
(`validate()`/`diagnose()`), since that convention is already the
community's problem to have gotten right, not this project's own. Real
cost: installing/compiling CheMPS2 on Laguna (untested, may hit its own
linking issues) and a short benchmark run — smaller than the validation
work already sunk into block2, but not zero. Not attempted tonight.
