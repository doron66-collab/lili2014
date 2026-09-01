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

### 2d. Queue SHCI from the browser — "▶ Queue SHCI" button (Rung 3)

A DMRG record with stored geometry/AVAS (any real `--geometry` run submitted
after this was added) shows a **▶ Queue SHCI** button on its own row in the
Rung 3 table. Clicking it copies that record's own geometry/AVAS/charge/spin —
nothing re-typed — and queues an SHCI job for the SAME classical agent that
already runs HPC/DMRG jobs (`solange_hpc.py --agent`, started once via
`agent_keepalive.sh start`). No terminal interaction needed per job once that
agent is running.

One-time migrations (Supabase SQL editor) — additive, safe to re-run:

```sql
-- dmrg_classifications: store the reproducibility inputs a later SHCI
-- cross-validation needs to copy (only populated for a real --geometry run,
-- not --compound demo mode).
alter table public.dmrg_classifications add column if not exists geometry text;
alter table public.dmrg_classifications add column if not exists avas text;
alter table public.dmrg_classifications add column if not exists charge int;
alter table public.dmrg_classifications add column if not exists spin int;

-- hpc_dispatch: carries an SHCI job's parameters from the "Queue SHCI" button
-- to the agent that picks it up.
alter table public.hpc_dispatch add column if not exists geometry text;
alter table public.hpc_dispatch add column if not exists avas text;
alter table public.hpc_dispatch add column if not exists charge int;
alter table public.hpc_dispatch add column if not exists spin int;
alter table public.hpc_dispatch add column if not exists sweep_eps text;
alter table public.hpc_dispatch add column if not exists dmrg_classification_id uuid;
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

All runs are re-verified and notarized by **LEON** on ingestion, independent of
Laguna/IBM connectivity afterward.
