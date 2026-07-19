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

Start the QPU agent (spends real quantum time on each claimed job):

```bash
python scripts/laguna/solange_qpu.py --agent --backend ibm_kingston
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

# QPU (Rung 4)
python scripts/laguna/solange_qpu.py --check-credentials
python scripts/laguna/solange_qpu.py --key <KEY> --side native --dry-run --submit      # free
python scripts/laguna/solange_qpu.py --agent --backend ibm_kingston                    # agent, then queue from UI
python scripts/laguna/solange_qpu.py --key <KEY> --side native --hardware --backend ibm_kingston --submit   # direct
```

---

## Troubleshooting (things we actually hit)

| Symptom | Cause | Fix |
|---------|-------|-----|
| `CondaError: Run 'conda init'` | conda hook not loaded | it's in the `solange` function now — just type `solange` |
| `command not found: solange` in a new terminal | login shell didn't source `~/.bashrc` | already patched into your login file; if it recurs, `cd ~/lili2014` works (base is active) |
| `cannot open libmkl_def.so.1` / `undefined symbol` | block2 packaging | use `run_dmrg.sh` — it builds the LD_PRELOAD automatically |
| agent dot red, `1 queued` | no agent running | start it: `agent_keepalive.sh start` (HPC) / `solange_qpu.py --agent` (QPU) |
| `Provided API key could not be found` | wrong/expired API key | create a **new** IBM Cloud API key, copy once, `save_account` again |
| `not a valid instance name` | hidden char in the CRN, or wrong value | re-copy the CRN with the copy button; `save_account` again |
| `REFUSING hardware without credentials` | no IBM token/account | do §3a (save_account), then retry |

All runs are re-verified and notarized by **LEON** on ingestion, independent of
Laguna/IBM connectivity afterward.
