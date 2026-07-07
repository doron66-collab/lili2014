# SOLANGE on USC Discovery (CARC) — step-by-step for first-timers

A beginner-friendly walkthrough to run the SOLANGE HPC pipeline on the
**USC Discovery** cluster via an **SSH terminal**. Each step says *why*, not just
*what*. Work through it top to bottom.

Cluster facts (for reference):
- Scheduler: **SLURM**. GPU partition: `gpu` (A100 80GB/40GB, L40S, A40, V100, P100).
- We target **A100-80GB** — best FP64, ~30-qubit statevector ceiling.
- Home directory: `/home1/<NetID>`.

---

## Step 0 — Have ready
- **USC NetID + password**, and your **phone** for Duo 2FA.
- A terminal: built in on macOS/Linux; on Windows use PowerShell or Windows Terminal.

---

## Step 1 — Log in (SSH)
```bash
ssh <your-NetID>@discovery.usc.edu
```
Enter your password, then approve the Duo prompt on your phone.

You are now on a **login node**. RULE: never run heavy compute here — it is shared.
We send the real work to a compute node in Step 5.

---

## Step 2 — Get oriented
```bash
whoami      # your username
pwd         # should be /home1/<NetID> — your home directory
```

---

## Step 3 — Get the code
```bash
cd ~
git clone <REPO-URL> lili2014      # first time
cd lili2014
# next time instead:  cd ~/lili2014 && git pull
```
This brings down `scripts/laguna/` with the runner + this guide.

---

## Step 4 — Build the Python environment (one time only)
```bash
module purge
module load conda
mamba create -n solange python=3.11 -y
mamba activate solange
pip install pyscf openfermion scipy numpy
pip install pennylane pennylane-lightning[gpu]
```
Takes a few minutes. Next time you only need:
`module load conda && mamba activate solange`.

If `module load conda` fails, run `module avail conda` and use the name it lists.

---

## Step 5 — Grab an interactive GPU node
This moves you off the login node onto a real **A100**:
```bash
salloc --partition=gpu --gres=gpu:a100:1 --cpus-per-task=8 --mem=32G --time=01:00:00
```
Wait for the allocation; the prompt changes to a compute-node name. Confirm the GPU:
```bash
nvidia-smi --query-gpu=name,memory.total --format=csv
# expect e.g.  NVIDIA A100-SXM4-80GB, 81920 MiB
```

---

## Step 6 — The dry-run (sanity check)
```bash
cd ~/lili2014
mamba activate solange        # if not already active
python scripts/laguna/solange_hpc.py --compound acetamide --basis sto-3g \
  --ncas 2 --nelecas 2 --key ARID2_LOF --side native --residue dry-run --out ./out
```

---

## Step 7 — What to check
- Banner: `GPU: NVIDIA A100-80GB ... statevector ceiling ≈ 30 qubits`
  → auto-detection works.
- `CASSCF E = -205.31950...` → **must** be ≈ -205.3195. If yes, the toolchain is
  correct and consistent with the live SOLANGE value. 🎉
- Ends with `WROTE .../jw_ARID2_LOF_native.json` and `provenance_...json`.

The script has built-in guards: it enforces CASSCF convergence and refuses to emit
a Hamiltonian whose active-space FCI disagrees with e_casscf — so a bad run errors
out loudly instead of producing a wrong number.

---

## Step 8 — Send the result back
Copy the Step-6 output (or just: "worked, CASSCF matches"). If anything failed,
send the full error text. Then we move to the first ladder run: **CAS(8,8)** on GPU:
```bash
python scripts/laguna/solange_hpc.py --compound acetamide --basis 6-31g \
  --ncas 8 --nelecas 8 --key ARID2_LOF --side native --residue "Gln1118 amide" \
  --vqe --out ./out
```

---

## Two safety rules
- ⚠️ Never run `python solange_hpc.py` on the login node — only after `salloc`.
- ⏱️ `salloc` gives you 1 hour. If it expires, just run `salloc` again.

## If you get stuck
Send me: the exact command you ran + the full output/error. Common first-time snags:
- `module load conda` name differs → `module avail conda`.
- No A100 free → try `--gres=gpu:a40:1` (48GB, ~29-qubit ceiling; the script
  auto-caps to fit) or wait and retry.
- `salloc` pending a long time → the `gpu` partition is busy; retry later or lower
  `--mem`/`--time`.
