# SOLANGE — HPC runner for the Laguna cluster

Scales the Phase-3A active space beyond the 4-qubit laptop proxy, on a real
SLURM + GPU cluster, and emits artifacts in the **exact schema SOLANGE already
consumes** — so plugging the results back in is a merge, not a rewrite.

## What it produces (the two "sockets")

| File | Socket in SOLANGE | What it is |
|------|-------------------|------------|
| `out/jw_<KEY>_<side>.json` | `jw_hamiltonians.json` | the scaled Hamiltonian (Pauli terms + energies) |
| `out/provenance_<KEY>_<side>.json` | Compliance tab / Supabase | a **signed P1–P9 record** (P8 = SHA-256, re-verifiable with no GPU) |

The P8 seal uses the *same* algorithm as the SOLANGE backend (`build_p8_seal`),
so the platform certifies the record's integrity without re-running the physics.
`p3_backend` and `provenance_source` state plainly that the run happened on
Laguna — an HPC record is never dressed up as a backend run.

## The honest scaling ladder

A minimal STO-3G model compound has too few orbitals for a large active space.
Growing the active space means growing **two axes together**:

| Tier | basis | CAS | qubits | runs where |
|------|-------|-----|--------|-----------|
| dry-run (validate) | `sto-3g` | (2,2) | 4 | anywhere — must match the live 4q result |
| intermediate | `6-31g` | (8,8) | 16 | any GPU |
| benchmark edge | `cc-pvdz` | (12,12)…(16,16) | 24…32 | 80 GB GPU |
| full site (CDKN2A/ARID2) | `cc-pvdz` | up to full | 40…56 | `largemem` + DMRG |

`--ncas/--nelecas/--basis` are explicit flags. Nothing scales silently.

## Steps on Laguna

1. **Clone the repo** (so `generate_expansion_jw.py` — the geometry table — is importable):
   ```bash
   git clone <repo-url> ~/lili2014 && cd ~/lili2014
   ```
2. **Environment** (one time): PySCF + openfermion for the Hamiltonian; PennyLane
   + lightning.gpu only if you pass `--vqe`:
   ```bash
   python -m pip install --user pyscf openfermion scipy numpy
   python -m pip install --user pennylane pennylane-lightning[gpu]
   ```
3. **Preflight** — confirm the GPU and the qubit ceiling the script will allow:
   ```bash
   srun -p gpu --gres=gpu:1 -t 00:02:00 \
       python scripts/laguna/solange_hpc.py --compound acetamide --basis sto-3g \
       --ncas 2 --nelecas 2 --key ARID2_LOF --side native --residue dry-run --out ./out
   ```
   The banner prints `statevector ceiling ≈ N qubits`. The (2,2) energy MUST match
   the value already in `jw_hamiltonians.json` for `ARID2_LOF/native` — that proves
   the toolchain before you trust any scaled number.
4. **Submit the scaled run**:
   ```bash
   sbatch --export=ALL,COMPOUND=acetamide,BASIS=6-31g,NCAS=8,NELECAS=8,KEY=ARID2_LOF,SIDE=native,RESIDUE="Gln1118 amide" \
       scripts/laguna/submit_gpu.sbatch
   ```
5. **Send both JSON files back** (`out/jw_*.json`, `out/provenance_*.json`).
   They get synced into `backend/jw_hamiltonians.json` (+ root) and the Compliance
   view, then pushed — exactly the ARID2 loop.

## Notes

- `nvidia-smi` drives the qubit ceiling automatically; the job refuses a CAS that
  won't fit VRAM rather than OOM-ing mid-run.
- Exact diagonalisation is used as the reference up to 16 qubits; beyond that,
  trust VQE / CASSCF (the script says so in its output).
- For the **full active sites** of CDKN2A (40q) and ARID2 (56q), statevector won't
  fit — those go through DMRG on `largemem`/`oneweek` (separate script, once the
  GPU ladder is validated).
- **De-identified data only.** These runs use model-compound geometries, never PHI.
