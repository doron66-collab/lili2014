# SOLANGE — Startup Guide

Bring the whole system up: SOLANGE in the browser, the classical (HPC + DMRG)
agent on Laguna, and the QPU agent on your Mac. Order doesn't matter — just
have all three running before you queue anything.

**0.** Open SOLANGE in your browser (`solange-platform.bio`).

---

## Chapter 1 — Laguna: HPC + DMRG agent

One agent handles **both** HPC (CASSCF/VQE) and DMRG classifications — you
don't start two.

**Where:** CARC OnDemand → JupyterLab (compute node) → Terminal.

```bash
solange
bash scripts/laguna/agent_keepalive.sh start
```

That's it — the agent runs detached (survives closing the terminal). SOLANGE's
**HPC/DMRG agent** dot turns green within ~15s.

```bash
# check on it any time
bash scripts/laguna/agent_keepalive.sh status

# stop it
bash scripts/laguna/agent_keepalive.sh stop
```

---

## Chapter 2 — Mac: QPU agent

Runs on the Mac (not Laguna) because it needs reliable internet to IBM Cloud.
Spends **real quantum time** on every job it claims — only claims jobs you
queued in SOLANGE.

**Where:** Terminal.app on your Mac.

```bash
cd ~/Desktop/lili2014-qpu
git pull origin claude/code-access-clarification-ab1W8
bash scripts/laguna/qpu_keepalive.sh start
```

Runs detached (survives closing the terminal) — same as Chapter 1. SOLANGE's
**QPU agent** dot turns green within ~15s. Default backend is `ibm_kingston`;
to use a different one:

```bash
SOLANGE_QPU_BACKEND=ibm_fez bash scripts/laguna/qpu_keepalive.sh start
```

```bash
# check on it any time
bash scripts/laguna/qpu_keepalive.sh status

# stop it
bash scripts/laguna/qpu_keepalive.sh stop
```

---

## Then — queue and watch

In SOLANGE (**Orchestration** tab): load an NGS report if the mutation list
isn't showing, tick a mutation, and:

| Button | Rung | Agent |
|---|---|---|
| ▶ Dispatch selected to HPC | 2 | Laguna |
| ▶ Queue DMRG | 3 | Laguna |
| ▶ Queue for QPU | 4 | Mac |

Results land in their rung automatically — no manual refresh needed while a
job is in flight (the ladder polls on its own).
