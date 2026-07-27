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

One word brings it up — the `solange` shortcut enters the repo, starts the
agent, and prints its status:

```bash
solange
```

SOLANGE's **QPU agent** dot turns green within ~15s. The agent runs detached
(survives closing the terminal). Look for `RUNNING` and
`IBM Quantum credentials: AVAILABLE` in the status.

**One-time setup** (only if `solange` isn't defined yet on this Mac — e.g. a
fresh machine). Paste this once, then `source ~/.bash_profile`:

```bash
cat >> ~/.bash_profile <<'EOF'

# SOLANGE — bring the QPU agent up in one word
solange() {
  cd ~/Desktop/lili2014-qpu && \
  bash scripts/laguna/qpu_keepalive.sh start && \
  bash scripts/laguna/qpu_keepalive.sh status
}
EOF
```

To pull the latest code first (e.g. new Hamiltonians), before `solange`:

```bash
cd ~/Desktop/lili2014-qpu && git pull origin main
```

Default backend is `ibm_kingston`. To use a different one, or to check/stop:

```bash
SOLANGE_QPU_BACKEND=ibm_fez bash scripts/laguna/qpu_keepalive.sh start
bash scripts/laguna/qpu_keepalive.sh status
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
