#!/usr/bin/env bash
# agent_keepalive.sh — keep the SOLANGE Laguna agent alive independently of the
# terminal. Start it once and forget it: it runs detached (nohup), and a tiny
# supervisor restarts the agent if it crashes. Closing the terminal, an SSH drop,
# or a transient crash no longer takes the agent down.
#
#   bash scripts/laguna/agent_keepalive.sh start    # launch detached + auto-restart
#   bash scripts/laguna/agent_keepalive.sh status    # is it running? recent log?
#   bash scripts/laguna/agent_keepalive.sh stop      # stop the supervisor + agent
#
# HONEST LIMIT (this is the security boundary, by design — §06.iii / DP2):
# nothing here — and nothing SOLANGE could ever do from the browser — keeps the
# agent alive past YOUR authenticated Laguna session. When Duo 2FA expires, the
# SLURM allocation ends, or the node reboots, the agent stops. That ceiling is the
# guarantee, not a bug: no credential lives outside the cluster and no inbound
# channel exists to restart it remotely. Re-authenticate and run `start` again.
set -uo pipefail

# ── EDIT these two to your environment (or export them before calling) ───────
CONDA_ENV="${SOLANGE_ENV:-solange}"           # conda env with pyscf + pennylane + numpy
REPO_DIR="${SOLANGE_REPO:-$HOME/lili2014}"    # where you cloned the repo on Laguna
# ─────────────────────────────────────────────────────────────────────────────

STATE_DIR="${SOLANGE_STATE:-$HOME/.solange}"
PID_FILE="$STATE_DIR/agent.supervisor.pid"
LOG_FILE="$STATE_DIR/agent.log"
mkdir -p "$STATE_DIR" "$REPO_DIR/out" 2>/dev/null || true

_running() { [[ -f "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE" 2>/dev/null)" 2>/dev/null; }

# The supervisor body — runs detached. Re-launches the agent whenever it exits,
# with a short backoff so a hard failure (e.g. login expired) doesn't hot-loop.
_supervisor() {
  module load conda 2>/dev/null || true
  # shellcheck disable=SC1091
  source activate "$CONDA_ENV" 2>/dev/null || conda activate "$CONDA_ENV" 2>/dev/null || true
  cd "$REPO_DIR" || { echo "[keepalive] repo not found: $REPO_DIR"; exit 1; }
  while true; do
    echo "[keepalive $(date '+%F %T')] starting agent (env=$CONDA_ENV)"
    python scripts/laguna/solange_hpc.py --agent --out ./out
    echo "[keepalive $(date '+%F %T')] agent exited (rc=$?); restarting in 10s"
    sleep 10
  done
}

case "${1:-start}" in
  start)
    if _running; then
      echo "Already running (supervisor PID $(cat "$PID_FILE")). Logs: $LOG_FILE"
      exit 0
    fi
    # Re-exec ourselves in --supervise mode, detached from this terminal.
    nohup bash "$0" --supervise >>"$LOG_FILE" 2>&1 &
    echo $! >"$PID_FILE"
    disown 2>/dev/null || true
    echo "Agent supervisor started (PID $(cat "$PID_FILE"))."
    echo "  logs:   tail -f $LOG_FILE"
    echo "  status: bash $0 status"
    echo "SOLANGE's agent dot should turn green within ~15s."
    ;;
  --supervise)   # internal entry point (called by nohup above)
    _supervisor
    ;;
  status)
    if _running; then
      echo "RUNNING — supervisor PID $(cat "$PID_FILE")"
    else
      echo "NOT running."
    fi
    if [[ -f "$LOG_FILE" ]]; then echo "── last log lines ──"; tail -n 8 "$LOG_FILE"; fi
    ;;
  stop)
    if _running; then
      pid="$(cat "$PID_FILE")"
      # kill the supervisor and any agent/python children it spawned
      pkill -P "$pid" 2>/dev/null || true
      kill "$pid" 2>/dev/null || true
      rm -f "$PID_FILE"
      echo "Stopped."
    else
      echo "Not running."; rm -f "$PID_FILE"
    fi
    ;;
  *)
    echo "usage: bash $0 {start|status|stop}"; exit 2;;
esac
