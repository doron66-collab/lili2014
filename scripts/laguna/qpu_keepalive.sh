#!/usr/bin/env bash
# qpu_keepalive.sh — keep the SOLANGE QPU agent alive independently of the
# terminal (Mac/local machine version of agent_keepalive.sh). Start it once and
# forget it: it runs detached (nohup), and a tiny supervisor restarts the agent
# if it crashes. Closing the terminal or a transient crash no longer takes the
# agent down.
#
#   bash scripts/laguna/qpu_keepalive.sh start    # launch detached + auto-restart
#   bash scripts/laguna/qpu_keepalive.sh status    # is it running? recent log?
#   bash scripts/laguna/qpu_keepalive.sh stop      # stop the supervisor + agent
#
# Each claimed job spends REAL quantum time on SOLANGE_QPU_BACKEND — this agent
# only claims jobs YOU queued in SOLANGE (job_type=qpu); it never runs on its own.
#
# HONEST LIMIT (same posture as agent_keepalive.sh): nothing here keeps the agent
# alive past your machine being on/awake and your IBM Quantum credentials being
# valid. If the Mac sleeps, or the saved IBM account expires/is revoked, the agent
# stops — re-run `start` after waking the machine / refreshing credentials.
set -uo pipefail

# ── EDIT these to your environment (or export them before calling) ───────────
VENV_DIR="${SOLANGE_QPU_VENV:-$HOME/solange-venv}"           # python venv with qiskit-ibm-runtime
REPO_DIR="${SOLANGE_QPU_REPO:-$HOME/Desktop/lili2014-qpu}"    # where you cloned the repo
BACKEND="${SOLANGE_QPU_BACKEND:-ibm_kingston}"                # ibm_kingston | ibm_fez | ibm_marrakesh
# ─────────────────────────────────────────────────────────────────────────────

STATE_DIR="${SOLANGE_STATE:-$HOME/.solange}"
PID_FILE="$STATE_DIR/qpu.supervisor.pid"
LOG_FILE="$STATE_DIR/qpu.log"
mkdir -p "$STATE_DIR" "$REPO_DIR/out" 2>/dev/null || true

_running() { [[ -f "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE" 2>/dev/null)" 2>/dev/null; }

# The supervisor body — runs detached. Re-launches the agent whenever it exits,
# with a short backoff so a hard failure (e.g. expired credentials) doesn't hot-loop.
_supervisor() {
  if [[ -f "$VENV_DIR/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "$VENV_DIR/bin/activate"
  else
    echo "[keepalive] venv not found at $VENV_DIR — set SOLANGE_QPU_VENV." >&2
    exit 1
  fi
  cd "$REPO_DIR" || { echo "[keepalive] repo not found: $REPO_DIR — set SOLANGE_QPU_REPO."; exit 1; }
  while true; do
    echo "[keepalive $(date '+%F %T')] starting QPU agent (backend=$BACKEND)"
    # -u = unbuffered, so the agent's login/poll/job lines reach the log file live.
    python -u scripts/laguna/solange_qpu.py --agent --backend "$BACKEND" --out ./out
    echo "[keepalive $(date '+%F %T')] QPU agent exited (rc=$?); restarting in 10s"
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
    echo "QPU agent supervisor started (PID $(cat "$PID_FILE"))."
    echo "  backend: $BACKEND"
    echo "  logs:    tail -f $LOG_FILE"
    echo "  status:  bash $0 status"
    echo "SOLANGE's QPU agent dot should turn green within ~15s."
    echo "You can now close this terminal — the agent keeps running."
    ;;
  --supervise)   # internal entry point (called by nohup above)
    _supervisor
    ;;
  status)
    if _running; then
      echo "RUNNING — supervisor PID $(cat "$PID_FILE") · backend=$BACKEND"
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
