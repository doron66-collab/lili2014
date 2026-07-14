#!/usr/bin/env bash
# One-command launcher for the SOLANGE Copilot.
# Frees the port if a previous server is still running, then starts the
# zero-dependency server. Open http://localhost:8000/app in your browser.
set -e
PORT="${PORT:-8000}"
DIR="$(cd "$(dirname "$0")" && pwd)"

# stop any old server holding the port (ignore if none)
lsof -ti:"$PORT" 2>/dev/null | xargs kill -9 2>/dev/null || true

echo "──────────────────────────────────────────────"
echo "  Starting SOLANGE Copilot on port $PORT"
echo "  Open:  http://localhost:$PORT/app"
echo "  Stop:  press Ctrl+C"
echo "──────────────────────────────────────────────"
cd "$DIR/backend"
exec python3 serve.py
