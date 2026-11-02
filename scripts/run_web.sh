#!/usr/bin/env zsh
# Build the Show Clusters site and serve it locally (http://127.0.0.1:PORT/).
# Uses the web-tools container, so nothing is installed on this machine.
#
# Usage: scripts/run_web.sh [--port 8000] [--snapshots DIR] [--dry-run]
set -euo pipefail

ROOT="${0:A:h:h}"
PORT=8000
SNAPSHOTS=""
DRY_RUN=0
while (( $# )); do
  case "$1" in
    --port) PORT="${2:?--port needs a value}"; shift ;;
    --snapshots) SNAPSHOTS="${2:?--snapshots needs a directory}"; shift ;;
    --dry-run) DRY_RUN=1 ;;
    -h|--help) sed -n '2,6p' "$0"; exit 0 ;;
    *) print -u2 "Unknown argument: $1"; exit 2 ;;
  esac
  shift
done
[[ "$PORT" == <1024-65535> ]] || { print -u2 "Port must be 1024-65535"; exit 2; }
site_args=(--no-snapshots)
[[ -n "$SNAPSHOTS" ]] && site_args=(--snapshots "$SNAPSHOTS")

vendor=(docker-compose run --rm -T web-tools python scripts/vendor_pyodide.py)
build=(docker-compose run --rm -T web-tools python scripts/build_site.py "${site_args[@]}")
serve=(docker-compose run --rm -T -p "127.0.0.1:$PORT:8000" web-tools
       python -m http.server 8000 --bind 0.0.0.0 -d _site)
cd "$ROOT"
if (( DRY_RUN )); then
  for cmd in vendor build serve; do print -r -- "Would run: ${(q)${(P)cmd}[@]}"; done
  exit 0
fi
"${vendor[@]}"
"${build[@]}"
print "Serving http://127.0.0.1:$PORT/ShowClusters/ (Ctrl-C to stop)"
"${serve[@]}"
