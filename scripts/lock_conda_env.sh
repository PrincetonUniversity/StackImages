#!/usr/bin/env zsh
# Regenerate ci/conda-linux-64.lock from environment.yml inside a pinned
# Linux container, so the CI desktop environment is fully reproducible.
# Usage: scripts/lock_conda_env.sh [--dry-run]
set -euo pipefail

IMAGE="continuumio/miniconda3:25.3.1-1"
ROOT="${0:A:h:h}"
DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    *) print -u2 "Unknown argument: $arg"; exit 2 ;;
  esac
done

cmd=(docker run --rm --platform linux/amd64
  -v "$ROOT/environment.yml:/env.yml:ro" -v "$ROOT/ci:/out"
  "$IMAGE" bash -euo pipefail -c
  'conda env create -q -n lock -f /env.yml >/dev/null
   conda list -n lock --explicit --md5 > /out/conda-linux-64.lock.tmp
   mv /out/conda-linux-64.lock.tmp /out/conda-linux-64.lock')

if (( DRY_RUN )); then
  print -- "Would run: ${(q)cmd[@]}"
  exit 0
fi
"${cmd[@]}"
print "Wrote $ROOT/ci/conda-linux-64.lock"
