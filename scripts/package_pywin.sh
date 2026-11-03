#!/usr/bin/env zsh
# Build PyWin.tar: the minimal desktop StackImages bundle students download.
#
# Only committed files on the allowlist are packaged (via `git archive`), under a
# top-level PyWin/ folder, matching the README's
# `tar -xvf StackImages.tar --strip-components=1 PyWin`.
# Web, CI, test and tooling files are never included.
#
# Usage: scripts/package_pywin.sh [--dry-run] [--ref <git-ref>] [--out <path>]
#        scripts/package_pywin.sh --fingerprint | --app-version [--ref <git-ref>]
#   --fingerprint  hash of the app and data files in the bundle (not README/LICENSE);
#                  unchanged fingerprint = nothing new to release
#   --app-version  the `version = N` number in ass.py (the release's major version)
set -euo pipefail

ROOT="${0:A:h:h}"
REF="HEAD"
OUT="$ROOT/dist/PyWin.tar"
DRY_RUN=0
MODE=build

while (( $# )); do
  case "$1" in
    --dry-run) DRY_RUN=1 ;;
    --fingerprint) MODE=fingerprint ;;
    --app-version) MODE=app-version ;;
    --ref) REF="${2:?--ref needs a value}"; shift ;;
    --out) OUT="${2:?--out needs a value}"; shift ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) print -u2 "Unknown argument: $1"; exit 2 ;;
  esac
  shift
done

# What the desktop app needs at runtime, and nothing else.
INCLUDE=(
  ass.py
  GetFit.py
  ShowClusters.py
  StackImages.py
  bias_file.txt
  flat_file.txt
  csv
  fit
  png
  txt
  LICENSE
  README.md
)
# Must never appear in the bundle (checked after building).
FORBIDDEN_PREFIXES=(web/ tests/ ci/ scripts/ docs/ .github/ node_modules/ dist/ .env)

cd "$ROOT"
git rev-parse --verify --quiet "$REF^{commit}" >/dev/null \
  || { print -u2 "Not a commit: $REF"; exit 1; }

for item in "${INCLUDE[@]}"; do
  if [[ -z "$(git ls-tree --name-only "$REF" -- "$item")" ]]; then
    print -u2 "Missing from $REF: $item"
    exit 1
  fi
done

if [[ $MODE == fingerprint ]]; then
  # Tree entries carry each file's content hash, so this changes exactly when a file does.
  git ls-tree -r "$REF" -- "${INCLUDE[@]:#(README.md|LICENSE)}" | git hash-object --stdin
  exit 0
fi
if [[ $MODE == app-version ]]; then
  version=$(git show "${REF}:ass.py" | tr -d '\r' | sed -nE 's/^version *= *([0-9]+) *$/\1/p' | head -n1)
  [[ -n "$version" ]] || { print -u2 "No 'version = N' line in ass.py at $REF"; exit 1; }
  print -- "$version"
  exit 0
fi

if (( DRY_RUN )); then
  print "Would write $OUT from $REF with:"
  git ls-tree -r --name-only "$REF" -- "${INCLUDE[@]}" | sed 's#^#  PyWin/#'
  exit 0
fi

mkdir -p "${OUT:h}"
tmp="$OUT.partial"
trap 'rm -f "$tmp"' EXIT
git archive --format=tar --prefix=PyWin/ -o "$tmp" "$REF" -- "${INCLUDE[@]}"

# Verify before publishing the file: every member under PyWin/, nothing forbidden,
# and each entry point present.
members=("${(@f)$(tar -tf "$tmp")}")
(( ${#members} > 0 )) || { print -u2 "Empty archive"; exit 1; }
for m in "${members[@]}"; do
  [[ "$m" == PyWin/* ]] || { print -u2 "Outside PyWin/: $m"; exit 1; }
  rel="${m#PyWin/}"
  for bad in "${FORBIDDEN_PREFIXES[@]}"; do
    [[ "$rel" == "$bad"* ]] && { print -u2 "Forbidden file packaged: $m"; exit 1; }
  done
done
for required in ass.py StackImages.py ShowClusters.py GetFit.py; do
  (( ${members[(Ie)PyWin/$required]} )) || { print -u2 "Missing entry point: $required"; exit 1; }
done

mv "$tmp" "$OUT"
trap - EXIT
print "Wrote $OUT (${#members} entries from $(git rev-parse --short "$REF"))"
