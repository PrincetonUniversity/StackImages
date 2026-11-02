"""Decide the next desktop release tag, or that no release is needed.

Scheme: v<ass.py version>.<n>, e.g. v21.0, v21.1, …; a new ass.py version (Bob bumps
`version = N`) starts again at .0. A release is only made when the bundle's app and data
files changed since the newest versioned release (compared by fingerprint, which
scripts/package_pywin.sh --fingerprint prints and each release's notes record).
Older v-<sha> and snapshots-* releases are ignored.

Usage (CI):
  gh api repos/OWNER/REPO/releases --paginate > releases.json
  python3 scripts/next_release.py --app-version 21 --fingerprint ABC --releases releases.json
Prints the tag to create, or nothing when the bundle is unchanged.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

TAG = re.compile(r"^v(\d+)\.(\d+)$")
FINGERPRINT_LINE = "Bundle fingerprint: "


def next_tag(releases: list[dict], app_version: int, fingerprint: str) -> str | None:
    versioned = []
    for release in releases:
        match = TAG.match(str(release.get("tag_name", "")))
        if match and not release.get("draft"):
            versioned.append(((int(match.group(1)), int(match.group(2))), release))
    if versioned:
        _, newest = max(versioned, key=lambda item: item[0])
        if f"{FINGERPRINT_LINE}{fingerprint}" in (newest.get("body") or ""):
            return None
    minors = [minor for (major, minor), _ in versioned if major == app_version]
    return f"v{app_version}.{max(minors) + 1 if minors else 0}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--app-version", type=int, required=True)
    parser.add_argument("--fingerprint", required=True)
    parser.add_argument("--releases", type=Path, required=True, help="JSON from the GitHub releases API")
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[0-9a-f]{40,64}", args.fingerprint):
        raise SystemExit(f"Not a fingerprint: {args.fingerprint!r}")
    pages = json.loads(args.releases.read_text() or "[]")
    # `gh api --paginate` concatenates pages; accept a list or a list of lists.
    releases = [r for page in pages for r in (page if isinstance(page, list) else [page])]
    tag = next_tag(releases, args.app_version, args.fingerprint)
    print(tag or "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
