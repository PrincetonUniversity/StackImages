"""The web site's version: web-MAJOR.N releases, their changelog entries, and the page label.

Scheme: web-1.0, web-1.1, … One version covers the whole site (home page, Show Clusters,
HR Diagram of Your Own Image). Each deploy that changes the site is the next N; a deploy of
an unchanged site keeps its version. "Changed" is decided by the site fingerprint
(scripts/zip_site.py --fingerprint), which each release's notes record and which leaves out
the build time and version.js, so stamping a version is not itself a change. MAJOR comes from
web/VERSION; raising it starts again at .0. Each version's notes are its section of
docs/CHANGELOG.md ("## 1.4 (2026-11-07)" and the lines below it).

Usage (CI, pages.yml):
  gh api "repos/OWNER/REPO/releases?per_page=100" --paginate --slurp > build/releases.json
  python scripts/web_version.py next --fingerprint FP --releases build/releases.json
      prints "1.4 new 2026-11-07" (a new version) or "1.3 same 2026-11-05" (unchanged site)
  python scripts/web_version.py notes --version 1.4      the changelog entry; exit 1 if none
  python scripts/web_version.py stamp --site _site --version 1.4 --date 2026-11-07
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAJOR_FILE = ROOT / "web" / "VERSION"
CHANGELOG = ROOT / "docs" / "CHANGELOG.md"
PAGE_FOLDERS = ("", "ShowClusters", "HRDiagram")  # each page folder has its own version.js
PLACEHOLDER = 'const SITE_VERSION = { version: "", date: "" };'
FINGERPRINT_LINE = "Site fingerprint: "
TAG = re.compile(r"^web-(\d+)\.(\d+)$")
NUMBER = re.compile(r"^\d+\.\d+$")
DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


@dataclass(frozen=True)
class Version:
    number: str  # "1.4"
    date: str  # release date, YYYY-MM-DD
    new: bool  # True: publish a release for it

    @property
    def tag(self) -> str:
        return f"web-{self.number}"


def read_major() -> int:
    text = MAJOR_FILE.read_text().strip()
    if not text.isdigit() or int(text) < 1:
        raise SystemExit(f"{MAJOR_FILE} must hold a whole number from 1, not {text!r}")
    return int(text)


def next_version(releases: list[dict], major: int, fingerprint: str, today: str) -> Version:
    versioned = []
    for release in releases:
        match = TAG.match(str(release.get("tag_name", "")))
        if match and not release.get("draft"):
            versioned.append(((int(match.group(1)), int(match.group(2))), release))
    if versioned:
        (newest_major, newest_minor), newest = max(versioned, key=lambda item: item[0])
        if f"{FINGERPRINT_LINE}{fingerprint}" in (newest.get("body") or ""):
            date = str(newest.get("published_at") or today)[:10]
            return Version(f"{newest_major}.{newest_minor}", date, new=False)
    minors = [minor for (m, minor), _ in versioned if m == major]
    return Version(f"{major}.{max(minors) + 1 if minors else 0}", today, new=True)


def changelog_entry(text: str, number: str) -> str | None:
    """The lines under "## <number>" (up to the next "## "), or None if absent or empty."""
    heading = re.compile(rf"^## {re.escape(number)}(?![\d.])[^\n]*$", re.MULTILINE)
    match = heading.search(text)
    if not match:
        return None
    rest = text[match.end() :]
    following = re.search(r"^## ", rest, re.MULTILINE)
    body = (rest[: following.start()] if following else rest).strip()
    return body or None


def stamp(site: Path, number: str, date: str) -> list[Path]:
    """Write the version into every page folder's version.js; returns the files written."""
    if not NUMBER.match(number):
        raise SystemExit(f"Not a version number: {number!r}")
    if not DATE.match(date):
        raise SystemExit(f"Not a YYYY-MM-DD date: {date!r}")
    written = []
    for folder in PAGE_FOLDERS:
        path = site / folder / "version.js"
        if not path.is_file():
            raise SystemExit(f"{path} is missing; run scripts/build_site.py first")
        text = path.read_text()
        if PLACEHOLDER not in text:
            raise SystemExit(f"{path} has no version placeholder (already stamped?)")
        stamped = f'const SITE_VERSION = {{ version: "{number}", date: "{date}" }};'
        path.write_text(text.replace(PLACEHOLDER, stamped))
        written.append(path)
    return written


def _releases(path: Path) -> list[dict]:
    pages = json.loads(path.read_text() or "[]")
    # `gh api --paginate --slurp` gives a list of pages; accept that or one flat list.
    return [r for page in pages for r in (page if isinstance(page, list) else [page])]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    nxt = commands.add_parser("next", help="the version for this build")
    nxt.add_argument("--fingerprint", required=True)
    nxt.add_argument("--releases", type=Path, required=True, help="JSON from the GitHub releases API")
    nxt.add_argument("--today", default=datetime.now(UTC).date().isoformat())
    notes = commands.add_parser("notes", help="print a version's changelog entry")
    notes.add_argument("--version", required=True)
    notes.add_argument("--changelog", type=Path, default=CHANGELOG)
    stp = commands.add_parser("stamp", help="write the version into a built site")
    stp.add_argument("--site", type=Path, required=True)
    stp.add_argument("--version", required=True)
    stp.add_argument("--date", required=True)
    args = parser.parse_args(argv)

    if args.command == "next":
        if not re.fullmatch(r"[0-9a-f]{64}", args.fingerprint):
            raise SystemExit(f"Not a fingerprint: {args.fingerprint!r}")
        if not DATE.match(args.today):
            raise SystemExit(f"Not a YYYY-MM-DD date: {args.today!r}")
        v = next_version(_releases(args.releases), read_major(), args.fingerprint, args.today)
        print(f"{v.number} {'new' if v.new else 'same'} {v.date}")
        return 0
    if args.command == "notes":
        entry = changelog_entry(args.changelog.read_text(), args.version)
        if entry is None:
            print(f"{args.changelog} has no entry for {args.version}", file=sys.stderr)
            return 1
        print(entry)
        return 0
    for path in stamp(args.site, args.version, args.date):
        print(f"Stamped {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
