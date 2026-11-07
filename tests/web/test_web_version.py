"""The web site's version (scripts/web_version.py): web-MAJOR.N releases, the changelog, and
the label every page shows."""

from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import web_version as wv  # noqa: E402

FP_A, FP_B = "a" * 64, "b" * 64


def release(tag: str, fingerprint: str = FP_A, draft: bool = False, date: str = "2026-11-01") -> dict:
    body = f"Notes.\n\nSite fingerprint: {fingerprint}\n"
    return {"tag_name": tag, "body": body, "draft": draft, "published_at": f"{date}T12:00:00Z"}


# --- numbering ----------------------------------------------------------------------------


def test_first_version_is_major_dot_zero() -> None:
    v = wv.next_version([], 1, FP_A, "2026-11-07")
    assert (v.number, v.tag, v.date, v.new) == ("1.0", "web-1.0", "2026-11-07", True)


def test_unchanged_site_keeps_its_version_and_date() -> None:
    releases = [release("web-1.0", FP_B), release("web-1.1", FP_A, date="2026-11-05")]
    v = wv.next_version(releases, 1, FP_A, "2026-11-07")
    assert (v.number, v.date, v.new) == ("1.1", "2026-11-05", False)


def test_changed_site_is_the_next_version() -> None:
    releases = [release("web-1.0", FP_B), release("web-1.1", FP_B)]
    assert wv.next_version(releases, 1, FP_A, "2026-11-07").number == "1.2"


def test_versions_compare_as_numbers() -> None:
    releases = [release(f"web-1.{n}", FP_B) for n in (9, 10, 2)]
    assert wv.next_version(releases, 1, FP_A, "2026-11-07").number == "1.11"


def test_new_major_starts_at_zero() -> None:
    releases = [release("web-1.4", FP_B)]
    assert wv.next_version(releases, 2, FP_A, "2026-11-07").number == "2.0"


def test_dated_releases_drafts_and_other_tags_are_ignored() -> None:
    releases = [
        release("web-20261107-1441", FP_A),
        release("hrdiagram-20261107-1441", FP_A),
        release("v21.0", FP_A),
        release("web-1.3", FP_A, draft=True),
        release("web-1.0", FP_B),
    ]
    v = wv.next_version(releases, 1, FP_A, "2026-11-07")
    assert (v.number, v.new) == ("1.1", True)


def test_major_comes_from_the_version_file() -> None:
    assert wv.read_major() == int((ROOT / "web" / "VERSION").read_text().strip()) >= 1


# --- changelog ----------------------------------------------------------------------------

CHANGELOG = """# Changelog

## 1.10 (2026-12-01)

Ten.

## 1.1 (2026-11-08)

- One.
- Also one.

## 1.0 (2026-11-07)

First.
"""


def test_changelog_entry_is_the_section_for_that_version() -> None:
    assert wv.changelog_entry(CHANGELOG, "1.1") == "- One.\n- Also one."
    assert wv.changelog_entry(CHANGELOG, "1.0") == "First."
    assert wv.changelog_entry(CHANGELOG, "1.10") == "Ten."


def test_missing_changelog_entry_is_none() -> None:
    assert wv.changelog_entry(CHANGELOG, "1.2") is None
    assert wv.changelog_entry("# Changelog\n\n## 1.3\n\n", "1.3") is None  # heading without text


def test_repository_changelog_has_the_first_version() -> None:
    text = (ROOT / "docs" / "CHANGELOG.md").read_text()
    assert wv.changelog_entry(text, f"{wv.read_major()}.0")


# --- command line -------------------------------------------------------------------------


def test_next_command(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = tmp_path / "releases.json"
    path.write_text(json.dumps([[release("web-1.0", FP_B)]]))  # gh api --paginate: list of pages
    assert wv.main(["next", "--fingerprint", FP_A, "--releases", str(path), "--today", "2026-11-07"]) == 0
    assert capsys.readouterr().out.strip() == "1.1 new 2026-11-07"


def test_next_command_rejects_a_bad_fingerprint(tmp_path: Path) -> None:
    path = tmp_path / "releases.json"
    path.write_text("[]")
    with pytest.raises(SystemExit):
        wv.main(["next", "--fingerprint", "nope", "--releases", str(path), "--today", "2026-11-07"])


def test_notes_command(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    log = tmp_path / "CHANGELOG.md"
    log.write_text(CHANGELOG)
    assert wv.main(["notes", "--version", "1.1", "--changelog", str(log)]) == 0
    assert capsys.readouterr().out.strip() == "- One.\n- Also one."
    assert wv.main(["notes", "--version", "1.2", "--changelog", str(log)]) == 1  # missing: fails


# --- stamping -----------------------------------------------------------------------------


@pytest.fixture()
def built(tmp_path: Path) -> Path:
    """The three page folders as build_site.py leaves them, each with the unstamped version.js."""
    site = tmp_path / "site"
    for folder in ("", "ShowClusters", "HRDiagram"):
        (site / folder).mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "web" / "version.js", site / folder / "version.js")
    return site


def test_stamp_writes_the_version_into_every_page_folder(built: Path) -> None:
    written = wv.stamp(built, "1.3", "2026-11-07")
    assert sorted(p.relative_to(built).as_posix() for p in written) == [
        "HRDiagram/version.js",
        "ShowClusters/version.js",
        "version.js",
    ]
    for path in written:
        text = path.read_text()
        assert 'version: "1.3"' in text and 'date: "2026-11-07"' in text


@pytest.mark.parametrize(
    ("number", "date"), [("1.x", "2026-11-07"), ("1.3", "7 Nov"), ('1.3"; alert(1); "', "2026-11-07")]
)
def test_stamp_refuses_bad_values(built: Path, number: str, date: str) -> None:
    with pytest.raises(SystemExit):
        wv.stamp(built, number, date)


def test_stamp_refuses_a_site_without_version_files(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        wv.stamp(tmp_path, "1.0", "2026-11-07")


def test_unstamped_source_says_development_build() -> None:
    source = (ROOT / "web" / "version.js").read_text()
    assert wv.PLACEHOLDER in source and "Development build" in source


@pytest.mark.parametrize("page", ["web/home.html", "web/index.html", "web/hrdiagram/index.html"])
def test_every_page_shows_the_version_and_links_the_changelog(page: str) -> None:
    html = (ROOT / page).read_text()
    assert re.search(r'<script src="version\.js" defer></script>', html)
    assert "data-site-version" in html
    assert "docs/CHANGELOG.md" in html
