"""Tests for scripts/next_release.py (desktop release versioning)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import next_release as nr  # noqa: E402

FP_A = "a" * 40
FP_B = "b" * 40


def release(tag: str, fingerprint: str = FP_A, draft: bool = False) -> dict:
    return {"tag_name": tag, "draft": draft, "body": f"Notes\n\n{nr.FINGERPRINT_LINE}{fingerprint}\n"}


def test_first_versioned_release_ignores_old_style_tags() -> None:
    old = [{"tag_name": "v-535325b2", "body": ""}, {"tag_name": "snapshots-20261102-2324", "body": ""}]
    assert nr.next_tag(old, 21, FP_A) == "v21.0"


def test_unchanged_bundle_means_no_release() -> None:
    assert nr.next_tag([release("v21.0", FP_A)], 21, FP_A) is None


def test_changed_bundle_increments_minor() -> None:
    releases = [release("v21.0", FP_A), release("v21.1", FP_A)]
    assert nr.next_tag(releases, 21, FP_B) == "v21.2"


def test_new_app_version_starts_at_zero() -> None:
    assert nr.next_tag([release("v21.3", FP_A)], 22, FP_B) == "v22.0"


def test_compares_with_newest_by_number_not_list_order() -> None:
    releases = [release("v21.10", FP_B), release("v21.9", FP_A)]
    assert nr.next_tag(releases, 21, FP_B) is None
    assert nr.next_tag(releases, 21, FP_A) == "v21.11"


def test_drafts_are_ignored() -> None:
    assert nr.next_tag([release("v21.0", FP_A, draft=True)], 21, FP_A) == "v21.0"


def run(tmp_path: Path, pages: object, fingerprint: str = FP_A) -> subprocess.CompletedProcess[str]:
    data = tmp_path / "releases.json"
    data.write_text(json.dumps(pages))
    return subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "next_release.py"),
            "--app-version",
            "21",
            "--fingerprint",
            fingerprint,
            "--releases",
            str(data),
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )


def test_cli_reads_slurped_pages(tmp_path: Path) -> None:
    result = run(tmp_path, [[release("v21.0", FP_B)], [release("v-abc")]])
    assert result.returncode == 0 and result.stdout.strip() == "v21.1"
    unchanged = run(tmp_path, [[release("v21.0", FP_A)]])
    assert unchanged.returncode == 0 and unchanged.stdout.strip() == ""


def test_cli_rejects_bad_fingerprint(tmp_path: Path) -> None:
    result = run(tmp_path, [], fingerprint="not-a-hash; rm -rf /")
    assert result.returncode != 0 and "Not a fingerprint" in result.stderr


@pytest.mark.parametrize("pages", [[], [[]]])
def test_cli_with_no_releases(tmp_path: Path, pages: object) -> None:
    assert run(tmp_path, pages).stdout.strip() == "v21.0"
