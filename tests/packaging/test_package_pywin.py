"""Tests for scripts/package_pywin.sh (the student PyWin.tar bundle)."""

import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "package_pywin.sh"


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["zsh", str(SCRIPT), *args], cwd=ROOT, capture_output=True, text=True, timeout=120
    )


def test_bundle_is_minimal_and_under_pywin(tmp_path: Path) -> None:
    out = tmp_path / "PyWin.tar"
    result = run("--out", str(out))
    assert result.returncode == 0, result.stderr
    with tarfile.open(out) as tar:
        names = tar.getnames()
    assert names, "archive is empty"
    assert all(n == "PyWin" or n.startswith("PyWin/") for n in names), names[:5]
    for entry in ("ass.py", "StackImages.py", "ShowClusters.py", "GetFit.py"):
        assert f"PyWin/{entry}" in names
    assert "PyWin/csv/NI2021_MessierObjects.csv" in names
    assert "PyWin/fit/sxv_psf4-RGB.fit" in names
    for forbidden in ("web/", "tests/", "ci/", "scripts/", "docs/", ".github/"):
        assert not any(n.startswith(f"PyWin/{forbidden}") for n in names), forbidden


def test_dry_run_writes_nothing(tmp_path: Path) -> None:
    out = tmp_path / "PyWin.tar"
    result = run("--dry-run", "--out", str(out))
    assert result.returncode == 0, result.stderr
    assert "PyWin/ass.py" in result.stdout
    assert not out.exists()


def test_rerun_overwrites_cleanly(tmp_path: Path) -> None:
    out = tmp_path / "PyWin.tar"
    assert run("--out", str(out)).returncode == 0
    first = out.read_bytes()
    assert run("--out", str(out)).returncode == 0
    assert out.read_bytes() == first
    assert not (tmp_path / "PyWin.tar.partial").exists()


def test_bad_ref_fails_loudly(tmp_path: Path) -> None:
    result = run("--ref", "no-such-ref", "--out", str(tmp_path / "x.tar"))
    assert result.returncode != 0
    assert "Not a commit" in result.stderr
    assert not (tmp_path / "x.tar").exists()


def test_unknown_argument_fails() -> None:
    result = run("--bogus")
    assert result.returncode == 2
