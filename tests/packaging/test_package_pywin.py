"""Tests for scripts/package_pywin.sh (the student PyWin.tar bundle)."""

import os
import shutil
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "package_pywin.sh"


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["zsh", str(SCRIPT), *args], cwd=ROOT, capture_output=True, text=True, timeout=120)


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


def make_repo(tmp_path: Path) -> Path:
    """A throwaway repo with the bundle's paths, so history can be changed freely."""
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    shutil.copy2(SCRIPT, repo / "scripts" / SCRIPT.name)
    (repo / "ass.py").write_bytes(b"global version\r\nversion = 21\r\n")
    for name in (
        "GetFit.py",
        "ShowClusters.py",
        "StackImages.py",
        "bias_file.txt",
        "flat_file.txt",
        "LICENSE",
    ):
        (repo / name).write_text(name)
    for folder in ("csv", "fit", "png", "txt"):
        (repo / folder).mkdir()
        (repo / folder / "a").write_text(folder)
    (repo / "README.md").write_text("readme")
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "initial")
    return repo


def git(repo: Path, *args: str) -> str:
    env = {
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@example.com",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@example.com",
        "PATH": os.environ["PATH"],
        "HOME": str(repo),
    }
    return subprocess.run(
        ["git", *args], cwd=repo, env=env, check=True, capture_output=True, text=True
    ).stdout


def fingerprint(repo: Path) -> str:
    result = subprocess.run(
        ["zsh", str(repo / "scripts" / SCRIPT.name), "--fingerprint"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_fingerprint_tracks_app_files_not_docs(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    first = fingerprint(repo)
    assert len(first) >= 40
    (repo / "README.md").write_text("new docs")
    (repo / "LICENSE").write_text("new license")
    git(repo, "commit", "-qam", "docs only")
    assert fingerprint(repo) == first
    (repo / "csv" / "a").write_text("new data")
    git(repo, "commit", "-qam", "data")
    assert fingerprint(repo) != first


def test_app_version_reads_ass_py_with_crlf(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    script = repo / "scripts" / SCRIPT.name
    result = subprocess.run(["zsh", str(script), "--app-version"], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0 and result.stdout.strip() == "21"
    (repo / "ass.py").write_text("no version here\n")
    git(repo, "commit", "-qam", "drop version")
    result = subprocess.run(["zsh", str(script), "--app-version"], capture_output=True, text=True, timeout=60)
    assert result.returncode == 1 and "No 'version = N' line" in result.stderr
