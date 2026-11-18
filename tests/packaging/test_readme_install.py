"""scripts/readme_install.py: finding the README's commands and checking what they leave behind.

The real commands run on Windows and macOS runners (install-commands.yml); here a fake shell
extracts a small bundle, so every way an update can go wrong is shown to be caught.
"""

from __future__ import annotations

import io
import shutil
import sys
import tarfile
from collections.abc import Callable
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import readme_install as ri  # noqa: E402

README = (ROOT / "README.md").read_text(encoding="utf-8")
BUNDLE = {"StackImages.py": b"gui", "GetFit.py": b"getfit", "ass.py": b"science", "csv/a.csv": b"1,2"}


def make_tar(path: Path, files: dict[str, bytes], extra_apple: bool = True) -> None:
    with tarfile.open(path, "w") as tar:
        names = {f"PyWin/{n}": d for n, d in files.items()}
        if extra_apple:
            names["PyWin/._ass.py"] = b"apple"
        for name, data in names.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def extract(tarball: Path, target: Path, *, keep_old: bool = False, strip: bool = True) -> None:
    with tarfile.open(tarball) as tar:
        for member in tar.getmembers():
            name = member.name.split("/", 1)[1] if strip else member.name
            if name.rsplit("/", 1)[-1].startswith("._"):
                continue
            dest = target / name
            if keep_old and dest.exists():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            data = tar.extractfile(member)
            assert data is not None
            dest.write_bytes(data.read())


def fake_shell(home: Path, source: Path, update: str = "good") -> ri.Runner:
    """A stand-in for cmd.exe/zsh: 'downloads' source, then extracts like the real flags would."""
    downloads = home / "Downloads"
    target, tarball = downloads / "StackImages", downloads / "StackImages.tar"

    def run(platform: str, cmd: str) -> None:
        assert "PyWin.tar" in cmd
        shutil.copy(source, tarball)
        if " -xvf " in cmd:
            target.mkdir(exist_ok=True)
            extract(tarball, target)
        elif update == "fails":
            raise ri.CheckError("command exited with 1")
        elif update == "keep_old":
            extract(tarball, target, keep_old=True)
        elif update == "wipe":
            shutil.rmtree(target)
            target.mkdir()
            extract(tarball, target)
        elif update == "nested":
            extract(tarball, target, strip=False)
        else:
            extract(tarball, target)

    return run


@pytest.fixture
def bundle(tmp_path: Path) -> Path:
    path = tmp_path / "PyWin.tar"
    make_tar(path, BUNDLE)
    return path


# --- reading the README ---------------------------------------------------------------------


@pytest.mark.parametrize("platform", ["macos", "windows"])
def test_real_readme_has_one_line_download_and_update(platform: str) -> None:
    download = ri.command(README, platform, "download")
    update = ri.command(README, platform, "update")
    assert "https://vanderbei.princeton.edu/tmp/PyWin.tar" in download
    assert "https://vanderbei.princeton.edu/tmp/PyWin.tar" in update
    assert "--strip-components=1" in download and "--strip-components=1" in update
    assert download != update


def test_windows_commands_are_cmd_syntax() -> None:
    for kind in ri.KINDS:
        cmd = ri.command(README, "windows", kind)
        assert "%USERPROFILE%\\Downloads\\StackImages" in cmd
        assert "$HOME" not in cmd and "\\\n" not in cmd and "'" not in cmd  # cmd keeps quotes literally


def test_macos_commands_use_home() -> None:
    for kind in ri.KINDS:
        assert "$HOME/Downloads/StackImages" in ri.command(README, "macos", kind)


def test_missing_section_is_an_error() -> None:
    with pytest.raises(ValueError, match="no 'Windows' section"):
        ri.command("Title\n=====\n\nmacOS\n-----\n", "windows", "update")


def test_duplicate_command_is_an_error() -> None:
    block = "```\ncurl x/PyWin.tar && tar -xvkf y\n```\n"
    with pytest.raises(ValueError, match="found 2"):
        ri.command(f"Windows\n-------\n{block}{block}", "windows", "update")


def test_missing_command_is_an_error() -> None:
    with pytest.raises(ValueError, match="found 0"):
        ri.command("Windows\n-------\n```\necho hi\n```\n", "windows", "download")


def test_multi_line_command_is_an_error() -> None:
    text = "macOS\n-----\n   ```bash\n   curl x/PyWin.tar && \\\n   tar -xvf y\n   ```\n"
    with pytest.raises(ValueError, match="single line"):
        ri.command(text, "macos", "download")


def test_unclosed_block_is_an_error() -> None:
    with pytest.raises(ValueError, match="unclosed"):
        ri.command("Windows\n-------\n```\ncurl x/PyWin.tar -xvf\n", "windows", "download")


def test_section_stops_at_next_heading() -> None:
    text = (
        "macOS\n-----\n```\ncurl a/PyWin.tar -xvf b\n```\n"
        "Windows\n-------\n```\ncurl c/PyWin.tar -xvf d\n```\n"
    )
    assert ri.command(text, "macos", "download") == "curl a/PyWin.tar -xvf b"
    assert ri.command(text, "windows", "download") == "curl c/PyWin.tar -xvf d"


def test_show_prints_the_command(capsys: pytest.CaptureFixture[str]) -> None:
    assert ri.main(["show", "windows", "update"]) == 0
    assert capsys.readouterr().out.strip() == ri.command(README, "windows", "update")


# --- checking the result ------------------------------------------------------------------


def test_check_passes_when_update_restores_and_keeps(
    tmp_path: Path, bundle: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    home = tmp_path / "home"
    ri.check(README, "windows", home, fake_shell(home, bundle))
    out = capsys.readouterr().out
    assert "download: 4 files match" in out
    assert "second update: 4 files match the bundle; ass.py restored, my_notes.txt kept" in out
    target = home / "Downloads" / "StackImages"
    assert (target / "ass.py").read_bytes() == b"science"
    assert (target / "my_notes.txt").read_text() == "mine\n"


def test_check_replaces_a_previous_install(tmp_path: Path, bundle: Path) -> None:
    home = tmp_path / "home"
    old = home / "Downloads" / "StackImages"
    old.mkdir(parents=True)
    (old / "ass.py").write_text("old")
    ri.check(README, "macos", home, fake_shell(home, bundle))
    assert (old / "ass.py").read_bytes() == b"science"


@pytest.mark.parametrize(
    ("update", "message"),
    [
        ("keep_old", r"update: 1 file\(s\) missing or different"),
        ("wipe", "removed or changed the student's own my_notes.txt"),
        ("nested", "landed in"),
        ("fails", "exited with 1"),
    ],
)
def test_check_catches_a_bad_update(tmp_path: Path, bundle: Path, update: str, message: str) -> None:
    home = tmp_path / "home"
    with pytest.raises(ri.CheckError, match=message):
        ri.check(README, "windows", home, fake_shell(home, bundle, update))


def test_check_catches_extracted_apple_double_files(tmp_path: Path, bundle: Path) -> None:
    home = tmp_path / "home"
    base = fake_shell(home, bundle)

    def leaky(platform: str, cmd: str) -> None:
        base(platform, cmd)
        (home / "Downloads" / "StackImages" / "._ass.py").write_bytes(b"apple")

    with pytest.raises(ri.CheckError, match="AppleDouble"):
        ri.check(README, "macos", home, leaky)


def test_check_catches_a_bundle_without_the_app(tmp_path: Path) -> None:
    source = tmp_path / "PyWin.tar"
    make_tar(source, {"README.md": b"docs"})
    home = tmp_path / "home"
    with pytest.raises(ri.CheckError, match="has no StackImages.py"):
        ri.check(README, "macos", home, fake_shell(home, source))


def test_check_catches_a_missing_download(tmp_path: Path) -> None:
    home = tmp_path / "home"
    noop: Callable[[str, str], None] = lambda platform, cmd: None  # noqa: E731
    with pytest.raises(ri.CheckError, match="was not downloaded"):
        ri.check(README, "windows", home, noop)


def test_main_reports_failure(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    readme = tmp_path / "README.md"
    readme.write_text("nothing here\n")
    assert ri.main(["--readme", str(readme), "show", "macos", "update"]) == 1
    assert "FAILED: README has no 'macOS' section" in capsys.readouterr().err


def test_check_refuses_outside_ci(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    called: list[str] = []
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(ri, "check", lambda *args: called.append("check"))
    assert ri.main(["check", "macos"]) == 1
    assert called == []
    assert "runs only where CI=true" in capsys.readouterr().err


def test_check_runs_in_ci(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    called: list[str] = []
    monkeypatch.setenv("CI", "true")
    monkeypatch.setattr(ri, "check", lambda text, platform, home: called.append(platform))
    assert ri.main(["check", "windows"]) == 0
    assert called == ["windows"]
    assert "OK: the README's Windows download and update commands work" in capsys.readouterr().out
