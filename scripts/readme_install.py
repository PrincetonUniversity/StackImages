"""Run the README's download and update commands as a student would, and check the result.

The commands are read from README.md itself (the macOS or Windows section), so CI tests the
text students copy. Windows runs them in cmd.exe (what Anaconda Prompt is), macOS in zsh.

    python scripts/readme_install.py show windows update   # print one command
    python scripts/readme_install.py check windows         # download, damage, update, verify

`check` uses ~/Downloads/StackImages and ~/Downloads/StackImages.tar, replacing what is there;
so it refuses to run unless CI=true (set on GitHub runners). Standard library only.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SECTIONS = {"macos": "macOS", "windows": "Windows"}
KINDS = {"download": " -xvf ", "update": " -xvkf "}  # tar flags tell the two commands apart
EXPECTED = ("StackImages.py", "GetFit.py", "ass.py")
LOCAL_FILE = "my_notes.txt"  # a student's own file: the update must keep it
DAMAGED = "ass.py"  # changed locally: the update must restore it

Runner = Callable[[str, str], None]


class CheckError(Exception):
    """The commands ran but left the wrong files, or a command failed."""


def section(text: str, title: str) -> str:
    """The body of a setext "Title\\n-----" section, up to the next such heading."""
    lines = text.splitlines()
    starts = [i for i in range(len(lines) - 1) if re.fullmatch(r"-{3,}|={3,}", lines[i + 1].strip())]
    for n, i in enumerate(starts):
        if lines[i].strip() == title:
            end = starts[n + 1] if n + 1 < len(starts) else len(lines)
            return "\n".join(lines[i + 2 : end])
    raise ValueError(f"README has no {title!r} section")


def code_blocks(text: str) -> list[list[str]]:
    """Fenced code blocks (also indented ones inside list items), as lists of stripped lines."""
    blocks: list[list[str]] = []
    current: list[str] | None = None
    for line in text.splitlines():
        if line.strip().startswith("```"):
            if current is None:
                current = []
            else:
                blocks.append(current)
                current = None
        elif current is not None and line.strip():
            current.append(line.strip())
    if current is not None:
        raise ValueError("unclosed code block")
    return blocks


def command(text: str, platform: str, kind: str) -> str:
    """The one-line download or update command from the platform's README section."""
    body = section(text, SECTIONS[platform])
    found = [b for b in code_blocks(body) if "PyWin.tar" in " ".join(b) and KINDS[kind] in f"{' '.join(b)} "]
    if len(found) != 1:
        raise ValueError(f"expected one {kind} command in {SECTIONS[platform]}, found {len(found)}")
    if len(found[0]) != 1:
        raise ValueError(f"the {SECTIONS[platform]} {kind} command must be a single line")
    return found[0][0]


def run_in_shell(platform: str, cmd: str) -> None:
    """Run cmd in the shell students paste into; raise CheckError on a non-zero exit."""
    if platform == "windows":
        # Windows' own tar.exe and curl.exe first, as in Anaconda Prompt: a runner's Git Bash
        # puts GNU tar ahead of them, which reads "C:" as a remote host.
        system32 = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32"
        env = {**os.environ, "PATH": f"{system32};{os.environ.get('PATH', '')}"}
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / "step.cmd"
            script.write_text(f"where tar curl\r\n{cmd}\r\nexit /b %errorlevel%\r\n", encoding="utf-8")
            result = subprocess.run(["cmd", "/d", "/c", str(script)], check=False, env=env)
    else:
        result = subprocess.run(["/bin/zsh", "-c", cmd], check=False)
    if result.returncode != 0:
        raise CheckError(f"command exited with {result.returncode}: {cmd}")


def bundle_files(tarball: Path) -> dict[str, bytes]:
    """Regular files under PyWin/ in the tar, by path below PyWin/ (AppleDouble ._ files left out)."""
    files: dict[str, bytes] = {}
    with tarfile.open(tarball) as tar:
        for member in tar.getmembers():
            parts = member.name.split("/")
            if not member.isfile() or parts[0] != "PyWin" or len(parts) < 2:
                continue
            if parts[-1].startswith("._"):
                continue
            data = tar.extractfile(member)
            assert data is not None
            files["/".join(parts[1:])] = data.read()
    if not files:
        raise CheckError(f"{tarball} has no files under PyWin/")
    return files


def verify(target: Path, tarball: Path, step: str) -> int:
    """Every bundle file is in target with the bundle's bytes; returns the number checked."""
    if not tarball.is_file():
        raise CheckError(f"{step}: {tarball} was not downloaded")
    files = bundle_files(tarball)
    for name in EXPECTED:
        if name not in files:
            raise CheckError(f"{step}: the bundle has no {name}")
    if (target / "PyWin").exists():
        raise CheckError(f"{step}: files landed in {target / 'PyWin'}, not {target}")
    wrong = [
        n for n, data in files.items() if not (target / n).is_file() or (target / n).read_bytes() != data
    ]
    if wrong:
        raise CheckError(f"{step}: {len(wrong)} file(s) missing or different, e.g. {sorted(wrong)[:5]}")
    apple = [p for p in target.rglob("._*")]
    if apple:
        raise CheckError(f"{step}: AppleDouble files were extracted, e.g. {apple[:3]}")
    return len(files)


def check(text: str, platform: str, home: Path, run: Runner = run_in_shell) -> None:
    """Download, damage the install, update twice; raise CheckError at the first problem."""
    download, update = command(text, platform, "download"), command(text, platform, "update")
    downloads = home / "Downloads"
    target, tarball = downloads / "StackImages", downloads / "StackImages.tar"
    downloads.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(target, ignore_errors=True)
    tarball.unlink(missing_ok=True)

    run(platform, download)
    print(f"download: {verify(target, tarball, 'download')} files match the bundle")

    (target / DAMAGED).write_text("stale\n")
    (target / LOCAL_FILE).write_text("mine\n")
    for attempt in ("update", "second update"):
        run(platform, update)
        count = verify(target, tarball, attempt)
        if not (target / LOCAL_FILE).is_file() or (target / LOCAL_FILE).read_text() != "mine\n":
            raise CheckError(f"{attempt}: removed or changed the student's own {LOCAL_FILE}")
        print(f"{attempt}: {count} files match the bundle; {DAMAGED} restored, {LOCAL_FILE} kept")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--readme", type=Path, default=ROOT / "README.md")
    sub = parser.add_subparsers(dest="action", required=True)
    show = sub.add_parser("show", help="print one command")
    show.add_argument("platform", choices=SECTIONS)
    show.add_argument("kind", choices=KINDS)
    run_check = sub.add_parser("check", help="run the commands and verify the files")
    run_check.add_argument("platform", choices=SECTIONS)
    args = parser.parse_args(argv)
    text = args.readme.read_text(encoding="utf-8")
    try:
        if args.action == "show":
            print(command(text, args.platform, args.kind))
        elif os.environ.get("CI") != "true":
            raise CheckError("check replaces ~/Downloads/StackImages; it runs only where CI=true")
        else:
            check(text, args.platform, Path.home())
            print(f"OK: the README's {SECTIONS[args.platform]} download and update commands work")
    except (ValueError, CheckError) as exc:
        print(f"FAILED: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
