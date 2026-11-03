"""Download the pinned Pyodide runtime and the wheels the page needs; verify every file.

Output goes to build/vendor/ (gitignored; never committed). Re-running reuses verified
files. Core runtime hashes are pinned below; wheel hashes come from the pinned
pyodide-lock.json. Any mismatch deletes the file and fails.

Usage: .venv/bin/python scripts/vendor_pyodide.py [--dest build/vendor]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

VERSION = "0.29.3"
BASE = f"https://cdn.jsdelivr.net/pyodide/v{VERSION}/full/"
ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ("numpy", "matplotlib", "astropy")
CORE_SHA256 = {
    "pyodide.js": "718d40f1c015dd25ec724cc8fc4e2325d6a45a92ae225121ff6953f224a16f72",
    "pyodide.asm.js": "1263f02b5b26099b96112378156f242dd98b39a8201ba7765e5fe3d455c5ce91",
    "pyodide.asm.wasm": "e2f4ee75b325e35eb31bfb8c613d4dd5098f5502c156a97847686875b5025480",
    "python_stdlib.zip": "4298b6ee445cb724c3973437da47789752b9e6ff4e26619026b283ec801fc46b",
    "pyodide-lock.json": "3256ffc76388de0e37f4b34d42ab484268d1afc675179ff97b2a5bb14f84ccac",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(dest: Path, name: str, expected: str) -> None:
    target = dest / name
    if target.exists() and sha256(target) == expected:
        return
    with urllib.request.urlopen(BASE + name, timeout=300) as response:  # noqa: S310 - fixed https URL
        content = response.read()
    if hashlib.sha256(content).hexdigest() != expected:
        raise SystemExit(f"Checksum mismatch for {name}; refusing to use it.")
    partial = target.with_suffix(target.suffix + ".partial")
    partial.write_bytes(content)
    partial.replace(target)
    print(f"Downloaded {name}", flush=True)


def wheels(lock: dict, names: tuple[str, ...]) -> dict[str, str]:
    """File name → sha256 for the packages and all their dependencies."""
    packages = lock["packages"]
    found: dict[str, str] = {}
    seen: set[str] = set()

    def visit(name: str) -> None:
        key = name.lower().replace("_", "-")
        if key in seen:
            return
        seen.add(key)
        entry = packages[key]
        for dependency in entry["depends"]:
            visit(dependency)
        found[entry["file_name"]] = entry["sha256"]

    for name in names:
        visit(name)
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dest", type=Path, default=ROOT / "build" / "vendor")
    args = parser.parse_args(argv)
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    for name, expected in CORE_SHA256.items():
        fetch(dest, name, expected)
    lock = json.loads((dest / "pyodide-lock.json").read_text())
    needed = wheels(lock, PACKAGES)
    for name, expected in needed.items():
        fetch(dest, name, expected)
    keep = set(CORE_SHA256) | set(needed)
    for stale in [p for p in dest.iterdir() if p.is_file() and p.name not in keep]:
        stale.unlink()
    print(f"Pyodide {VERSION}: {len(CORE_SHA256)} runtime files, {len(needed)} wheels verified in {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
