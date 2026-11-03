"""Pack the built site (_site/) into one zip for hosting on any static web server.

Unzipping gives a StackImages/ folder: the home page, ShowClusters/ (the page, Pyodide
and the saved snapshots) and HOSTING.txt. Every link inside is relative, so the folder
works from any directory of any web server. The zip is reproducible: the same site
gives the same bytes.

Usage:
  .venv/bin/python scripts/zip_site.py [--site _site] [--out _site/StackImages-web.zip] [--dry-run]
  .venv/bin/python scripts/zip_site.py --fingerprint _site/StackImages-web.zip

--fingerprint prints a hash of the zip's contents, leaving out the build time recorded in
manifest.json; pages.yml publishes a web-* release only when it changes.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FOLDER = "StackImages"
NAME = "StackImages-web.zip"
SKIP = {".nojekyll", NAME}  # GitHub Pages only; the zip never contains itself
STORED = {".whl", ".zip", ".gz"}  # already compressed
TIMESTAMP = (1980, 1, 1, 0, 0, 0)  # fixed, so rebuilding an unchanged site gives the same zip
BUILD_RECORD = f"{FOLDER}/ShowClusters/manifest.json"  # holds the build time; its file list is
# made of hashes of the other files, so leaving it out loses nothing

HOSTING = """\
Show Clusters (StackImages): a copy to put on your own web site
================================================================

1. Upload this whole StackImages folder to any directory of your web site, unchanged.
   The home page is StackImages/index.html; Show Clusters is StackImages/ShowClusters/.
   Nothing to install or run on the server: these are ordinary static files.

2. Open it through the web server (https://...), not from your disk (file:///...):
   browsers do not run the page's background worker from local files.

3. No server settings are needed. The page works whatever content types the server
   sends for its .wasm, .whl, .npy and data files (application/octet-stream is fine),
   as long as .js, .html and .css have their usual types.

The page reads star data from GAVO (https://dc.g-vo.org) in the visitor's browser, so a
Content-Security-Policy must allow connect-src https://dc.g-vo.org (https: covers it).

Built from https://github.com/PrincetonUniversity/StackImages
(the same page as https://princetonuniversity.github.io/StackImages/ShowClusters/).
Data credits and licences are listed at the bottom of the page.
"""


def site_files(site: Path) -> list[Path]:
    if not (site / "index.html").is_file() or not (site / "ShowClusters" / "index.html").is_file():
        raise SystemExit(f"{site} is not a built site; run scripts/build_site.py first")
    return sorted(p for p in site.rglob("*") if p.is_file() and p.name not in SKIP)


def entry(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(f"{FOLDER}/{name}", TIMESTAMP)
    info.external_attr = 0o644 << 16  # rw-r--r--, whatever the build machine's umask
    info.create_system = 3  # Unix, so those permissions are honoured
    info.compress_type = zipfile.ZIP_STORED if Path(name).suffix in STORED else zipfile.ZIP_DEFLATED
    return info


def write_zip(site: Path, out: Path) -> int:
    files = site_files(site)
    partial = out.with_name(out.name + ".partial")
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(partial, "w", compresslevel=9) as archive:
        archive.writestr(entry("HOSTING.txt"), HOSTING)
        for path in files:
            archive.writestr(entry(path.relative_to(site).as_posix()), path.read_bytes())
    partial.replace(out)  # never leave a half-written zip under the real name
    return len(files) + 1


def fingerprint(archive: Path) -> str:
    """Hash of every file in the zip (names and contents) except the build record."""
    digest = hashlib.sha256()
    with zipfile.ZipFile(archive) as z:
        for name in sorted(z.namelist()):
            if name != BUILD_RECORD:
                digest.update(f"{name}\0{hashlib.sha256(z.read(name)).hexdigest()}\n".encode())
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--site", type=Path, default=ROOT / "_site")
    parser.add_argument("--out", type=Path, help=f"default: <site>/{NAME}")
    parser.add_argument("--dry-run", action="store_true", help="list what would be packed")
    parser.add_argument("--fingerprint", type=Path, metavar="ZIP", help="print a zip's content hash")
    args = parser.parse_args(argv)
    if args.fingerprint:
        print(fingerprint(args.fingerprint))
        return 0
    out: Path = args.out or args.site / NAME
    if args.dry_run:
        files = site_files(args.site)
        for path in files:
            print(f"{FOLDER}/{path.relative_to(args.site).as_posix()}")
        print(f"{FOLDER}/HOSTING.txt\n{len(files) + 1} files would go into {out}")
        return 0
    count = write_zip(args.site, out)
    print(f"Wrote {out}: {count} files, {out.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
