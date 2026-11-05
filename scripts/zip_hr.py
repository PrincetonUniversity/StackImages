"""Pack "HR Diagram of Your Own Image" (_site/HRDiagram/) into one zip of its own.

Unzipping gives an HRDiagram/ folder: the page, its own copy of Pyodide (vendor/) and
HOSTING.txt. On the site the page shares Show Clusters' Pyodide (../ShowClusters/vendor/),
so the zip points the worker at vendor/ and the StackImages link at the public site;
nothing else changes, and none of Show Clusters' data is included. The zip is
reproducible: the same site gives the same bytes.

Usage:
  .venv/bin/python scripts/zip_hr.py [--site _site] [--out _site/HRDiagram-web.zip] [--dry-run]
  .venv/bin/python scripts/zip_hr.py --fingerprint _site/HRDiagram-web.zip

--fingerprint prints a hash of the zip's contents; pages.yml publishes an hrdiagram-*
release only when it changes.
"""

from __future__ import annotations

import argparse
import sys
import zipfile
from pathlib import Path

from zip_site import entry, fingerprint

ROOT = Path(__file__).resolve().parents[1]
FOLDER = "HRDiagram"
NAME = "HRDiagram-web.zip"
HOME = "https://princetonuniversity.github.io/StackImages/"
# Each must occur exactly once, so a change to the page can't silently leave the zip
# pointing outside its folder.
REWRITES = {
    "worker.js": [('"../ShowClusters/vendor/"', '"vendor/"')],
    "index.html": [('href="../"', f'href="{HOME}"')],
}

HOSTING = """\
HR Diagram of Your Own Image (StackImages): a copy to put on your own web site
==============================================================================

1. Upload this whole HRDiagram folder to any directory of your web site, unchanged.
   The page is HRDiagram/index.html. Nothing to install or run on the server: these
   are ordinary static files.

2. Open it through the web server (https://...), not from your disk (file:///...):
   browsers do not run the page's background worker from local files.

3. No server settings are needed. The page works whatever content types the server
   sends for its .wasm and .whl files (application/octet-stream is fine), as long as
   .js, .html and .css have their usual types.

Students' images are read by Python running in their browser and are never uploaded;
the page fetches nothing from other sites. The diagram is shown as a data: image, so
a Content-Security-Policy must allow img-src data:.

Built from https://github.com/PrincetonUniversity/StackImages
(the same page as https://princetonuniversity.github.io/StackImages/HRDiagram/).
"""


def page_files(site: Path) -> list[tuple[str, Path]]:
    """(name in the folder, file in the site) for the page and Pyodide."""
    page, vendor = site / "HRDiagram", site / "ShowClusters" / "vendor"
    if not (page / "index.html").is_file() or not (vendor / "pyodide.js").is_file():
        raise SystemExit(f"{site} is not a built site; run scripts/build_site.py first")
    files = [(p.relative_to(page).as_posix(), p) for p in sorted(page.rglob("*")) if p.is_file()]
    files += [
        (f"vendor/{p.relative_to(vendor).as_posix()}", p) for p in sorted(vendor.rglob("*")) if p.is_file()
    ]
    return files


def contents(name: str, path: Path) -> bytes:
    data = path.read_bytes()
    for old, new in REWRITES.get(name, []):
        text = data.decode()
        if text.count(old) != 1:
            raise SystemExit(
                f"{FOLDER}/{name}: expected {old} once, found it {text.count(old)} times; "
                "update REWRITES in scripts/zip_hr.py"
            )
        data = text.replace(old, new).encode()
    return data


def write_zip(site: Path, out: Path) -> int:
    files = page_files(site)
    partial = out.with_name(out.name + ".partial")
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(partial, "w", compresslevel=9) as archive:
            archive.writestr(entry("HOSTING.txt", FOLDER), HOSTING)
            for name, path in files:
                archive.writestr(entry(name, FOLDER), contents(name, path))
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    partial.replace(out)  # never leave a half-written zip under the real name
    return len(files) + 1


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
        files = page_files(args.site)
        for name, _ in files:
            print(f"{FOLDER}/{name}")
        print(f"{FOLDER}/HOSTING.txt\n{len(files) + 1} files would go into {out}")
        return 0
    count = write_zip(args.site, out)
    print(f"Wrote {out}: {count} files, {out.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
