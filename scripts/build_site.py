"""Assemble the GitHub Pages site in _site/ from web/, build/vendor/ and repo data.

Generated data (never committed):
  data/csv/            Messier, NGC and IC position catalogs (copied)
  data/bright_stars.npy  rows of csv/hygdata_v3.csv that ShowClusters can use
  data/psf.npy         the star image from fit/sxv_psf4-RGB.fit, scaled for the browser
  data/snapshots/      saved ESA results (from scripts/build_snapshots.py), if given

Usage:
  .venv/bin/python scripts/build_site.py --snapshots DIR   # normal build
  .venv/bin/python scripts/build_site.py --no-snapshots    # local testing only
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "web"))

import showclusters as sc  # noqa: E402

WEB_FILES = ("index.html", "app.js", "worker.js", "style.css", "showclusters.py")
CATALOGS = ("NI2021_MessierObjects.csv", "NI2021_NGCObjects.csv", "ICObj.csv")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_snapshots(source: Path) -> dict[str, dict]:
    """Validate a snapshot directory before publishing it. Fails loudly on any problem."""
    index_path = source / "index.json"
    if not index_path.is_file():
        raise SystemExit(f"{index_path} is missing")
    index = json.loads(index_path.read_text())
    if not isinstance(index, dict) or not index:
        raise SystemExit("Snapshot index is empty")
    for name, entry in index.items():
        missing = {"file", "ra", "dec", "fov", "maxmag", "date", "rows"} - set(entry)
        if missing:
            raise SystemExit(f"Snapshot {name} lacks {sorted(missing)}")
        path = source / entry["file"]
        if path.name != entry["file"] or not path.is_file():
            raise SystemExit(f"Snapshot {name}: bad or missing file {entry['file']!r}")
        table = sc.parse_gaia_table(path.read_bytes(), path.name)
        if len(table["ra"]) != entry["rows"]:
            raise SystemExit(f"Snapshot {name}: {len(table['ra'])} rows, index says {entry['rows']}")
    return index


def build(out: Path, vendor: Path, snapshots: Path | None) -> None:
    if not (vendor / "pyodide.mjs").is_file():
        raise SystemExit(f"{vendor} has no Pyodide runtime; run scripts/vendor_pyodide.py first")
    index = check_snapshots(snapshots) if snapshots else {}
    if out.exists():
        shutil.rmtree(out)
    (out / "data" / "csv").mkdir(parents=True)

    for name in WEB_FILES:
        shutil.copy2(ROOT / "web" / name, out / name)
    shutil.copytree(vendor, out / "vendor", ignore=shutil.ignore_patterns("*.partial"))
    for name in CATALOGS:
        shutil.copy2(ROOT / "csv" / name, out / "data" / "csv" / name)

    bright = np.array(sc.bright_star_rows(ROOT / "csv" / "hygdata_v3.csv"), dtype=np.float64)
    np.save(out / "data" / "bright_stars.npy", bright)
    psf = sc.scale_psf(sc.prepare_psf(fits.getdata(ROOT / "fit" / "sxv_psf4-RGB.fit")), sc.RENDER_SIZE)
    np.save(out / "data" / "psf.npy", psf)

    snapshot_dir = out / "data" / "snapshots"
    snapshot_dir.mkdir()
    for entry in index.values():
        shutil.copy2(snapshots / entry["file"], snapshot_dir / entry["file"])  # type: ignore[operator]
    (snapshot_dir / "index.json").write_text(json.dumps(index, indent=1, sort_keys=True) + "\n")

    # Cache-bust our own files so a redeploy is picked up immediately.
    page = (out / "index.html").read_text()
    for name in ("app.js", "style.css"):
        page = page.replace(f'"{name}"', f'"{name}?v={digest(out / name)[:12]}"')
    (out / "index.html").write_text(page)
    app = (out / "app.js").read_text()
    worker_version = digest(out / "worker.js")[:12]
    (out / "app.js").write_text(app.replace('new URL("worker.js"', f'new URL("worker.js?v={worker_version}"'))
    (out / ".nojekyll").write_text("")

    files = {
        p.relative_to(out).as_posix(): digest(p)
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name != "manifest.json"
    }
    manifest = {
        "built": datetime.now(UTC).isoformat(timespec="seconds"),
        "snapshots": len(index),
        "files": files,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f"Built {out}: {len(files)} files, {size / 1e6:.1f} MB, {len(index)} snapshots")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=ROOT / "_site")
    parser.add_argument("--vendor", type=Path, default=ROOT / "build" / "vendor")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--snapshots", type=Path, help="directory from scripts/build_snapshots.py")
    group.add_argument("--no-snapshots", action="store_true", help="omit the saved-copy fallback")
    args = parser.parse_args(argv)
    build(args.out, args.vendor, None if args.no_snapshots else args.snapshots)
    return 0


if __name__ == "__main__":
    sys.exit(main())
