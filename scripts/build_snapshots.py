"""Save ESA Gaia DR3 results for common clusters: the page's fallback when GAVO is down.

Each snapshot is the exact query ShowClusters sends (ESA's gaiadr3.gaia_source), centred
on the catalog position as the page displays it, for a generous field and magnitude
limit. The page filters it down to the student's field. Output (not committed):
  <out>/index.json and <out>/<NAME>.csv

Re-runnable: valid existing snapshots are kept. Uses ESA's asynchronous job queue
(the synchronous endpoint rejects work when busy) with retries.

Usage:
  .venv/bin/python scripts/build_snapshots.py --out build/snapshots [--only M13 M44] [--dry-run]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "web"))

import showclusters as sc  # noqa: E402

COLUMNS = list(sc.REQUIRED)
MAX_SNAPSHOT_ROWS = 60_000  # keeps each file to a few MB; magnitude limit drops until it fits
GLOBULAR = (2, 3, 4, 5, 9, 10, 12, 13, 14, 15, 19, 22, 28, 30, 53, 54, 55, 56, 62, 68, 69, 70, 71, 72)
GLOBULAR += (75, 79, 80, 92, 107)
OPEN = (6, 7, 11, 18, 21, 23, 25, 26, 29, 34, 35, 36, 37, 38, 39, 41, 46, 47, 48, 50, 52, 67, 93, 103)
OTHER = ("NGC104", "NGC5139", "NGC6397", "NGC6752", "NGC5466", "NGC188", "NGC869", "NGC884", "NGC2158")
OTHER += ("NGC2516", "NGC3532", "NGC6791", "NGC6819", "NGC7789", "NGC752", "NGC457", "IC4665", "IC2602")
OTHER += ("IC2391",)

# name → (field of view in degrees, magnitude limit). Big, bright clusters get wide fields.
PLAN: dict[str, tuple[float, float]] = {f"M{n}": (0.5, 18.0) for n in GLOBULAR}
PLAN.update({f"M{n}": (1.0, 17.0) for n in OPEN})
PLAN.update({name: (1.0, 17.0) for name in OTHER})
PLAN.update({"M44": (2.0, 18.0), "M45": (2.0, 18.0), "NGC5139": (1.0, 17.0), "NGC104": (0.6, 17.0)})


def query_esa(adql: str, attempts: int = 5) -> list[dict[str, str]]:
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = -1
    for attempt in range(1, attempts + 1):
        try:
            job = Gaia.launch_job_async(adql, output_format="csv", dump_to_file=False)
            table = job.get_results()
            buffer = io.StringIO()
            table.write(buffer, format="ascii.csv")
            return list(csv.DictReader(io.StringIO(buffer.getvalue())))
        except Exception as error:  # ESA outages surface as many exception types
            if attempt == attempts:
                raise
            wait = 30 * attempt
            print(f"  ESA error ({type(error).__name__}: {error}); retrying in {wait}s", flush=True)
            time.sleep(wait)
    raise AssertionError("unreachable")


def build_one(name: str, out: Path, today: str) -> dict:
    fov, maxmag = PLAN[name]
    position = sc.lookup_object(name, ROOT / "csv")
    ra_text, dec_text = sc.format_sexagesimal(position.ra_hours), sc.format_sexagesimal(position.dec_deg)
    ra, dec = sc.parse_sexagesimal(ra_text), sc.parse_sexagesimal(dec_text)
    while True:
        rows = query_esa(sc.gaia_adql(ra, dec, fov, maxmag, "esa"))
        if len(rows) <= MAX_SNAPSHOT_ROWS or maxmag <= 12:
            break
        print(f"  {len(rows):,} rows at G < {maxmag}; lowering the limit", flush=True)
        maxmag -= 1
    path = out / f"{position.name}.csv"
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    parsed = sc.parse_gaia_table(path.read_bytes(), path.name)  # same reader the page uses
    return {
        "file": path.name,
        "ra": ra_text,
        "dec": dec_text,
        "fov": fov,
        "maxmag": maxmag,
        "rows": len(parsed["ra"]),
        "date": today,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=ROOT / "build" / "snapshots")
    parser.add_argument("--only", nargs="+", metavar="NAME", help="subset of clusters, e.g. M13 M44")
    parser.add_argument("--dry-run", action="store_true", help="list what would be queried")
    args = parser.parse_args(argv)
    names = [n.upper() for n in args.only] if args.only else list(PLAN)
    unknown = [n for n in names if n not in PLAN]
    if unknown:
        raise SystemExit(f"Not in the snapshot plan: {', '.join(unknown)}")
    if args.dry_run:
        for name in names:
            print(f"{name}: FOV {PLAN[name][0]}°, G < {PLAN[name][1]}")
        return 0

    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    index_path = out / "index.json"
    index: dict[str, dict] = json.loads(index_path.read_text()) if index_path.exists() else {}
    today = datetime.now(UTC).date().isoformat()
    failed: list[str] = []
    for name in names:
        entry = index.get(name)
        if entry and (out / entry["file"]).is_file():
            print(f"{name}: kept ({entry['rows']:,} rows from {entry['date']})")
            continue
        print(f"{name}: querying ESA…", flush=True)
        try:
            index[name] = build_one(name, out, today)
        except Exception as error:
            print(f"{name}: FAILED ({type(error).__name__}: {error})", flush=True)
            failed.append(name)
            continue
        index_path.write_text(json.dumps(index, indent=1, sort_keys=True) + "\n")  # resume point
        print(f"{name}: {index[name]['rows']:,} rows", flush=True)
    total = sum((out / e["file"]).stat().st_size for e in index.values())
    print(f"{len(index)} snapshots, {total / 1e6:.1f} MB in {out}")
    if failed:
        print(f"Failed: {', '.join(failed)}. Re-run to retry only these.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
