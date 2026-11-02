"""Save Gaia DR3 results for common clusters: the page's fallback when GAVO is down.

Each snapshot is the query ShowClusters sends, centred on the catalog position exactly as
the page displays it, for a generous field and magnitude limit; the page filters it down
to the student's field. Output (never committed; published as a snapshots-* release):
  <out>/index.json and <out>/<NAME>.csv  (columns: showclusters.REQUIRED)

Sources: ESA's archive (gaiadr3.gaia_source, what the desktop app uses) through its
asynchronous job queue, or GAVO's copy (gaia.dr3lite). `auto` tries ESA and falls back
to GAVO per cluster. Re-runnable: valid existing snapshots are kept.

Usage:
  .venv/bin/python scripts/build_snapshots.py [--out DIR] [--source auto|esa|gavo]
                                              [--only M13 M44] [--refresh] [--dry-run]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from collections.abc import Callable, Mapping
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "web"))

import showclusters as sc  # noqa: E402

GAVO_TAP = "https://dc.g-vo.org/tap/sync"
MAX_SNAPSHOT_ROWS = 60_000  # keeps each file to a few MB; the magnitude limit drops until it fits
ATTEMPTS = 4
GLOBULAR: tuple[int, ...] = (2, 3, 4, 5, 9, 10, 12, 13, 14, 15, 19, 22, 28, 30, 53, 54, 55, 56, 62, 68)
GLOBULAR += (69, 70, 71, 72, 75, 79, 80, 92, 107)
OPEN: tuple[int, ...] = (6, 7, 11, 18, 21, 23, 25, 26, 29, 34, 35, 36, 37, 38, 39, 41, 46, 47, 48, 50)
OPEN += (52, 67, 93, 103)
OTHER: tuple[str, ...] = ("NGC104", "NGC5139", "NGC6397", "NGC6752", "NGC5466", "NGC188", "NGC869")
OTHER += ("NGC884", "NGC2158", "NGC2516", "NGC3532", "NGC6791", "NGC6819", "NGC7789", "NGC752", "NGC457")
OTHER += ("IC4665", "IC2602", "IC2391")

# name → (field of view in degrees, magnitude limit). Big, bright clusters get wide fields.
PLAN: dict[str, tuple[float, float]] = {f"M{n}": (0.5, 18.0) for n in GLOBULAR}
PLAN.update({f"M{n}": (1.0, 17.0) for n in OPEN})
PLAN.update({name: (1.0, 17.0) for name in OTHER})
PLAN.update({"M44": (2.0, 18.0), "M45": (2.0, 18.0), "NGC104": (0.6, 17.0)})

Fetch = Callable[[str], bytes]


def fetch_esa(adql: str) -> bytes:
    """Run the ESA query as an asynchronous job; return CSV bytes."""
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = -1
    table = Gaia.launch_job_async(adql, dump_to_file=False).get_results()
    buffer = io.StringIO()
    table.write(buffer, format="ascii.csv")
    return buffer.getvalue().encode()


def fetch_gavo(adql: str) -> bytes:
    body = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": str(sc.MAX_ROWS), "QUERY": adql}
    ).encode()
    with urllib.request.urlopen(GAVO_TAP, body, timeout=300) as response:  # noqa: S310 - fixed https URL
        if "csv" not in response.headers.get("content-type", ""):
            raise RuntimeError(f"GAVO returned {response.headers.get('content-type')}")
        return bytes(response.read())


FETCHERS: dict[str, tuple[str, Fetch]] = {"esa": ("esa", fetch_esa), "gavo": ("gavo", fetch_gavo)}


def with_retries(fetch: Fetch, adql: str, sleep: Callable[[float], None] = time.sleep) -> bytes:
    for attempt in range(1, ATTEMPTS + 1):
        try:
            return fetch(adql)
        except Exception as error:  # outages surface as many exception types
            if attempt == ATTEMPTS:
                raise
            wait = 30 * attempt
            print(f"  {type(error).__name__}: {str(error)[:200]}; retrying in {wait}s", flush=True)
            sleep(wait)
    raise AssertionError("unreachable")


def build_one(
    name: str,
    out: Path,
    today: str,
    sources: list[str],
    fetchers: Mapping[str, tuple[str, Fetch]] = FETCHERS,
    sleep: Callable[[float], None] = time.sleep,
) -> dict:
    fov, maxmag = PLAN[name]
    position = sc.lookup_object(name, ROOT / "csv")
    ra_text, dec_text = sc.format_sexagesimal(position.ra_hours), sc.format_sexagesimal(position.dec_deg)
    ra, dec = sc.parse_sexagesimal(ra_text), sc.parse_sexagesimal(dec_text)
    errors: list[str] = []
    for source in sources:
        service, fetch = fetchers[source]
        limit = maxmag
        try:
            while True:
                data = with_retries(fetch, sc.gaia_adql(ra, dec, fov, limit, service), sleep)
                table = sc.parse_gaia_table(data, f"{name}.csv")
                if len(table["ra"]) <= MAX_SNAPSHOT_ROWS or limit <= 12:
                    break
                print(f"  {len(table['ra']):,} rows at G < {limit}; lowering the limit", flush=True)
                limit -= 1
        except Exception as error:
            errors.append(f"{source}: {type(error).__name__}: {str(error)[:200]}")
            print(f"  {source} failed; {'trying next source' if source != sources[-1] else 'giving up'}")
            continue
        path = out / f"{position.name}.csv"
        write_csv(path, table)
        return {
            "file": path.name,
            "ra": ra_text,
            "dec": dec_text,
            "fov": fov,
            "maxmag": limit,
            "rows": len(table["ra"]),
            "date": today,
            "source": source.upper(),
        }
    raise RuntimeError("; ".join(errors))


def write_csv(path: Path, table: dict[str, np.ndarray]) -> None:
    columns = list(sc.REQUIRED)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row in zip(*(table[c] for c in columns), strict=True):
            writer.writerow(["" if np.isnan(v) else repr(float(v)) for v in row])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=ROOT / "build" / "snapshots")
    parser.add_argument("--source", choices=("auto", "esa", "gavo"), default="auto")
    parser.add_argument("--only", nargs="+", metavar="NAME", help="subset of clusters, e.g. M13 M44")
    parser.add_argument("--refresh", action="store_true", help="query again even if a snapshot exists")
    parser.add_argument("--dry-run", action="store_true", help="list what would be queried")
    args = parser.parse_args(argv)
    names = [n.upper() for n in args.only] if args.only else list(PLAN)
    unknown = [n for n in names if n not in PLAN]
    if unknown:
        raise SystemExit(f"Not in the snapshot plan: {', '.join(unknown)}")
    sources = ["esa", "gavo"] if args.source == "auto" else [args.source]
    if args.dry_run:
        for name in names:
            print(f"{name}: FOV {PLAN[name][0]}°, G < {PLAN[name][1]}, from {' then '.join(sources)}")
        return 0

    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    index_path = out / "index.json"
    index: dict[str, dict] = json.loads(index_path.read_text()) if index_path.exists() else {}
    today = datetime.now(UTC).date().isoformat()
    failed: list[str] = []
    for name in names:
        entry = index.get(name)
        if entry and not args.refresh and (out / entry["file"]).is_file():
            print(f"{name}: kept ({entry['rows']:,} rows from {entry.get('source', 'ESA')}, {entry['date']})")
            continue
        print(f"{name}: querying {' then '.join(sources)}…", flush=True)
        try:
            index[name] = build_one(name, out, today, sources)
        except Exception as error:
            print(f"{name}: FAILED ({error})", flush=True)
            failed.append(name)
            continue
        index_path.write_text(json.dumps(index, indent=1, sort_keys=True) + "\n")  # resume point
        print(f"{name}: {index[name]['rows']:,} rows from {index[name]['source']}", flush=True)
    total = sum((out / e["file"]).stat().st_size for e in index.values())
    print(f"{len(index)} snapshots, {total / 1e6:.1f} MB in {out}")
    if failed:
        print(f"Failed: {', '.join(failed)}. Re-run to retry only these.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
