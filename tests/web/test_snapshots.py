"""Tests for scripts/build_snapshots.py and build_site's snapshot validation (no network)."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "web"))

import build_site  # noqa: E402
import build_snapshots as bs  # noqa: E402
import showclusters as sc  # noqa: E402

ESA_CSV = (FIXTURES / "m13_esa.csv").read_bytes()
GAVO_CSV = (FIXTURES / "m13_gavo.csv").read_bytes()


def limited(data: bytes, scratch: Path):  # type: ignore[no-untyped-def]
    """A fake service that honours the query's magnitude limit, like the real one."""
    table = sc.parse_gaia_table(data, "x.csv")

    def fetch(adql: str) -> bytes:
        limit = float(re.search(r"phot_g_mean_mag < +([0-9.]+)", adql).group(1))  # type: ignore[union-attr]
        keep = table["phot_g_mean_mag"] < limit
        path = scratch / f"limit_{limit}.csv"
        bs.write_csv(path, {k: v[keep] for k, v in table.items()})
        return path.read_bytes()

    return fetch


def failing(adql: str) -> bytes:
    raise ConnectionError("Server overloaded")


def no_sleep(seconds: float) -> None:
    pass


def test_esa_snapshot_matches_source_and_index(tmp_path: Path) -> None:
    fetchers = {"esa": ("esa", lambda q: ESA_CSV), "gavo": ("gavo", failing)}
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa", "gavo"], fetchers, no_sleep)
    assert entry == {
        "file": "M13.csv",
        "ra": "+16 : 41 : 42",
        "dec": "+36 : 27 : 39",
        "fov": 0.5,
        "maxmag": 18.0,
        "rows": 1498,
        "date": "2026-11-02",
        "source": "ESA",
    }
    saved = sc.parse_gaia_table((tmp_path / "M13.csv").read_bytes(), "M13.csv")
    original = sc.parse_gaia_table(ESA_CSV, "x.csv")
    for name in sc.REQUIRED:
        np.testing.assert_array_equal(saved[name], original[name], err_msg=name)


def test_falls_back_to_gavo_when_esa_fails(tmp_path: Path) -> None:
    calls: list[str] = []

    def esa(adql: str) -> bytes:
        calls.append("esa")
        raise ConnectionError("503")

    fetchers = {"esa": ("esa", esa), "gavo": ("gavo", lambda q: GAVO_CSV)}
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa", "gavo"], fetchers, no_sleep)
    assert entry["source"] == "GAVO" and entry["rows"] == 1498
    assert calls == ["esa"] * bs.ATTEMPTS  # retried before giving up on ESA
    assert "gaia.dr3lite" not in (tmp_path / "M13.csv").read_text()


def test_all_sources_failing_raises(tmp_path: Path) -> None:
    fetchers = {"esa": ("esa", failing), "gavo": ("gavo", failing)}
    with pytest.raises(RuntimeError, match="esa: ConnectionError.*gavo: ConnectionError"):
        bs.build_one("M13", tmp_path, "2026-11-02", ["esa", "gavo"], fetchers, no_sleep)
    assert not (tmp_path / "M13.csv").exists()


def test_magnitude_limit_drops_until_rows_fit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bs, "MAX_SNAPSHOT_ROWS", 600)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    fetchers = {"esa": ("esa", limited(ESA_CSV, scratch)), "gavo": ("gavo", failing)}
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa"], fetchers, no_sleep)
    assert entry["maxmag"] < 16 and 0 < entry["rows"] <= 600
    saved = sc.parse_gaia_table((tmp_path / "M13.csv").read_bytes(), "M13.csv")
    assert saved["phot_g_mean_mag"].max() < entry["maxmag"]


def test_main_resumes_and_reports_failures(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bs, "FETCHERS", {"esa": ("esa", lambda q: ESA_CSV), "gavo": ("gavo", failing)})
    monkeypatch.setattr(bs.time, "sleep", no_sleep)
    original_build = bs.build_one

    def build(name, out, today, sources):  # type: ignore[no-untyped-def]
        if name == "M44":
            raise RuntimeError("esa: 500")
        return original_build(name, out, today, sources, bs.FETCHERS, no_sleep)

    monkeypatch.setattr(bs, "build_one", build)
    assert bs.main(["--out", str(tmp_path), "--only", "M13", "M44", "--source", "esa"]) == 1
    index = json.loads((tmp_path / "index.json").read_text())
    assert list(index) == ["M13"]
    stamp = (tmp_path / "M13.csv").stat().st_mtime_ns
    assert bs.main(["--out", str(tmp_path), "--only", "M13"]) == 0  # kept, not re-queried
    assert (tmp_path / "M13.csv").stat().st_mtime_ns == stamp


def test_unknown_cluster_is_rejected() -> None:
    with pytest.raises(SystemExit, match="Not in the snapshot plan: M999"):
        bs.main(["--only", "M999", "--dry-run"])


def test_every_planned_cluster_is_in_the_catalogs() -> None:
    for name in bs.PLAN:
        assert sc.lookup_object(name, ROOT / "csv").name == name
    assert len(bs.PLAN) > 70


def snapshot_dir(tmp_path: Path) -> Path:
    fetchers = {"esa": ("esa", lambda q: ESA_CSV), "gavo": ("gavo", failing)}
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa"], fetchers, no_sleep)
    (tmp_path / "index.json").write_text(json.dumps({"M13": entry}))
    return tmp_path


def test_site_build_accepts_a_valid_snapshot_directory(tmp_path: Path) -> None:
    assert build_site.check_snapshots(snapshot_dir(tmp_path))["M13"]["rows"] == 1498


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda e: e.update(rows=5), "1498 rows, index says 5"),
        (lambda e: e.update(file="../../etc/passwd"), "bad or missing file"),
        (lambda e: e.pop("date"), "lacks \\['date'\\]"),
    ],
)
def test_site_build_rejects_bad_snapshots(tmp_path: Path, change, message: str) -> None:  # type: ignore[no-untyped-def]
    directory = snapshot_dir(tmp_path)
    index = json.loads((directory / "index.json").read_text())
    change(index["M13"])
    (directory / "index.json").write_text(json.dumps(index))
    with pytest.raises(SystemExit, match=message):
        build_site.check_snapshots(directory)


def test_site_build_rejects_empty_index(tmp_path: Path) -> None:
    (tmp_path / "index.json").write_text("{}")
    with pytest.raises(SystemExit, match="empty"):
        build_site.check_snapshots(tmp_path)
