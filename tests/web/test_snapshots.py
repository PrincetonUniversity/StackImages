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
sys.path.insert(0, str(Path(__file__).parent))

import build_site  # noqa: E402
import build_snapshots as bs  # noqa: E402
import showclusters as sc  # noqa: E402
from fakegaia import FakeGaia  # noqa: E402

ESA_CSV = (FIXTURES / "m13_esa.csv").read_bytes()
M13_FOV = 0.2  # the recorded fixtures cover M13 out to this field of view


@pytest.fixture(autouse=True)
def small_m13_plan(monkeypatch: pytest.MonkeyPatch) -> None:
    """Snapshot M13 at the size the fixtures cover, so fakes can answer it completely."""
    monkeypatch.setitem(bs.PLAN, "M13", (M13_FOV, 16.0))


# The snapshot's field: M13's displayed position, M13_FOV, G < 16.
M13_FIELD = (sc.parse_sexagesimal("+16 : 41 : 42"), sc.parse_sexagesimal("+36 : 27 : 39"), M13_FOV, 16.0)


def m13_query_rows() -> int:
    """Stars the snapshot's own (unsplit) M13 query selects from the fixture."""
    body = FakeGaia().answer(sc.gaia_adql(*M13_FIELD, "gavo"))[1]
    return len(sc.parse_gaia_table(body, "x.csv")["ra"])


def no_sleep(seconds: float) -> None:
    pass


def service(fake: FakeGaia):  # type: ignore[no-untyped-def]
    def fetch(adql: str) -> bytes:
        status, body = fake.answer(adql)
        if status == 500 and b"timed out" in body:
            raise bs.GavoTimeoutError("GAVO query timed out")
        if status != 200:
            raise ConnectionError(f"HTTP {status}")
        return body

    return fetch


def failing(adql: str) -> bytes:
    raise ConnectionError("Server overloaded")


def fields(esa=failing, gavo=failing):  # type: ignore[no-untyped-def]
    return {"esa": bs.esa_field(esa, no_sleep), "gavo": bs.gavo_field(gavo, no_sleep)}


def test_esa_snapshot_matches_source_and_index(tmp_path: Path) -> None:
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa", "gavo"], fields(esa=lambda q: ESA_CSV))
    assert entry == {
        "file": "M13.csv",
        "ra": "+16 : 41 : 42",
        "dec": "+36 : 27 : 39",
        "fov": M13_FOV,
        "maxmag": 16.0,
        "rows": m13_query_rows(),
        "date": "2026-11-02",
        "source": "ESA",
    }
    saved = sc.parse_gaia_table((tmp_path / "M13.csv").read_bytes(), "M13.csv")
    expected = sc.assemble([sc.parse_gaia_table(ESA_CSV, "x.csv")], *M13_FIELD)
    for name in sc.REQUIRED:
        np.testing.assert_array_equal(saved[name], expected[name], err_msg=name)


def test_falls_back_to_gavo_when_esa_fails(tmp_path: Path) -> None:
    calls: list[str] = []

    def esa(adql: str) -> bytes:
        calls.append(adql)
        raise ConnectionError("503")

    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa", "gavo"], fields(esa, service(FakeGaia())))
    assert entry["source"] == "GAVO" and entry["rows"] == m13_query_rows() > 1400
    assert len(calls) == bs.ESA_ATTEMPTS  # retried, then moved on


def test_gavo_timeouts_are_split_not_retried(tmp_path: Path) -> None:
    fake = FakeGaia(timeout_radius=0.11)  # the single 0.141° query times out; its 0.102° quarters don't
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["gavo"], fields(gavo=service(fake)))
    assert len(fake.queries) == 1 + 4  # one timeout, then four pieces; nothing retried
    split = sc.parse_gaia_table((tmp_path / "M13.csv").read_bytes(), "M13.csv")
    whole = sc.assemble([sc.parse_gaia_table(FakeGaia().answer(fake.queries[0])[1], "x.csv")], *M13_FIELD)
    assert entry["rows"] == len(whole["ra"]) == m13_query_rows()
    for name in sc.REQUIRED:  # star for star, same values, same order as one query
        np.testing.assert_array_equal(split[name], whole[name], err_msg=name)


def test_field_too_dense_even_when_split_lowers_the_limit(tmp_path: Path) -> None:
    attempts: list[float] = []

    def dense(adql: str) -> bytes:
        attempts.append(float(re.search(r"phot_g_mean_mag < +([0-9.]+)", adql).group(1)))  # type: ignore[union-attr]
        if attempts[-1] > 12.5:
            raise bs.GavoTimeoutError("timed out")
        return FakeGaia().answer(adql)[1]

    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["gavo"], fields(gavo=dense))
    assert entry["maxmag"] == 12.0 and entry["source"] == "GAVO"
    assert sorted(set(attempts), reverse=True) == [16.0, 14.0, 12.0]


def test_all_sources_failing_raises(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match=r"esa \(G < 16.0\): ConnectionError.*gavo \(G < 12.0\)"):
        bs.build_one("M13", tmp_path, "2026-11-02", ["esa", "gavo"], fields())
    assert not (tmp_path / "M13.csv").exists()


def test_magnitude_limit_drops_until_rows_fit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bs, "MAX_SNAPSHOT_ROWS", 600)
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["gavo"], fields(gavo=service(FakeGaia())))
    assert entry["maxmag"] < 16 and 0 < entry["rows"] <= 600
    saved = sc.parse_gaia_table((tmp_path / "M13.csv").read_bytes(), "M13.csv")
    assert saved["phot_g_mean_mag"].max() < entry["maxmag"]


def test_main_resumes_and_reports_failures(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(bs.FIELDS, "esa", lambda: bs.esa_field(lambda q: ESA_CSV, no_sleep))
    original_build = bs.build_one

    def build(name, out, today, sources):  # type: ignore[no-untyped-def]
        if name == "M44":
            raise RuntimeError("esa: 500")
        return original_build(name, out, today, sources)

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
    entry = bs.build_one("M13", tmp_path, "2026-11-02", ["esa"], fields(esa=lambda q: ESA_CSV))
    (tmp_path / "index.json").write_text(json.dumps({"M13": entry}))
    return tmp_path


def test_site_build_accepts_a_valid_snapshot_directory(tmp_path: Path) -> None:
    assert build_site.check_snapshots(snapshot_dir(tmp_path))["M13"]["rows"] == m13_query_rows()


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda e: e.update(rows=5), "rows, index says 5"),
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
