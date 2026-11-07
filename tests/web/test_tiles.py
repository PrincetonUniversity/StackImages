"""Large fields are fetched in pieces; the end result must equal one query's, exactly."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(ROOT / "web"))
sys.path.insert(0, str(Path(__file__).parent))

import showclusters as sc  # noqa: E402
from fakegaia import FakeGaia  # noqa: E402


def angle(ra1: float, dec1: float, ra2: Any, dec2: Any) -> Any:
    a, b = np.radians(dec1), np.radians(dec2)
    h = np.sin((b - a) / 2) ** 2 + np.cos(a) * np.cos(b) * np.sin(np.radians(ra2 - ra1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1))))


def centre(tile: dict[str, Any], ra_hours: float, dec: float) -> tuple[float, float]:
    return sc._sky(ra_hours * 15, dec, tile["xi"], tile["eta"])


def random_points_in_circle(
    ra_deg: float, dec: float, radius: float, n: int
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(7)
    r = radius * np.sqrt(rng.random(n)) * 0.9999
    theta = rng.random(n) * 2 * np.pi
    points = [sc._sky(ra_deg, dec, x, y) for x, y in zip(r * np.cos(theta), r * np.sin(theta), strict=True)]
    ras, decs = zip(*points, strict=True)
    return np.array(ras), np.array(decs)


# --- geometry ---------------------------------------------------------------------------


def test_small_field_is_one_unchanged_query() -> None:
    tiles = sc.field_tiles(16.695, 36.46, 0.3, 16, "gavo")
    assert len(tiles) == 1
    assert tiles[0]["query"] == sc.gaia_adql(16.695, 36.46, 0.3, 16, "gavo")


@pytest.mark.parametrize(
    ("ra_hours", "dec", "fov"),
    [(8.67, 19.67, 2.0), (8.67, 19.67, 5.0), (0.02, -5.0, 3.0), (2.5, 87.0, 4.0), (12.0, -88.5, 2.5)],
)
def test_pieces_cover_the_whole_field(ra_hours: float, dec: float, fov: float) -> None:
    tiles = sc.field_tiles(ra_hours, dec, fov, 16, "gavo")
    assert len(tiles) > 1 and max(t["radius"] for t in tiles) <= sc.TILE_RADIUS + 1e-9
    ras, decs = random_points_in_circle(ra_hours * 15, dec, fov / math.sqrt(2), 4000)
    covered = np.zeros(len(ras), dtype=bool)
    for tile in tiles:
        ra_c, dec_c = centre(tile, ra_hours, dec)
        covered |= angle(ra_c, dec_c, ras, decs) < tile["radius"]
    assert covered.all(), f"{(~covered).sum()} uncovered points"


def in_square(tile: dict[str, Any], xi: np.ndarray, eta: np.ndarray) -> np.ndarray:
    return (np.abs(xi - tile["xi"]) <= tile["half"]) & (np.abs(eta - tile["eta"]) <= tile["half"])


def plane_points(radius: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Random tangent-plane points (degrees) inside the projection of a sky circle."""
    rng = np.random.default_rng(11)
    r = math.degrees(math.tan(math.radians(radius))) * np.sqrt(rng.random(n)) * 0.9999
    theta = rng.random(n) * 2 * np.pi
    return r * np.cos(theta), r * np.sin(theta)


@pytest.mark.parametrize("fov", [2.0, 5.0])
def test_piece_squares_cover_the_field(fov: float) -> None:
    """What splitting relies on: the pieces' squares (not just their circles) cover the field."""
    tiles = sc.field_tiles(8.67, 19.67, fov, 16, "gavo")
    xi, eta = plane_points(fov / math.sqrt(2), 4000)
    covered = np.zeros(len(xi), dtype=bool)
    for tile in tiles:
        covered |= in_square(tile, xi, eta)
    assert covered.all(), f"{(~covered).sum()} uncovered points"


@pytest.mark.parametrize("fov", [0.2, 3.0])
def test_split_pieces_cover_their_parent(fov: float) -> None:
    ra_hours, dec = 8.67, 19.67
    parent = sc.field_tiles(ra_hours, dec, fov, 16, "gavo")[0]
    children = sc.split_tile(parent, ra_hours, dec, 16, "gavo")
    assert len(children) == 4 and all(c["radius"] < parent["radius"] for c in children)
    rng = np.random.default_rng(3)
    xi = parent["xi"] + (rng.random(3000) * 2 - 1) * parent["half"]
    eta = parent["eta"] + (rng.random(3000) * 2 - 1) * parent["half"]
    covered = np.zeros(len(xi), dtype=bool)
    for child in children:
        covered |= in_square(child, xi, eta)
    assert covered.all()
    # and each child's query circle really contains its square (sky distances)
    for child in children:
        ra_c, dec_c = centre(child, ra_hours, dec)
        corners = [
            sc._sky(ra_hours * 15, dec, child["xi"] + dx * child["half"], child["eta"] + dy * child["half"])
            for dx in (-1, 1)
            for dy in (-1, 1)
        ]
        assert max(angle(ra_c, dec_c, r, d) for r, d in corners) < child["radius"]


def test_splitting_stops_at_the_minimum_size() -> None:
    tile = sc.field_tiles(16.695, 36.46, 0.2, 16, "gavo")[0]
    depth = 0
    while smaller := sc.split_tile(tile, 16.695, 36.46, 16, "gavo"):
        tile, depth = smaller[0], depth + 1
    assert tile["radius"] >= sc.MIN_TILE_RADIUS and depth >= 1


ORION = (5.6036, -1.2019)  # Alnilam, the middle of Orion's belt (HW: FOV 60, magnitude 7)


@pytest.mark.parametrize(("fov", "maxmag"), [(60.0, 7.0), (60.0, sc.WIDE_MAXMAG), (90.0, 6.0), (8.0, 8.0)])
def test_wide_bright_field_is_one_magnitude_index_query(fov: float, maxmag: float) -> None:
    """GAVO answers 60° at G < 8.5 in ~2 s through its magnitude index (measured 2026-11-07);
    pieces would need thousands of queries."""
    tiles = sc.field_tiles(*ORION, fov, maxmag, "gavo")
    assert len(tiles) == 1
    query = tiles[0]["query"]
    assert "+ 0" not in query and f"< {maxmag:5f}" in query
    assert tiles[0]["radius"] == pytest.approx(fov / math.sqrt(2))
    # Same circle as the desktop's query; splitting it would only make things slower.
    assert sc.split_tile(tiles[0], *ORION, maxmag, "gavo") == []


@pytest.mark.parametrize(("fov", "maxmag"), [(45.0, 10.0), (30.0, 12.0), (60.0, 9.0), (20.0, 16.0)])
def test_wide_deep_field_is_one_async_job(session: sc.Session, fov: float, maxmag: float) -> None:
    """Too many pieces to fetch in time and too deep for the magnitude index's quick answer: the
    page sends the whole field as one GAVO async job instead (measured 2026-11-07: FOV 45/G 10 in
    112 s, FOV 30/G 12 in 279 s, 120,593 stars; ESA's anonymous downloads stop at 50,000)."""
    assert sc.field_tiles(*ORION, fov, maxmag, "gavo") == []
    ra, dec = (sc.format_sexagesimal(v) for v in ORION)
    plan = session.plan(ra, dec, fov, maxmag)
    assert plan["gavo_tiles"] == []
    query = plan["gavo_async"]
    assert "gaia.dr3lite" in query and "+ 0" not in query and f"< {maxmag:5f}" in query
    assert f"{fov / math.sqrt(2):14f}".strip() in query


@pytest.mark.parametrize(("fov", "maxmag"), [(0.3, 16.0), (5.0, 16.0), (60.0, 7.0)])
def test_fields_gavo_answers_quickly_need_no_async_job(
    session: sc.Session, fov: float, maxmag: float
) -> None:
    ra, dec = (sc.format_sexagesimal(v) for v in ORION)
    plan = session.plan(ra, dec, fov, maxmag)
    assert plan["gavo_tiles"] and plan["gavo_async"] is None


def test_largest_field_that_is_still_fetched_in_pieces() -> None:
    tiles = sc.field_tiles(*ORION, 5.0, 16, "gavo")
    assert 1 < len(tiles) <= sc.MAX_TILES
    assert all("+ 0" in t["query"] for t in tiles)


def test_wide_bright_field_esa_query_is_the_desktop_query() -> None:
    query = sc.gaia_adql(*ORION, 60, 7, "esa")
    assert "gaiadr3.gaia_source" in query and "42.426407" in query and "< 7.000000" in query


def test_merge_removes_overlap_duplicates_only() -> None:
    table = sc.parse_gaia_table((Path(__file__).parent / "fixtures" / "m13_gavo.csv").read_bytes(), "x.csv")
    first = {k: v[:900] for k, v in table.items()}
    second = {k: v[600:] for k, v in table.items()}  # 300 rows in both
    merged = sc.merge_tables([first, second])
    assert len(merged["ra"]) == len(table["ra"])
    np.testing.assert_array_equal(np.sort(merged["ra"]), np.sort(table["ra"]))


# --- the whole pipeline: one query vs. pieces -------------------------------------------


RA_TEXT, DEC_TEXT, FOV, MAXMAG = "+16 : 41 : 42", "+36 : 27 : 39", 0.2, 16.0  # M13, the fixtures' field


def fetch_like_the_page(session: sc.Session, fake: FakeGaia) -> list[bytes]:
    """What app.js does: fetch each piece; a piece that times out is split and re-queued."""
    queue = list(session.plan(RA_TEXT, DEC_TEXT, FOV, MAXMAG)["gavo_tiles"])
    parts = []
    while queue:
        tile = queue.pop(0)
        status, body = fake.answer(tile["query"])
        if status == 200:
            parts.append(body)
            continue
        smaller = session.split(tile)
        assert smaller, "field too dense for the fake service"
        queue.extend(smaller)
    return parts


def run(session: sc.Session, parts: list[bytes]) -> dict[str, Any]:
    result = session.show_stars(parts, "gavo.csv", "test", RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13")
    assert session.frame is not None
    return {
        "result": result,
        "stars": session.frame.stars,
        "hr": {m: sc.hr_arrays(session.frame.stars, m) for m in ("Apparent", "Absolute")},
        "png": session.hr_diagram("Absolute"),
    }


@pytest.mark.parametrize(("timeout_radius", "pieces"), [(0.11, 4), (0.06, 16)])
def test_pieces_give_exactly_the_one_query_result(
    session: sc.Session, timeout_radius: float, pieces: int
) -> None:
    whole_fake, split_fake = FakeGaia(), FakeGaia(timeout_radius=timeout_radius)
    whole = run(session, fetch_like_the_page(session, whole_fake))
    split_parts = fetch_like_the_page(session, split_fake)
    assert len(whole_fake.queries) == 1 and len(split_parts) == pieces
    split = run(session, split_parts)

    for key in ("count", "brightest", "pct_lo", "pct_hi", "size"):
        assert split["result"][key] == whole["result"][key], key
    assert split["result"]["count"] > 1000
    assert split["result"]["image"] == whole["result"]["image"]  # the picture, byte for byte
    for name, values in whole["stars"].__dict__.items():
        np.testing.assert_array_equal(getattr(split["stars"], name), values, err_msg=name)
    for mode, arrays in whole["hr"].items():
        for key, values in arrays.items():
            np.testing.assert_array_equal(split["hr"][mode][key], values, err_msg=f"{mode} {key}")
    assert split["png"] == whole["png"]  # the downloaded HR diagram, byte for byte


def test_dense_field_is_one_quick_query(session: sc.Session) -> None:
    """A service that times out on a bare magnitude condition (GAVO in M23's field) must
    still answer the page's query at once, without splitting."""
    fake = FakeGaia(slow_magnitude_index=True)
    parts = fetch_like_the_page(session, fake)
    assert len(parts) == 1 and len(fake.queries) == 1
    whole = run(session, fetch_like_the_page(session, FakeGaia()))
    assert run(session, parts)["result"]["count"] == whole["result"]["count"] > 1000


@pytest.mark.parametrize("uploaded", [True, False])
def test_outside_note_only_for_uploaded_files(session: sc.Session, uploaded: bool) -> None:
    data = (Path(__file__).parent / "fixtures" / "m13_esa.csv").read_bytes()
    # A field centred 22 s of RA away from the file's centre: only ~40% of its stars are inside.
    result = session.show_stars(
        data, "m13.csv", "test", "+16 : 41 : 20", DEC_TEXT, 0.1, MAXMAG, "M13", 0, uploaded
    )
    notes = result["notes"]
    if uploaded:
        assert notes == ["Only 604 of 1,498 stars in the file are inside the field."]
    else:
        assert notes == []


ORION_TEXT = ("+05 : 36 : 13", "-01 : 12 : 07")


def test_orion_at_60_degrees_shows_stars_with_absolute_magnitudes(session: sc.Session) -> None:
    """Homework: Orion, FOV 60, magnitude 7 (recorded from GAVO's magnitude-index query)."""
    data = (FIXTURES / "orion60_gavo.csv").read_bytes()
    result = session.show_stars([data], "gavo.csv", "test", *ORION_TEXT, 60.0, 7.0, "Orion")
    assert result["rows"] == 3259
    assert 2000 < result["count"] <= 3259 + 200  # plus Hipparcos stars brighter than Gaia's limit
    assert result["magnitude_type"] == "Absolute"  # ass.py above 5°
    assert result["pct_hi"] > 0
    assert session.hr_diagram("Absolute")[:8] == b"\x89PNG\r\n\x1a\n"


def test_orion_gavo_and_esa_give_the_same_stars() -> None:
    gavo = sc.parse_gaia_table((FIXTURES / "orion60_gavo.csv").read_bytes(), "gavo.csv")
    esa = sc.parse_gaia_table((FIXTURES / "orion60_esa.vot").read_bytes(), "esa.vot")
    order_g, order_e = np.lexsort((gavo["dec"], gavo["ra"])), np.lexsort((esa["dec"], esa["ra"]))
    for name in ("ra", "dec", "phot_g_mean_mag", "parallax"):
        np.testing.assert_allclose(gavo[name][order_g], esa[name][order_e], rtol=1e-6, err_msg=name)


def test_narrow_field_starts_with_apparent_magnitudes(session: sc.Session) -> None:
    result = session.show_stars(
        [(FIXTURES / "m13_gavo.csv").read_bytes()], "gavo.csv", "test", RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13"
    )
    assert result["magnitude_type"] == "Apparent"


# --- ESA downloads that stopped early ---------------------------------------------------


def m13_rows(n: int) -> bytes:
    """An ESA-style CSV of n stars inside the M13 field (the fixture's rows, repeated)."""
    lines = (FIXTURES / "m13_esa.csv").read_text().splitlines()
    body = [lines[1 + i % (len(lines) - 1)] for i in range(n)]
    return ("\n".join([lines[0], *body]) + "\n").encode()


def test_esa_download_cut_at_50000_stars_is_refused(session: sc.Session) -> None:
    """ESA stops anonymous results at 50,000 rows without a sign in CSV (measured 2026-11-07:
    Orion FOV 30/G 12 came back as exactly 50,000 of ~120,000)."""
    with pytest.raises(sc.InputError, match=r"50,000 stars.*incomplete"):
        session.show_stars(
            m13_rows(50_000), "esa.csv", "test", RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13", 0, True
        )


def test_overflow_flag_in_votable_is_refused(session: sc.Session) -> None:
    data = (FIXTURES / "m13_esa.vot").read_bytes()
    marker = b'<RESOURCE type="results">'
    assert marker in data
    flagged = data.replace(marker, marker + b'\n<INFO name="QUERY_STATUS" value="OVERFLOW"/>', 1)
    with pytest.raises(sc.InputError, match="incomplete"):
        session.show_stars(flagged, "esa.vot", "test", RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13", 0, True)


@pytest.mark.parametrize("n", [49_999, 50_001])
def test_other_large_uploads_are_accepted(session: sc.Session, n: int) -> None:
    result = session.show_stars(
        m13_rows(n), "esa.csv", "test", RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13", 0, True
    )
    assert result["rows"] == n


def test_gavo_results_of_50000_rows_are_not_refused(session: sc.Session) -> None:
    """Only a student's ESA file is checked; a GAVO piece can hold exactly 50,000 stars."""
    result = session.show_stars([m13_rows(50_000)], "gavo.csv", "test", RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13")
    assert result["rows"] == 50_000
