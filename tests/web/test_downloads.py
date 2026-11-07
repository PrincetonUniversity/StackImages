"""The CSV downloads: star data (re-uploadable) and HR diagram points, with credit."""

from __future__ import annotations

import csv
import io
from pathlib import Path

import numpy as np
import pytest

import showclusters as sc
from fakegaia import FakeGaia

FIXTURES = Path(__file__).parent / "fixtures"
RA_TEXT, DEC_TEXT, FOV, MAXMAG = "+16 : 41 : 42", "+36 : 27 : 39", 0.2, 16.0
GAVO = "Gaia DR3 via GAVO"


def show(session: sc.Session, data: object, source: str = GAVO, uploaded: bool = False) -> dict:
    return session.show_stars(data, "x.csv", source, RA_TEXT, DEC_TEXT, FOV, MAXMAG, "M13", 0, uploaded)


def gavo_bytes(timeout_radius: float | None = None) -> list[bytes]:
    """The page's fetch: one query, or pieces when the fake times out."""
    fake = FakeGaia(timeout_radius=timeout_radius)
    probe = sc.Session(Path("."), np.empty((0, 7)), np.zeros((2, 2, 3), np.float32))
    queue = list(probe.plan(RA_TEXT, DEC_TEXT, FOV, MAXMAG)["gavo_tiles"])
    parts = []
    while queue:
        tile = queue.pop(0)
        status, body = fake.answer(tile["query"])
        if status == 200:
            parts.append(body)
        else:
            queue.extend(probe.split(tile))
    return parts


def comments(data: bytes) -> list[str]:
    return [line[2:] for line in data.decode().splitlines() if line.startswith("# ")]


def table(data: bytes) -> list[dict[str, str]]:
    lines = [line for line in data.decode().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines))))


# --- star data ----------------------------------------------------------------------------


def test_star_data_reuploads_to_the_identical_result(session: sc.Session) -> None:
    first = show(session, gavo_bytes())
    used = {k: v.copy() for k, v in session.gaia.items()}  # type: ignore[union-attr]
    first_hr = {m: sc.hr_arrays(session.frame.stars, m) for m in ("Apparent", "Absolute")}  # type: ignore[union-attr]
    download = session.star_data("2026-11-03")
    assert download["name"] == "M13_FOV0.2_G16_GAVO_2026-11-03_stars.csv"

    parsed = sc.parse_gaia_table(download["data"], download["name"])
    for name in sc.REQUIRED:  # exactly the rows used, value for value
        np.testing.assert_array_equal(parsed[name], used[name], err_msg=name)

    again = show(session, download["data"], "Gaia DR3 from uploaded ESA file", uploaded=True)
    for key in ("count", "brightest", "image", "pct_lo", "pct_hi"):
        assert again[key] == first[key], key
    assert again["notes"] == []  # the file is exactly the field
    assert session.frame is not None
    for mode, points in first_hr.items():  # the HR diagram's points, identical
        for key, values in points.items():
            np.testing.assert_array_equal(sc.hr_arrays(session.frame.stars, mode)[key], values, err_msg=key)


def test_star_data_is_the_same_whether_fetched_whole_or_in_pieces(session: sc.Session) -> None:
    show(session, gavo_bytes())
    whole = session.star_data("2026-11-03")["data"]
    show(session, gavo_bytes(timeout_radius=0.06))  # 16 pieces
    assert session.star_data("2026-11-03")["data"] == whole


def test_star_data_credit_and_columns(session: sc.Session) -> None:
    show(session, gavo_bytes())
    data = session.star_data("2026-11-03")["data"]
    text = " ".join(comments(data))
    assert "ESA/Gaia/DPAC" in text and "CC BY-NC 3.0 IGO" in text
    assert "This work has made use of data from the European Space Agency (ESA) mission Gaia" in text
    assert "Gaia Multilateral Agreement" in text
    assert "Prusti et al. 2016, A&A 595, A1" in text and "Vallenari et al. 2023, A&A 674, A1" in text
    assert "GAVO" in text and "M13, RA +16 : 41 : 42, Dec +36 : 27 : 39, FOV 0.2 deg, G < 16" in text
    assert list(table(data)[0]) == list(sc.REQUIRED)
    assert "HYG" not in text  # the star data holds Gaia rows only


def test_esa_sourced_star_data_has_no_gavo_line(session: sc.Session) -> None:
    show(session, (FIXTURES / "m13_esa.csv").read_bytes(), "Gaia DR3 saved copy from ESA (2026-11-02)")
    download = session.star_data("2026-11-03")
    assert download["name"].endswith("_saved-ESA_2026-11-03_stars.csv")
    assert not any("GAVO" in line for line in comments(download["data"])[2:])


# --- HR points ----------------------------------------------------------------------------


@pytest.mark.parametrize("magnitude_type", ["Apparent", "Absolute"])
def test_hr_points_match_the_diagram(session: sc.Session, magnitude_type: str) -> None:
    show(session, gavo_bytes())
    stars = session.frame.stars  # type: ignore[union-attr]
    download = session.hr_points(magnitude_type, "2026-11-03")
    assert download["name"] == f"M13_HR_{magnitude_type}_points.csv"
    rows = table(download["data"])
    assert len(rows) == len(stars.mag) > 1000

    drawn = {(np.float32(r["x_plotted"]), np.float32(r["y_plotted"])) for r in rows if r["drawn"] == "true"}
    plotted = sc.hr_arrays(stars, magnitude_type)
    shown = {(x, y) for x, y, s in zip(plotted["x"], plotted["y"], plotted["s"], strict=True) if s > 0}
    assert drawn == shown  # the CSV's drawn points are exactly the diagram's points

    hyg = [r for r in rows if r["catalog"] == "HYG"]
    gaia = [r for r in rows if r["catalog"] == "Gaia DR3"]
    assert len(hyg) == int(stars.hyg.sum()) and len(gaia) == len(rows) - len(hyg)
    assert all(r["parallax_mas"] == "" for r in hyg) and all(r["parallax_mas"] for r in gaia)
    assert rows[0]["drawn"] == "false"
    for r in rows[:50]:
        assert float(r["x_plotted"]) == pytest.approx(1.25 * (float(r["color_index"]) - 1), abs=1e-5)


def test_hr_points_credit_includes_hyg(session: sc.Session) -> None:
    show(session, gavo_bytes())
    text = " ".join(comments(session.hr_points("Apparent", "2026-11-03")["data"]))
    assert "ESA/Gaia/DPAC" in text and "CC BY-NC 3.0 IGO" in text
    assert "HYG database v3" in text and "CC BY-SA 2.5" in text
    assert "drawn = false" in text


def test_hr_points_include_bright_hyg_stars_when_in_the_picture(session: sc.Session) -> None:
    pleiades = session.lookup("m45")
    empty = b"ra,dec,parallax,phot_g_mean_mag,bp_rp,pmra,pmdec\n"
    session.show_stars(empty, "x.csv", GAVO, pleiades["ra"], pleiades["dec"], 2.0, 16.0, "M45", 0)
    rows = table(session.hr_points("Apparent", "2026-11-03")["data"])
    assert rows and all(r["catalog"] == "HYG" for r in rows)  # Alcyone and friends


def test_downloads_need_stars_first() -> None:
    fresh = sc.Session(Path("."), np.empty((0, 7)), np.zeros((2, 2, 3), np.float32))
    with pytest.raises(sc.InputError, match="Show Stars first"):
        fresh.star_data("2026-11-03")
    with pytest.raises(sc.InputError, match="Show Stars first"):
        fresh.hr_points("Apparent", "2026-11-03")


# --- reading files with comment lines ------------------------------------------------------


def test_comment_lines_are_skipped_and_line_numbers_stay_true() -> None:
    data = b"# a note\n# another\nra,dec,parallax,phot_g_mean_mag,bp_rp,pmra,pmdec\n1,2,3,4,5,6,7\n1,2\n"
    with pytest.raises(sc.InputError, match="Line 5 of the CSV has 2 fields"):
        sc.parse_gaia_table(data, "x.csv")
    ok = sc.parse_gaia_table(data.rsplit(b"1,2\n", 1)[0], "x.csv")
    assert len(ok["ra"]) == 1


def test_custom_fields_get_a_generic_file_name(session: sc.Session) -> None:
    session.show_stars(
        gavo_bytes(), "x.csv", GAVO, RA_TEXT, DEC_TEXT, FOV, MAXMAG, "RA +16 : 41 : 42, Dec x", 0
    )
    assert session.star_data("2026-11-03")["name"].startswith("field_FOV0.2_G16_GAVO_")


# --- picture settings recorded in the star data -------------------------------------------

DISPLAY = (
    "colsat=7.25; gamma=3.5; mini=0.04123105625617661; maxi=0.15 "
    "(shown as color saturation 7.25, nonlinear stretch 44.77, min 0.0017, max 0.0225)"
)


def test_star_data_records_and_restores_the_slider_positions(session: sc.Session) -> None:
    show(session, gavo_bytes())
    data = session.star_data("2026-11-03", DISPLAY)["data"]
    assert f"Display sliders: {DISPLAY}" in comments(data)
    expected = {"colsat": "7.25", "gamma": "3.5", "mini": "0.04123105625617661", "maxi": "0.15"}
    assert sc.display_from_file(data) == expected
    assert show(session, data, "uploaded", uploaded=True)["display"] == expected
    assert show(session, data, "saved copy", uploaded=False)["display"] is None  # only student uploads


def test_star_data_without_settings_still_works(session: sc.Session) -> None:
    show(session, gavo_bytes())
    data = session.star_data("2026-11-03")["data"]
    assert not any(line.startswith("Display sliders") for line in comments(data))
    assert show(session, data, "uploaded", uploaded=True)["display"] is None


@pytest.mark.parametrize(
    "line",
    [
        "colsat=5; gamma=2; mini=x; maxi=1",
        "colsat=5; gamma=2",
        "colsat=5; gamma=2; mini=0.1; maxi=0.2; extra=1",
        "colsat=5; gamma=2; mini=0.1; maxi=0.2</script><script>alert(1)",
        "colsat=nan; gamma=2; mini=0.1; maxi=0.2",
        "colsat=1e999; gamma=2; mini=0.1; maxi=0.2",
    ],
)
def test_malformed_or_tampered_settings_are_ignored(session: sc.Session, line: str) -> None:
    assert sc.display_settings(line) is None
    show(session, gavo_bytes())
    data = session.star_data("2026-11-03")["data"]
    tampered = data.replace(b"# Show Clusters", f"# Display sliders: {line}\n# Show Clusters".encode(), 1)
    assert show(session, tampered, "uploaded", uploaded=True)["display"] is None
    with pytest.raises(sc.InputError, match="could not be recorded"):
        session.star_data("2026-11-03", line)


# --- the field (the five form inputs) recorded in the downloads ----------------------------

FIELD_LINE = (
    "Field: Object Catalog Name=M13; Right Ascension (RA)=+16 : 41 : 42; Declination (Dec)=+36 : 27 : 39; "
    "Field of View (FOV)=0.2; Maximum Magnitude=16"
)
FIELD = {"object": "M13", "ra": "+16 : 41 : 42", "dec": "+36 : 27 : 39", "fov": 0.2, "maxmag": 16.0}


def test_field_is_recorded_in_both_csvs_and_read_back(session: sc.Session) -> None:
    show(session, gavo_bytes())
    stars = session.star_data("2026-11-03")["data"]
    points = session.hr_points("Apparent", "2026-11-03")["data"]
    assert FIELD_LINE in comments(stars) and FIELD_LINE in comments(points)
    assert sc.field_from_file(stars) == FIELD
    assert session.file_field(stars) == FIELD
    assert session.file_field(gavo_bytes()[0]) is None  # a plain GAVO/ESA file has no record


def test_hr_png_carries_the_field(session: sc.Session) -> None:
    from PIL import Image

    show(session, gavo_bytes())
    text = Image.open(io.BytesIO(session.hr_diagram("Absolute"))).text  # type: ignore[attr-defined]
    assert text["Object Catalog Name"] == "M13" and text["Right Ascension (RA)"] == "+16 : 41 : 42"
    assert text["Declination (Dec)"] == "+36 : 27 : 39"
    assert text["Field of View (FOV)"] == "0.2" and text["Maximum Magnitude"] == "16"
    assert "Absolute" in text["Title"] and "ESA/Gaia/DPAC" in text["Copyright"]


def test_custom_field_has_a_blank_object_name(session: sc.Session) -> None:
    title = f"RA {RA_TEXT}, Dec {DEC_TEXT}"
    session.show_stars(gavo_bytes(), "x.csv", GAVO, RA_TEXT, DEC_TEXT, FOV, MAXMAG, title, 0)
    record = sc.field_from_file(session.star_data("2026-11-03")["data"])
    assert record is not None and record["object"] == "" and record["ra"] == RA_TEXT


@pytest.mark.parametrize(
    "text",
    [
        FIELD_LINE.replace("+16 : 41 : 42", "+25 : 00 : 00"),  # RA out of range
        FIELD_LINE.replace("(FOV)=0.2", "(FOV)=91"),  # FOV out of range
        FIELD_LINE.replace("=M13", "=andromeda"),  # not a catalog name
        FIELD_LINE.replace("Maximum Magnitude", "Max Mag"),  # wrong label
        FIELD_LINE.rsplit("; ", 1)[0],  # missing value
        FIELD_LINE + "; Extra=1",
    ],
)
def test_bad_field_records_are_ignored(text: str) -> None:
    assert sc.field_settings(text.removeprefix("Field: ")) is None
    data = f"# {text}\nra,dec,parallax,phot_g_mean_mag,bp_rp,pmra,pmdec\n".encode()
    assert sc.field_from_file(data) is None
