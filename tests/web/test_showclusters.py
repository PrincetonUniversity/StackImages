"""Unit tests for web/showclusters.py (native Python; the same code runs in Pyodide)."""

from __future__ import annotations

import io
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(ROOT / "web"))

import showclusters as sc  # noqa: E402

CSV_DIR = ROOT / "csv"
M13 = sc.lookup_object("m13", CSV_DIR)
# The desktop app queries with the RA/Dec as shown in its text boxes (whole seconds).
RA = sc.parse_sexagesimal(sc.format_sexagesimal(M13.ra_hours))
DEC = sc.parse_sexagesimal(sc.format_sexagesimal(M13.dec_deg))


@pytest.fixture(scope="module")
def bright() -> np.ndarray:
    return np.array(sc.bright_star_rows(CSV_DIR / "hygdata_v3.csv"))


def load(name: str) -> dict[str, np.ndarray]:
    return sc.parse_gaia_table((FIXTURES / name).read_bytes(), name)


# --- lookup and coordinates ------------------------------------------------------------


@pytest.mark.parametrize(
    ("name", "canonical", "ra", "dec"),
    [
        ("m13", "M13", 16.6949, 36.4608),
        ("M 44", "M44", 8.6700, 19.6667),
        ("NGC5466", "NGC5466", 14.0920, 28.5342),
        ("ic4665", "IC4665", 17.7700, 5.7167),
    ],
)
def test_lookup_known_objects(name: str, canonical: str, ra: float, dec: float) -> None:
    found = sc.lookup_object(name, CSV_DIR)
    assert found.name == canonical
    assert found.ra_hours == pytest.approx(ra, abs=0.01)
    assert found.dec_deg == pytest.approx(dec, abs=0.01)


@pytest.mark.parametrize("name", ["", "m", "m13; DROP", "andromeda", "m12345", "ngc-5", "x13"])
def test_lookup_rejects_bad_names(name: str) -> None:
    with pytest.raises(sc.InputError, match="Messier, NGC or IC"):
        sc.lookup_object(name, CSV_DIR)


def test_lookup_unknown_number() -> None:
    with pytest.raises(sc.InputError, match="M999 is not in the catalog"):
        sc.lookup_object("m999", CSV_DIR)


def test_negative_declination_uses_sign_of_text() -> None:
    m4 = sc.lookup_object("m4", CSV_DIR)
    assert m4.dec_deg == pytest.approx(-26.53, abs=0.01)


@pytest.mark.parametrize("value", [0.0, 16.694861, -26.5258, 89.9, -0.25])
def test_sexagesimal_round_trip(value: float) -> None:
    assert sc.parse_sexagesimal(sc.format_sexagesimal(value)) == pytest.approx(value, abs=1 / 3600)


@pytest.mark.parametrize("text", ["", "16:41", "a : b : c", "+16 : 75 : 00", "+16 : 41 : 61", "1:2:3:4"])
def test_sexagesimal_rejects_bad_text(text: str) -> None:
    with pytest.raises(sc.InputError):
        sc.parse_sexagesimal(text)


@pytest.mark.parametrize(
    ("ra", "dec", "fov", "maxmag", "message"),
    [
        (24.0, 0, 0.3, 16, "Right ascension"),
        (-1, 0, 0.3, 16, "Right ascension"),
        (1, 91, 0.3, 16, "Declination"),
        (1, 0, 0.0, 16, "Field of view"),
        (1, 0, 6.0, 16, "Field of view"),
        (1, 0, 0.3, 22, "Maximum magnitude"),
        (1, 0, 0.3, 3, "Maximum magnitude"),
        (float("nan"), 0, 0.3, 16, "numbers"),
        (1, 0, float("inf"), 16, "numbers"),
    ],
)
def test_query_inputs_are_validated(ra: float, dec: float, fov: float, maxmag: float, message: str) -> None:
    with pytest.raises(sc.InputError, match=message):
        sc.gaia_adql(ra, dec, fov, maxmag, "gavo")


def test_adql_for_both_services() -> None:
    esa = sc.gaia_adql(RA, DEC, 0.2, 16, "esa")
    gavo = sc.gaia_adql(RA, DEC, 0.2, 16, "gavo")
    assert "FROM gaiadr3.gaia_source AS gaia" in esa and "gaia.bp_rp" in esa
    assert "FROM gaia.dr3lite AS gaia" in gavo and "phot_bp_mean_mag" in gavo
    for query in (esa, gavo):
        assert "250.425000" in query and "36.460833" in query and "0.141421" in query
        assert "phot_g_mean_mag < 16.000000" in query
    with pytest.raises(ValueError):
        sc.gaia_adql(RA, DEC, 0.2, 16, "other")


# --- reading Gaia results --------------------------------------------------------------


@pytest.mark.parametrize("name", ["m13_gavo.csv", "m13_esa.csv", "m13_esa.vot", "m13_esa.fits"])
def test_every_format_reads_all_rows(name: str) -> None:
    table = load(name)
    assert set(table) == set(sc.REQUIRED)
    assert len(table["ra"]) == 1498
    assert np.nanmax(table["phot_g_mean_mag"]) < 16


def test_esa_formats_agree_exactly() -> None:
    # Each file came from a separate query, and ESA returns rows in no fixed order.
    tables = [load(name) for name in ("m13_esa.csv", "m13_esa.vot", "m13_esa.fits")]
    orders = [np.lexsort((t["dec"], t["ra"])) for t in tables]
    csv_t, *others = [{k: v[o] for k, v in t.items()} for t, o in zip(tables, orders, strict=True)]
    for other in others:
        for name in sc.REQUIRED:
            np.testing.assert_allclose(other[name], csv_t[name], rtol=1e-6, equal_nan=True, err_msg=name)


def test_gavo_gives_same_answer_as_esa(bright: np.ndarray) -> None:
    results = []
    for name in ("m13_gavo.csv", "m13_esa.vot"):
        stars = sc.select_stars(load(name), bright, RA, DEC, 0.2)
        results.append(sc.frame_stars(stars, RA, DEC, 0.2))
    gavo, esa = results
    assert gavo.count == esa.count == 1369
    assert gavo.brightest == pytest.approx(esa.brightest, abs=1e-4)
    order = np.lexsort((gavo.stars.dec, gavo.stars.ra)), np.lexsort((esa.stars.dec, esa.stars.ra))
    np.testing.assert_allclose(gavo.stars.mag[order[0]], esa.stars.mag[order[1]], rtol=1e-6)
    np.testing.assert_allclose(gavo.stars.coloridx[order[0]], esa.stars.coloridx[order[1]], atol=1e-5)


@pytest.mark.parametrize(
    ("data", "message"),
    [
        (b"", "empty"),
        (b"   \n", "empty"),
        (b"ra,dec\n1,2\n", "missing these Gaia columns: parallax"),
        (b"ra,dec,parallax,phot_g_mean_mag,bp_rp,pmra,pmdec\n1,2,3\n", "Line 2 of the CSV has 3 fields"),
        (b"\xff\xfe\x00garbage", "not a CSV"),
        (b"<?xml version='1.0'?><nope/>", "could not be read as VOTABLE"),
        (b"SIMPLE  = T broken", "could not be read as FITS"),
    ],
)
def test_bad_uploads_are_rejected_with_a_reason(data: bytes, message: str) -> None:
    with pytest.raises(sc.InputError, match=message):
        sc.parse_gaia_table(data, "upload.csv" if not data.startswith((b"<", b"S")) else "upload")


def test_oversized_upload_is_rejected() -> None:
    with pytest.raises(sc.InputError, match="larger than 50 MB"):
        sc.parse_gaia_table(b"x" * (sc.MAX_UPLOAD_BYTES + 1), "big.csv")


def test_too_many_rows_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sc, "MAX_ROWS", 100)
    with pytest.raises(sc.InputError, match="more than 100 stars"):
        load("m13_esa.csv")


def test_missing_values_drop_the_star(bright: np.ndarray) -> None:
    text = (FIXTURES / "m13_esa.csv").read_text().splitlines()
    header, rows = text[0], text[1:6]
    fields = rows[0].split(",")
    fields[5] = ""  # pmra missing → ass.py skips the star
    data = "\n".join([header, ",".join(fields), *rows[1:]]).encode()
    table = sc.parse_gaia_table(data, "x.csv")
    stars = sc.select_stars(table, np.empty((0, 7)), RA, DEC, 0.2)
    assert len(stars.mag) == 4


# --- selection, picture and HR diagram -------------------------------------------------


def test_bright_star_subset_is_small_and_in_range(bright: np.ndarray) -> None:
    assert 2000 < len(bright) < 4000
    assert bright[:, 2].min() > -5 and bright[:, 2].max() < 5.5
    assert ((bright[:, 0] >= 0) & (bright[:, 0] < 24)).all()


def test_bright_stars_join_when_in_field(bright: np.ndarray) -> None:
    vega = sc.lookup_object("m57", CSV_DIR)  # Ring Nebula, about 4.6° from Vega
    gaia = {name: np.empty(0) for name in sc.REQUIRED}
    near = sc.select_stars(gaia, bright, vega.ra_hours, vega.dec_deg, 5.0)
    far = sc.select_stars(gaia, bright, vega.ra_hours, vega.dec_deg, 0.1)
    assert len(near.mag) > len(far.mag)
    assert (near.mag < 5.5).all()


def test_render_and_display(bright: np.ndarray) -> None:
    from astropy.io import fits

    stars = sc.select_stars(load("m13_gavo.csv"), bright, RA, DEC, 0.2)
    frame = sc.frame_stars(stars, RA, DEC, 0.2)
    psf = sc.prepare_psf(fits.getdata(ROOT / "fit" / "sxv_psf4-RGB.fit"))
    assert psf.shape == (500, 500, 3) and psf.max() == pytest.approx(1)
    image = sc.render_field(frame, psf, 600)
    assert image.shape == (600, 600, 3) and image.dtype == np.float32
    assert image.max() > 0.005 and (image >= 0).all()
    lit = (image.sum(axis=2) > 1e-4).mean()
    assert 0.01 < lit < 0.5, lit  # stars, not a blank or saturated frame
    base, lo, hi = sc.display_base(image)
    assert base.shape == image.shape and 0 < lo < hi <= 1


def test_scaled_psf_keeps_peak_and_shape() -> None:
    psf = np.zeros((500, 500, 3), dtype=np.float32)
    psf[245:255, 245:255] = 1
    small = sc.scale_psf(psf, 1200)
    assert small.shape == (150, 150, 3)
    assert small.max() == pytest.approx(1)
    assert small[75, 75, 0] > 0.5 and small[0, 0, 0] == 0


@pytest.mark.parametrize("magnitude_type", ["Apparent", "Absolute"])
def test_hr_diagram_png(bright: np.ndarray, magnitude_type: str) -> None:
    from PIL import Image

    stars = sc.select_stars(load("m13_gavo.csv"), bright, RA, DEC, 0.2)
    kept = sc.frame_stars(stars, RA, DEC, 0.2).stars
    data = sc.hr_arrays(kept, magnitude_type)
    assert len(data["x"]) == len(kept.mag)
    assert data["s"][0] == 0 and data["s"].max() == pytest.approx(12)
    png = sc.hr_png(kept, magnitude_type, "M13", "Gaia DR3 via GAVO")
    image = Image.open(io.BytesIO(png))
    assert image.format == "PNG" and image.size == (1000, 800)
    colors = image.convert("RGB").getcolors(maxcolors=1_000_000) or []
    assert len(colors) > 50  # points were drawn, not just a frame


def test_hr_diagram_input_errors(bright: np.ndarray) -> None:
    empty = sc.select_stars({n: np.empty(0) for n in sc.REQUIRED}, np.empty((0, 7)), RA, DEC, 0.2)
    with pytest.raises(sc.InputError, match="No stars"):
        sc.hr_arrays(empty, "Apparent")
    stars = sc.select_stars(load("m13_gavo.csv"), bright, RA, DEC, 0.2)
    with pytest.raises(sc.InputError, match="Apparent or Absolute"):
        sc.hr_arrays(stars, "Relative")
