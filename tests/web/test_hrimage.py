"""web/hrimage.py natively: file checks, the page's session, downloads and the picture.

Agreement with ass.py itself is tested in tests/desktop/test_hr_equivalence.py.
"""

from __future__ import annotations

import io
from typing import Any

import numpy as np
import pytest
from PIL import Image

import hrfixtures
import hrimage


@pytest.fixture(scope="module")
def cluster() -> bytes:
    return hrfixtures.fits_bytes(hrfixtures.cluster_rgb(), {"OBJECT": "M 2"})


@pytest.fixture(scope="module")
def opened(cluster: bytes) -> tuple[hrimage.Session, dict[str, Any]]:
    session = hrimage.Session()
    return session, session.open_image(cluster, "M2_cropped.fit")


# --- Reading the file -----------------------------------------------------------------


def test_rgb_planes_become_float32_channels() -> None:
    data = hrfixtures.cluster_rgb(80, 100, 10)
    image = hrimage.read_rgb_fits(hrfixtures.fits_bytes(data), "x.fit")
    assert image.rgb.shape == (80, 100, 3) and image.rgb.dtype == np.float32
    for c in range(3):
        np.testing.assert_array_equal(image.rgb[:, :, c], data[c])


def test_mono_image_fills_all_three_channels() -> None:
    data = hrfixtures.cluster_rgb(80, 100, 10)[0]
    image = hrimage.read_rgb_fits(hrfixtures.fits_bytes(data), "x.fit")
    for c in range(3):
        np.testing.assert_array_equal(image.rgb[:, :, c], data)


def test_image_in_an_extension_is_found() -> None:
    data = hrfixtures.cluster_rgb(80, 100, 10)
    image = hrimage.read_rgb_fits(hrfixtures.fits_bytes(data, {"OBJECT": "M 13"}, extension=True), "x.fit")
    assert image.rgb.shape == (80, 100, 3) and image.name == "M 13"


@pytest.mark.parametrize(
    ("header", "filename", "name"),
    [
        ({"OBJECT": "M 2"}, "anything.fit", "M 2"),
        ({"OBJECT": "NGC7089_session-2"}, "anything.fit", "NGC7089"),  # ass.py: split('-'), split('_')
        ({}, "M2_Globular_Cluster_cropped.fit", "M2"),
        ({}, "folder/NGC6205-stack.fits", "NGC6205"),
    ],
)
def test_name_follows_get_rgb_image(header: dict[str, object], filename: str, name: str) -> None:
    data = hrfixtures.cluster_rgb(80, 100, 10)
    assert hrimage.read_rgb_fits(hrfixtures.fits_bytes(data, header), filename).name == name


def test_header_position_is_read() -> None:
    data = hrfixtures.cluster_rgb(80, 100, 10)
    image = hrimage.read_rgb_fits(hrfixtures.fits_bytes(data, {"CRVAL1": 323.4, "CRVAL2": -0.82}), "x.fit")
    assert image.ra_hours == pytest.approx(323.4 / 15) and image.dec_deg == pytest.approx(-0.82)


def _table_only() -> bytes:
    from astropy.io import fits

    out = io.BytesIO()
    table = fits.BinTableHDU.from_columns([fits.Column(name="a", format="E", array=np.zeros(3))])
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(out)
    return out.getvalue()


@pytest.mark.parametrize(
    ("data", "message"),
    [
        (b"not a fits file at all", "not a FITS file"),
        (b"SIMPLE  =                    T / damaged", "could not be read"),
        (_table_only(), "no image"),
        (hrfixtures.fits_bytes(np.zeros((2, 80, 100), np.float32)), "fewer than three"),
        (hrfixtures.fits_bytes(np.zeros((3, 40, 100), np.float32)), "too small"),
    ],
)
def test_unusable_files_are_explained(data: bytes, message: str) -> None:
    with pytest.raises(hrimage.InputError, match=message):
        hrimage.read_rgb_fits(data, "x.fit")


def test_missing_pixels_are_explained() -> None:
    data = hrfixtures.cluster_rgb(80, 100, 10)
    data[1, 40, 50] = np.nan
    with pytest.raises(hrimage.InputError, match="NaN"):
        hrimage.read_rgb_fits(hrfixtures.fits_bytes(data), "x.fit")


def test_size_limits(monkeypatch: pytest.MonkeyPatch) -> None:
    data = hrfixtures.fits_bytes(hrfixtures.cluster_rgb(80, 100, 10))
    monkeypatch.setattr(hrimage, "MAX_PIXELS", 80 * 100 - 1)
    with pytest.raises(hrimage.InputError, match="megapixels"):
        hrimage.read_rgb_fits(data, "x.fit")
    monkeypatch.setattr(hrimage, "MAX_PIXELS", 80 * 100)
    assert hrimage.read_rgb_fits(data, "x.fit").rgb.shape[:2] == (80, 100)
    monkeypatch.setattr(hrimage, "MAX_FILE_BYTES", len(data) - 1)
    with pytest.raises(hrimage.InputError, match="larger than"):
        hrimage.read_rgb_fits(data, "x.fit")


# --- Magnitude boxes --------------------------------------------------------------------


@pytest.mark.parametrize(("text", "value"), [("", 10.0), ("  ", 10.0), ("8.5", 8.5), ("12", 12.0)])
def test_empty_box_means_the_default(text: str, value: float) -> None:
    assert hrimage.parse_magnitude(text, 10, "Brightest star (mag)") == value


@pytest.mark.parametrize(
    ("text", "message"), [("ten", "must be a number"), ("1e9", "between"), ("nan", "between")]
)
def test_bad_magnitudes_are_explained(text: str, message: str) -> None:
    with pytest.raises(hrimage.InputError, match=message):
        hrimage.parse_magnitude(text, 10, "Brightest star (mag)")


def test_faintest_must_be_fainter(opened: tuple[hrimage.Session, dict[str, Any]]) -> None:
    session, _ = opened
    with pytest.raises(hrimage.InputError, match="larger than"):
        session.hr_diagram("12", "11")


def test_diagram_needs_an_image() -> None:
    with pytest.raises(hrimage.InputError, match="Open an RGB FITS image first"):
        hrimage.Session().hr_diagram("", "20")


# --- The session the page uses ----------------------------------------------------------


def test_open_returns_the_picture(opened: tuple[hrimage.Session, dict[str, Any]]) -> None:
    _, result = opened
    assert result["name"] == "M 2"
    assert result["pixels"] == [800, 600]
    assert (result["width"], result["height"], result["factor"]) == (800, 600, 1)
    assert len(result["image"]) == 800 * 600 * 3 * 4
    assert 0 < result["pct_lo"] < result["pct_hi"] <= 1


def test_large_pictures_are_shown_smaller() -> None:
    rgb = np.ones((1000, 3300, 3), np.float32)
    rgb[:, :, 1] = 2
    shown = hrimage.display_image(rgb, size=1600)
    assert (shown["factor"], shown["width"], shown["height"]) == (3, 1100, 333)
    assert len(shown["image"]) == 1100 * 333 * 3 * 4


def test_starting_stretch_uses_the_full_image() -> None:
    import showclusters as sc

    rgb = hrfixtures.cluster_rgb(400, 500, 300).transpose(1, 2, 0).copy()
    _, pct_lo, pct_hi = sc.display_base(hrimage.slidershow_input(rgb))
    shown = hrimage.display_image(rgb, size=100)
    assert shown["factor"] == 5
    assert (shown["pct_lo"], shown["pct_hi"]) == (pct_lo, pct_hi)


def test_hr_diagram_png(opened: tuple[hrimage.Session, dict[str, Any]]) -> None:
    session, _ = opened
    result = session.hr_diagram("", "")
    image = Image.open(io.BytesIO(result["png"]))
    assert image.size == (1000, 880)
    assert image.info["Title"] == "Hertzsprung-Russell Diagram for M 2"
    assert image.info["Comment"] == "Brightest star (mag)=10; Faintest star (mag)=20"
    assert 100 < result["plotted"] <= result["detected"]


def test_color_diagram_png(opened: tuple[hrimage.Session, dict[str, Any]]) -> None:
    session, _ = opened
    result = session.color_diagram("", "20")
    assert Image.open(io.BytesIO(result["png"])).info["Title"] == "Stellar Locus for M 2"
    assert result["plotted"] == session.hr_diagram("", "20")["plotted"]


def test_measured_once_per_image(cluster: bytes, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[int] = []
    original = hrimage.measure

    def counting(rgb: np.ndarray) -> hrimage.Measurements:
        calls.append(1)
        return original(rgb)

    monkeypatch.setattr(hrimage, "measure", counting)
    session = hrimage.Session()
    session.open_image(cluster, "a.fit")
    counts = [session.hr_diagram("", faint)["plotted"] for faint in ("14", "16", "24")]
    session.color_diagram("", "20")
    assert len(calls) == 1
    assert counts == sorted(counts) and counts[0] < counts[-1]  # fainter limit, more stars
    session.open_image(cluster, "b.fit")
    session.hr_diagram("", "20")
    assert len(calls) == 2  # a new image is measured again


def test_failed_open_forgets_the_previous_image(cluster: bytes) -> None:
    session = hrimage.Session()
    session.open_image(cluster, "a.fit")
    with pytest.raises(hrimage.InputError):
        session.open_image(b"junk", "b.fit")
    with pytest.raises(hrimage.InputError, match="Open an RGB FITS image first"):
        session.hr_diagram("", "20")


def test_no_stars_is_explained() -> None:
    flat = np.full((3, 100, 120), 50, np.float32)
    session = hrimage.Session()
    session.open_image(hrfixtures.fits_bytes(flat), "flat.fit")
    with pytest.raises(hrimage.InputError, match="No stars were found"):
        session.hr_diagram("", "20")


def test_points_csv_matches_the_diagram(opened: tuple[hrimage.Session, dict[str, Any]]) -> None:
    session, _ = opened
    result = session.points("", "18")
    assert result["name"] == "M_2_HR_points.csv"
    lines = result["data"].decode().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    assert "# Brightest star (mag)=10; Faintest star (mag)=18" in comments
    rows = [line for line in lines if not line.startswith("#")]
    assert rows[0] == "x,y,R_flux,G_flux,B_flux,color_index_B_minus_R,apparent_magnitude"
    assert len(rows) - 1 == session.hr_diagram("", "18")["plotted"]
    points = hrimage.plot_points(hrimage.measure(session.image.rgb), 10, 18)  # type: ignore[union-attr]
    first = rows[1].split(",")
    assert (int(first[0]), int(first[1])) == (int(points.x[0]), int(points.y[0]))
    assert float(first[6]) == pytest.approx(float(points.fluxmag[0]), abs=1e-4)
    mags = [float(row.split(",")[6]) for row in rows[1:]]
    assert max(mags) <= 18 and min(mags) == pytest.approx(10)
