"""The web port (web/showclusters.py) must reproduce the original ass.py.

Runs in the desktop container (wxPython, astroquery, Xvfb). The original functions are
called unmodified; only their GUI widgets and the Gaia network call are replaced with
fakes that return a recorded ESA result: M13 (0.2°, magnitude < 16), and Orion (60°,
magnitude < 7, the wide-field homework) centred on Alnilam, whose coordinates students type.
If Bob changes the math in ass.py, these tests fail and the port must be updated.
"""

from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "web" / "fixtures"
# name: (recorded ESA result, RA and Dec as typed or None to look up the name, FOV, magnitude)
CASES = {
    "m13": (FIXTURES / "m13_esa.vot", None, 0.2, 16.0),
    "orion60": (FIXTURES / "orion60_esa.vot", ("+05 : 36 : 13", "-01 : 12 : 07"), 60.0, 7.0),
}

sys.path.insert(0, str(ROOT / "web"))
sys.path.insert(0, str(ROOT))


class Widget:
    """Stands in for the wx controls show_gaia_image reads and writes."""

    def __init__(self, value: str = "") -> None:
        self.value, self.label, self.enabled = value, "", False

    def GetValue(self) -> str:  # noqa: N802 - wx API
        return self.value

    def SetValue(self, value: str) -> None:  # noqa: N802
        self.value = value

    def SetLabel(self, label: str) -> None:  # noqa: N802
        self.label = label

    def Enable(self) -> None:  # noqa: N802
        self.enabled = True


@pytest.fixture(scope="module", params=list(CASES))
def original(request: pytest.FixtureRequest) -> SimpleNamespace:
    """Run ass.py's Show Stars on a recorded ESA result; capture what it produces."""
    fixture, typed, fov, maxmag = CASES[request.param]
    from astropy.table import Table

    os.chdir(ROOT)  # ass.py opens csv/ and fit/ relative to the working directory
    from astroquery.gaia import Gaia

    ass: Any = importlib.import_module("ass")  # Bob's module; not type-checked

    table = Table.read(fixture, format="votable")
    queries: list[str] = []

    def fake_launch(query: str, *args: Any, **kwargs: Any) -> SimpleNamespace:
        queries.append(query)
        return SimpleNamespace(get_results=lambda: table)

    Gaia.launch_job_async = fake_launch  # type: ignore[method-assign]
    if typed is None:
        position = ass.getRaDec(request.param)
        typed = (ass.decimal_to_min_sec(position[0]), ass.decimal_to_min_sec(position[1]))
    widgets = {
        "objectid": Widget(""),
        "objectra": Widget(typed[0]),
        "objectdec": Widget(typed[1]),
        "objectfov": Widget(str(fov)),
        "objectmaxmag": Widget(str(maxmag)),
        "objectdt": Widget("0.0"),
        "lbl4": Widget(),
        "lbl10": Widget(),
        "btn14": Widget(),
        "btn20": Widget(),
    }
    for name, widget in widgets.items():
        setattr(ass, name, widget)
    captured: dict[str, Any] = {}
    ass.slidershow = lambda img, **kwargs: captured.setdefault("image", np.array(img))
    ass.suptitle = lambda *a, **k: None
    ass.close = lambda *a, **k: None
    ass.image_type = "2D"
    ass.objname = request.param

    ass.show_gaia_image(None)

    return SimpleNamespace(
        ass=ass,
        query=queries[0],
        typed=typed,
        fixture=fixture,
        fov=fov,
        maxmag=maxmag,
        magnitude_type=ass.abs_vs_app,
        stars={k: np.array(getattr(ass, k)) for k in ("ra", "dec", "mag", "absmag", "coloridx")},
        count_label=widgets["lbl4"].label,
        brightest_label=widgets["lbl10"].label,
        image=captured["image"],
    )


@pytest.fixture(scope="module")
def port(original: SimpleNamespace) -> SimpleNamespace:
    import showclusters as sc

    ra_mid = sc.parse_sexagesimal(original.typed[0])
    dec_mid = sc.parse_sexagesimal(original.typed[1])
    gaia = sc.parse_gaia_table(original.fixture.read_bytes(), original.fixture.name)
    bright = np.array(sc.bright_star_rows(ROOT / "csv" / "hygdata_v3.csv"))
    stars = sc.select_stars(gaia, bright, ra_mid, dec_mid, original.fov)
    frame = sc.frame_stars(stars, ra_mid, dec_mid, original.fov)
    return SimpleNamespace(sc=sc, ra_mid=ra_mid, dec_mid=dec_mid, frame=frame)


def test_lookup_matches_getRaDec(original: SimpleNamespace, port: SimpleNamespace) -> None:  # noqa: N802
    sc, ass = port.sc, original.ass
    for name in ("m13", "m44", "M92", "ngc5466", "NGC 6205", "ic4665"):
        ours = sc.lookup_object(name, ROOT / "csv")
        theirs = ass.getRaDec(name.replace(" ", ""))
        assert (ours.ra_hours, ours.dec_deg) == pytest.approx(theirs, abs=1e-12), name


def test_sexagesimal_round_trip_matches(original: SimpleNamespace, port: SimpleNamespace) -> None:
    for value in (0.0, 16.694861, -5.3912, 36.4608, -0.5, 23.99):
        text = original.ass.decimal_to_min_sec(value)
        assert port.sc.format_sexagesimal(value) == text
        assert port.sc.parse_sexagesimal(text) == pytest.approx(original.ass.min_sec_to_decimal(text))


def test_query_is_identical(original: SimpleNamespace, port: SimpleNamespace) -> None:
    ours = port.sc.gaia_adql(port.ra_mid, port.dec_mid, original.fov, original.maxmag, "esa")
    assert " ".join(ours.split()) == " ".join(original.query.split())


def test_same_stars_in_same_order(original: SimpleNamespace, port: SimpleNamespace) -> None:
    ours = port.frame.stars
    assert len(ours.mag) == len(original.stars["mag"]) > 1000
    for name in ("ra", "dec", "mag", "coloridx"):
        np.testing.assert_allclose(getattr(ours, name), original.stars[name], rtol=1e-6, err_msg=name)
    finite = np.isfinite(original.stars["absmag"])
    np.testing.assert_array_equal(np.isfinite(ours.absmag), finite)
    np.testing.assert_allclose(ours.absmag[finite], original.stars["absmag"][finite], rtol=1e-5)


def test_show_stars_magnitude_type(original: SimpleNamespace) -> None:
    """ass.py picks absolute magnitudes above 5° (its 'Abs'); the page shows that as Absolute."""
    assert original.magnitude_type == ("Abs" if original.fov > 5 else "Apparent")


def test_star_count_and_brightest_labels(original: SimpleNamespace, port: SimpleNamespace) -> None:
    assert int(original.count_label.split()[-1]) == port.frame.count
    assert float(original.brightest_label.split()[-1]) == pytest.approx(port.frame.brightest, abs=0.005)


def test_picture_matches_at_full_resolution(original: SimpleNamespace, port: SimpleNamespace) -> None:
    from astropy.io import fits

    psf = port.sc.prepare_psf(fits.getdata(ROOT / "fit" / "sxv_psf4-RGB.fit"))
    ours = port.sc.render_field(port.frame, psf, size=port.sc.ORIGINAL_SIZE)
    assert ours.shape == original.image.shape == (4000, 4000, 3)
    assert original.image.max() > 0
    np.testing.assert_allclose(ours, original.image, atol=1e-5)


@pytest.mark.parametrize(
    "magnitude_type", ["Apparent", "Absolute"]
)  # Absolute last: ass.py mutates its arrays
def test_hr_diagram_points_match(
    original: SimpleNamespace, port: SimpleNamespace, magnitude_type: str
) -> None:
    from matplotlib.figure import Figure

    ass = original.ass
    figures: list[Figure] = []

    def add(name: str) -> Figure:
        figures.append(Figure())
        return figures[-1]

    ass.plotter = SimpleNamespace(add=add)
    ass.GUI = True
    ass.abs_vs_app = magnitude_type
    ass.GUI_makeHRdiag0(None)
    collection: Any = figures[0].axes[0].collections[0]

    ours = port.sc.hr_arrays(port.frame.stars, magnitude_type)
    theirs_xy = np.asarray(collection.get_offsets())
    theirs_s = np.asarray(collection.get_sizes())
    keep = np.isfinite(theirs_xy).all(axis=1) & np.isfinite(theirs_s)
    np.testing.assert_allclose(ours["x"], theirs_xy[keep, 0], rtol=1e-6)
    np.testing.assert_allclose(ours["y"], theirs_xy[keep, 1], rtol=1e-6)
    np.testing.assert_allclose(ours["s"], theirs_s[keep], rtol=1e-5, atol=1e-6)
    theirs_rgb = np.asarray(collection.get_facecolors())[keep, :3]
    np.testing.assert_allclose(ours["c"][:, :3], theirs_rgb, atol=1e-6)
    ylim = figures[0].axes[0].get_ylim()
    assert sorted(ylim) == pytest.approx(sorted((ours["limits"][0] - 0.2, ours["limits"][1])), abs=1e-5)
