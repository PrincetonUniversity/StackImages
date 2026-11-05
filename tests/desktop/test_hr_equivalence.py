"""The HR-diagram page's Python (web/hrimage.py) must reproduce the original ass.py.

Runs in the desktop container. Bob's functions are called unmodified, the way his window
calls them: get_RGB_image ("Open an RGB fits file: Select"), then GUI_HRdiag ("HR
diagram") and GUI_ColorDiag ("Color Diagram") with the two magnitude boxes. Only the wx
widgets, the file dialog and the plot window are replaced with stand-ins. The images are
synthetic star clusters (tests/web/hrfixtures.py).

Plotted points must be identical, not just close: positions, sizes, face and edge colours,
axis limits, titles and labels. If Bob changes the math in ass.py, these tests fail and
web/hrimage.py must be updated to match.
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
sys.path.insert(0, str(ROOT / "web"))
sys.path.insert(0, str(ROOT / "tests" / "web"))
sys.path.insert(0, str(ROOT))

import hrfixtures  # noqa: E402
import hrimage  # noqa: E402

GET_RGB_WIDGETS = (
    "lbl7", "lbl6a", "objectid", "objectra", "objectdec", "lbl_ra", "lbl_dec", "objectfov",
    "btn18", "lbl4", "lbl5", "lbl6", "pixres", "btn16", "btn12", "btn13", "btn14", "btn20",
)  # fmt: skip


class Widget:
    """Stands in for the wx controls these functions read and write."""

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


@pytest.fixture(scope="module")
def ass() -> Any:
    os.chdir(ROOT)  # ass.py opens csv/ relative to the working directory at import
    module: Any = importlib.import_module("ass")  # Bob's module; not type-checked
    return module


def open_with_ass(ass: Any, path: Path) -> SimpleNamespace:
    """Run get_RGB_image on `path`; return RGB, name and what slidershow was given."""
    for name in GET_RGB_WIDGETS:
        setattr(ass, name, Widget())
    captured: dict[str, Any] = {}
    ass.askopenfilenames = lambda **kwargs: str(path)
    ass.show_fwhm = lambda image: None  # fills the FWHM labels only
    ass.slidershow = lambda img, **kwargs: captured.setdefault("image", np.array(img))
    ass.close = lambda *a, **k: None
    ass.get_RGB_image(None)
    return SimpleNamespace(rgb=np.array(ass.RGB), name=ass.name, shown=captured["image"])


def diagram_with_ass(ass: Any, which: str, brightest: str, faintest: str, workdir: Path) -> Any:
    """Click HR diagram or Color Diagram; return the axes ass.py drew on."""
    from matplotlib.figure import Figure

    figures: list[Figure] = []

    def add(name: str) -> Figure:
        figures.append(Figure())
        return figures[-1]

    ass.plotter = SimpleNamespace(add=add)
    ass.GUI = True
    ass.mag0, ass.mag1 = Widget(brightest), Widget(faintest)
    (workdir / "txt").mkdir(exist_ok=True)  # makeHRdiag_new writes txt/color_mag.txt
    here = os.getcwd()
    os.chdir(workdir)  # the repository is mounted read-only in the container
    try:
        (ass.GUI_HRdiag if which == "hr" else ass.GUI_ColorDiag)(None)
    finally:
        os.chdir(here)
    return figures[-1].axes[-1]


def assert_same_axes(ours: Any, theirs: Any) -> int:
    """Identical scatter and axes; returns the number of points plotted."""
    mine, bobs = ours.collections[0], theirs.collections[0]
    for getter in ("get_offsets", "get_sizes", "get_facecolors", "get_edgecolors"):
        a, b = np.asarray(getattr(mine, getter)()), np.asarray(getattr(bobs, getter)())
        assert a.shape == b.shape, getter
        np.testing.assert_array_equal(a, b, err_msg=getter)
    assert ours.axis() == theirs.axis()
    for getter in ("get_title", "get_xlabel", "get_ylabel"):
        assert getattr(ours, getter)() == getattr(theirs, getter)(), getter
    return len(np.asarray(bobs.get_offsets()))


@pytest.fixture(scope="module")
def m2(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("fits") / "M2_Globular_Cluster_cropped.fit"
    header = {"OBJECT": "M 2", "CRVAL1": 323.3625, "CRVAL2": -0.8233}
    path.write_bytes(hrfixtures.fits_bytes(hrfixtures.cluster_rgb(450, 600, 700), header))
    return path


def test_open_matches_get_RGB_image(ass: Any, m2: Path) -> None:  # noqa: N802
    theirs = open_with_ass(ass, m2)
    ours = hrimage.read_rgb_fits(m2.read_bytes(), m2.name)
    assert ours.rgb.dtype == theirs.rgb.dtype
    np.testing.assert_array_equal(ours.rgb, theirs.rgb)
    assert ours.name == theirs.name == "M 2"
    assert ours.ra_hours == pytest.approx(323.3625 / 15) and ours.dec_deg == pytest.approx(-0.8233)
    np.testing.assert_array_equal(hrimage.slidershow_input(ours.rgb), theirs.shown)


@pytest.mark.parametrize("which", ["hr", "color"])
@pytest.mark.parametrize(
    ("brightest", "faintest"),
    [("", "20"), ("", "24"), ("8.5", "15"), ("12", "")],  # "" is an empty box: 10 and 20
)
def test_diagrams_match(
    ass: Any, m2: Path, tmp_path: Path, which: str, brightest: str, faintest: str
) -> None:
    open_with_ass(ass, m2)
    theirs = diagram_with_ass(ass, which, brightest, faintest, tmp_path)

    session = hrimage.Session()
    session.open_image(m2.read_bytes(), m2.name)
    points, minmag0, maxmag0 = session._points(brightest, faintest)
    assert session.image is not None
    if which == "hr":
        ours = hrimage.hr_figure(points, session.image.name, minmag0, maxmag0).axes[0]
    else:
        ours = hrimage.color_figure(points, session.image.name).axes[0]
    plotted = assert_same_axes(ours, theirs)
    assert plotted > (100 if (brightest, faintest) == ("", "20") else 0)  # a real diagram


@pytest.mark.parametrize(
    ("case", "expected_name"),
    [("mono", "NGC7089"), ("extension", "NGC7089"), ("filename", "NGC7089")],
)
def test_other_files_open_the_same(ass: Any, tmp_path: Path, case: str, expected_name: str) -> None:
    """A mono (2-D) image, an image in an extension HDU, and a name taken from the file name."""
    rgb = hrfixtures.cluster_rgb(300, 400, 200, seed=5)
    header: dict[str, object] = {} if case == "filename" else {"OBJECT": "NGC7089_session-2"}
    data = rgb[1] if case == "mono" else rgb
    path = tmp_path / "NGC7089-stack_final.fits"
    path.write_bytes(hrfixtures.fits_bytes(data, header, extension=case == "extension"))
    ass.objectra, ass.objectdec = Widget(), Widget()
    theirs = open_with_ass(ass, path)
    ours = hrimage.read_rgb_fits(path.read_bytes(), path.name)
    np.testing.assert_array_equal(ours.rgb, theirs.rgb)
    assert ours.name == theirs.name == expected_name
    np.testing.assert_array_equal(hrimage.slidershow_input(ours.rgb), theirs.shown)
