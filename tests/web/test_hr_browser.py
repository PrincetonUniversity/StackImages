"""End-to-end tests of "HR Diagram of Your Own Image" in Chromium (real Pyodide).

Builds the site, serves it locally and drives the page as a student would. The numbers
the browser produces are compared with web/hrimage.py run natively on the same file, so
the in-browser numpy/scipy (Pyodide) and the tested Python agree.
"""

from __future__ import annotations

import contextlib
import csv
import functools
import http.server
import io
import socketserver
import sys
import threading
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).parent))

playwright_api = pytest.importorskip("playwright.sync_api")
if not (ROOT / "build" / "vendor" / "pyodide.js").is_file():
    pytest.fail("build/vendor is missing: run scripts/vendor_pyodide.py before the browser tests")

import build_site  # noqa: E402
import hrfixtures  # noqa: E402
import hrimage  # noqa: E402

LOAD_TIMEOUT = 240_000  # Pyodide with scipy, matplotlib and astropy; slow in CI
WORK_TIMEOUT = 240_000
SEEN: dict[str, list[str]] = {"outside": [], "errors": []}  # filled by hr_page's listeners


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args: object) -> None:
        pass


@contextlib.contextmanager
def serve(directory: Path) -> Iterator[str]:
    handler = functools.partial(Quiet, directory=str(directory))
    with socketserver.ThreadingTCPServer(("127.0.0.1", 0), handler) as server:
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_address[1]}/"
        finally:
            server.shutdown()


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Iterator[str]:
    out = tmp_path_factory.mktemp("site") / "site"
    build_site.build(out, ROOT / "build" / "vendor", None)
    with serve(out) as url:
        yield url


@pytest.fixture(scope="module")
def m2_file(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("fits") / "M2_Globular_Cluster_cropped.fit"
    header = {"OBJECT": "M 2", "CRVAL1": 323.3625, "CRVAL2": -0.8233}
    path.write_bytes(hrfixtures.fits_bytes(hrfixtures.cluster_rgb(), header))
    return path


@pytest.fixture(scope="module")
def hr_page(browser: Any, site: str) -> Iterator[Any]:
    """One page for the whole module, as a student would use it (pytest-playwright's browser)."""
    context = browser.new_context(accept_downloads=True)
    page = context.new_page()
    page.on("pageerror", lambda error: SEEN["errors"].append(str(error)))
    page.on(
        "request",
        lambda request: None if request.url.startswith(site) else SEEN["outside"].append(request.url),
    )
    page.goto(f"{site}HRDiagram/")
    page.wait_for_function(
        "document.querySelector('#status').textContent.startsWith('Ready')", timeout=LOAD_TIMEOUT
    )
    yield page
    context.close()


def wait_idle(page: Any) -> None:
    page.wait_for_function("!document.querySelector('#file').disabled", timeout=WORK_TIMEOUT)


def open_file(page: Any, path: Path) -> None:
    page.set_input_files("#file", str(path))
    page.wait_for_function(
        "document.querySelector('#status').textContent.startsWith('Picture ready')"
        " || !document.querySelector('#error').hidden",
        timeout=WORK_TIMEOUT,
    )
    wait_idle(page)


def click_diagram(page: Any, button: str) -> None:
    page.click(button)
    page.wait_for_function(
        "!document.querySelector('#diagram').hidden || !document.querySelector('#error').hidden",
        timeout=WORK_TIMEOUT,
    )
    wait_idle(page)


def test_home_page_links_to_the_hr_page(site: str, hr_page: Any) -> None:
    page = hr_page
    page.goto(site)
    assert page.get_attribute("#open-hr", "href") == "HRDiagram/"
    page.click("#open-hr")
    page.wait_for_function(
        "document.querySelector('#status').textContent.startsWith('Ready')", timeout=LOAD_TIMEOUT
    )


def test_open_shows_the_picture(hr_page: Any, m2_file: Path) -> None:
    page = hr_page
    assert page.is_disabled("#hr-button")
    open_file(page, m2_file)
    assert page.is_hidden("#error")
    assert page.text_content("#object-name") == "M 2"
    assert page.text_content("#image-size") == "800 × 600 pixels"
    assert page.text_content("#position") == "+21 : 33 : 27 h, -00 : 49 : 24°"
    assert (page.eval_on_selector("#picture", "c => [c.width, c.height]")) == [800, 600]
    # The sliders start where slidershow starts them.
    native = hrimage.display_image(hrimage.read_rgb_fits(m2_file.read_bytes(), m2_file.name).rgb)
    mini = float(page.input_value("#mini"))
    maxi = float(page.input_value("#maxi"))
    assert mini == pytest.approx(np.sqrt(native["pct_lo"]), rel=1e-6)
    assert maxi == pytest.approx(np.sqrt(native["pct_hi"]), rel=1e-6)
    assert page.input_value("#colsat") == "5" and page.input_value("#gamma") == "2"
    lit = page.evaluate(
        """() => { const c = document.querySelector('#picture');
                   const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
                   let lit = 0; for (let i = 0; i < d.length; i += 4) if (d[i] + d[i+1] + d[i+2] > 60) lit++;
                   return lit; }"""
    )
    assert lit > 1000  # stars are drawn, not a blank canvas
    assert page.is_enabled("#hr-button") and page.is_enabled("#color-button")


def test_sliders_change_and_reset_the_picture(hr_page: Any) -> None:
    page = hr_page

    def pixels() -> str:
        return str(page.eval_on_selector("#picture", "c => c.toDataURL()"))

    start = pixels()
    page.fill("#gamma", "6")
    page.dispatch_event("#gamma", "input")
    page.wait_for_function("document.querySelector('#gamma-value').textContent !== '7.00'")
    assert pixels() != start
    page.click("#reset-sliders")
    assert pixels() == start


def browser_points(page: Any) -> list[list[float]]:
    with page.expect_download() as info:
        page.click("#points-download")
    text = Path(info.value.path()).read_text()
    rows = [line for line in text.splitlines() if not line.startswith("#")]
    return [[float(v) for v in row] for row in list(csv.reader(rows))[1:]]


def native_points(path: Path, brightest: str, faintest: str) -> list[list[float]]:
    session = hrimage.Session()
    session.open_image(path.read_bytes(), path.name)
    text = session.points(brightest, faintest)["data"].decode()
    rows = [line for line in text.splitlines() if not line.startswith("#")]
    return [[float(v) for v in row] for row in list(csv.reader(rows))[1:]]


def test_hr_diagram_matches_native_python(hr_page: Any, m2_file: Path) -> None:
    page = hr_page
    click_diagram(page, "#hr-button")
    assert page.is_hidden("#error")
    assert "Hertzsprung" in (page.get_attribute("#diagram", "alt") or "")
    summary = page.text_content("#diagram-summary") or ""
    ours = native_points(m2_file, "", "20")
    assert summary.startswith(f"{len(ours):,} stars plotted, of ")
    theirs = browser_points(page)
    assert len(theirs) == len(ours) > 100
    np.testing.assert_allclose(np.array(theirs), np.array(ours), rtol=1e-5, atol=2e-4)
    with page.expect_download() as info:
        page.click("#diagram-download")
    assert info.value.suggested_filename == "M_2_HR.png"
    from PIL import Image

    png = Image.open(io.BytesIO(Path(info.value.path()).read_bytes()))
    assert png.size == (1000, 880) and png.info["Title"] == "Hertzsprung-Russell Diagram for M 2"


def test_new_magnitudes_redraw_without_measuring_again(hr_page: Any, m2_file: Path) -> None:
    page = hr_page
    page.fill("#brightest", "8.5")
    page.fill("#faintest", "15")
    page.click("#hr-button")
    page.wait_for_function("document.querySelector('#status').textContent.length > 0")
    assert "Finding and measuring" not in (page.text_content("#status") or "")
    page.wait_for_function("!document.querySelector('#diagram').hidden", timeout=WORK_TIMEOUT)
    wait_idle(page)
    theirs = browser_points(page)
    ours = native_points(m2_file, "8.5", "15")
    assert len(theirs) == len(ours) > 0
    assert max(row[6] for row in theirs) <= 15 and min(row[6] for row in theirs) == pytest.approx(8.5)


def test_color_diagram(hr_page: Any) -> None:
    page = hr_page
    page.fill("#brightest", "")
    page.fill("#faintest", "20")
    click_diagram(page, "#color-button")
    assert page.is_hidden("#error")
    assert "Color diagram" in (page.get_attribute("#diagram", "alt") or "")
    with page.expect_download() as info:
        page.click("#diagram-download")
    assert info.value.suggested_filename == "M_2_Color.png"


def test_bad_magnitudes_are_explained(hr_page: Any) -> None:
    page = hr_page
    page.fill("#brightest", "12")
    page.fill("#faintest", "11")
    click_diagram(page, "#hr-button")
    assert "must be larger than" in (page.text_content("#error") or "")
    assert page.is_hidden("#diagram")
    page.fill("#brightest", "")
    page.fill("#faintest", "20")


def test_not_a_fits_file_is_explained(hr_page: Any, tmp_path: Path) -> None:
    page = hr_page
    junk = tmp_path / "notes.fits"
    junk.write_text("these are not pixels")
    open_file(page, junk)
    assert "not a FITS file" in (page.text_content("#error") or "")
    assert page.is_visible("#picture-empty") and page.is_disabled("#hr-button")


def test_too_large_file_is_refused_before_reading(hr_page: Any, tmp_path: Path) -> None:
    page = hr_page
    big = tmp_path / "huge.fits"
    with big.open("wb") as handle:
        handle.truncate(251 * 1024 * 1024)
    page.set_input_files("#file", str(big))
    page.wait_for_function("!document.querySelector('#error').hidden")
    assert "larger than 250 MB" in (page.text_content("#error") or "")


def test_nothing_left_the_page(hr_page: Any) -> None:
    page = hr_page
    assert page.url.endswith("/HRDiagram/")
    assert SEEN["outside"] == []  # nothing is fetched from, or sent to, another site
    assert SEEN["errors"] == []
