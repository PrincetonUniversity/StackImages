"""End-to-end tests of the built site in Chromium (real Pyodide, mocked Gaia services).

Needs build/vendor/ (scripts/vendor_pyodide.py) and Playwright's Chromium. The tests
build their own sites into a temporary directory, serve them locally and intercept
every request to GAVO, so no external service is contacted.
"""

from __future__ import annotations

import contextlib
import functools
import http.server
import json
import shutil
import socketserver
import sys
import threading
from collections.abc import Iterator
from pathlib import Path
from urllib.parse import parse_qs

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).parent))

playwright_api = pytest.importorskip("playwright.sync_api")
if not (ROOT / "build" / "vendor" / "pyodide.mjs").is_file():
    pytest.fail("build/vendor is missing: run scripts/vendor_pyodide.py before the browser tests")

import build_site  # noqa: E402

LOAD_TIMEOUT = 180_000  # first Pyodide + astropy start can be slow in CI
GAVO = "https://dc.g-vo.org/tap/sync"


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


@pytest.fixture(scope="session")
def sites(tmp_path_factory: pytest.TempPathFactory) -> Iterator[dict[str, str]]:
    base = tmp_path_factory.mktemp("sites")
    snapshots = base / "snapshots"
    snapshots.mkdir()
    shutil.copy2(FIXTURES / "m13_esa.csv", snapshots / "M13.csv")
    (snapshots / "index.json").write_text(
        json.dumps(
            {
                "M13": {
                    "file": "M13.csv",
                    "ra": "+16 : 41 : 42",
                    "dec": "+36 : 27 : 39",
                    "fov": 0.2,
                    "maxmag": 16,
                    "date": "2026-11-02",
                    "rows": 1498,
                }
            }
        )
    )
    build_site.build(base / "with", ROOT / "build" / "vendor", snapshots)
    build_site.build(base / "without", ROOT / "build" / "vendor", None)
    with serve(base / "with") as with_url, serve(base / "without") as without_url:
        yield {
            "with": f"{with_url}ShowClusters/",
            "without": f"{without_url}ShowClusters/",
            "root": without_url,
        }


def gavo_ok(requests: list[dict[str, list[str]]]):  # type: ignore[no-untyped-def]
    body = (FIXTURES / "m13_gavo.csv").read_bytes()

    def handle(route) -> None:  # type: ignore[no-untyped-def]
        requests.append(parse_qs(route.request.post_data or ""))
        route.fulfill(
            status=200, body=body, headers={"content-type": "text/csv", "access-control-allow-origin": "*"}
        )

    return handle


def gavo_down(route) -> None:  # type: ignore[no-untyped-def]
    route.fulfill(
        status=503,
        body="overloaded",
        headers={"content-type": "text/plain", "access-control-allow-origin": "*"},
    )


def open_page(browser, url: str, gavo_handler):  # type: ignore[no-untyped-def]
    page = browser.new_page(accept_downloads=True)
    errors: list[str] = []
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.route(f"{GAVO}*", gavo_handler)
    page.route(
        "**/*",
        lambda route: (
            route.abort()
            if not route.request.url.startswith(("http://127.0.0.1", GAVO))
            else route.fallback()
        ),
    )
    page.goto(url)
    page.wait_for_function(
        "document.getElementById('status').textContent.startsWith('Ready')", timeout=LOAD_TIMEOUT
    )
    page.errors = errors  # type: ignore[attr-defined]
    return page


def enter_m13(page) -> None:  # type: ignore[no-untyped-def]
    page.fill("#object", "m13")
    page.press("#object", "Enter")
    page.wait_for_function("document.getElementById('ra').value === '+16 : 41 : 42'")
    assert page.input_value("#dec") == "+36 : 27 : 39"
    page.fill("#fov", "0.2")
    page.fill("#maxmag", "16")


def canvas_lit_fraction(page) -> float:  # type: ignore[no-untyped-def]
    return page.evaluate(
        """() => {
            const c = document.getElementById('sky');
            const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
            let lit = 0;
            for (let i = 0; i < d.length; i += 4) if (d[i] + d[i + 1] + d[i + 2] > 30) lit++;
            return lit / (d.length / 4);
        }"""
    )


# `browser` is pytest-playwright's fixture (Chromium by default).


def test_homework_flow_with_live_gavo(browser, sites) -> None:  # type: ignore[no-untyped-def]
    requests: list[dict[str, list[str]]] = []
    page = open_page(browser, sites["without"], gavo_ok(requests))
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert page.text_content("#count") == "1,369"
    assert page.text_content("#brightest") == "11.36"
    assert page.text_content("#source") == "Gaia DR3 via GAVO"
    assert len(requests) == 1
    query = requests[0]["QUERY"][0]
    assert "gaia.dr3lite" in query and "250.425000" in query and "< 16.000000" in query
    assert requests[0]["MAXREC"] == ["300000"]
    lit = canvas_lit_fraction(page)
    assert 0.005 < lit < 0.6, lit

    before = canvas_lit_fraction(page)
    page.fill("#gamma", "6")
    page.dispatch_event("#gamma", "input")
    page.wait_for_timeout(300)
    assert canvas_lit_fraction(page) != before  # sliders change the picture

    for magnitude in ("Apparent", "Absolute"):
        page.select_option("#magtype", magnitude)
        page.click("#hr-button")
        page.wait_for_function(
            "(m) => document.getElementById('hr-download').download.endsWith(m + '.png')"
            " && !document.getElementById('hr').hidden && document.getElementById('hr').complete",
            arg=magnitude,
        )
        assert page.evaluate("document.getElementById('hr').naturalWidth") == 1000
        with page.expect_download() as info:
            page.click("#hr-download")
        download = info.value
        assert download.suggested_filename == f"M13_HR_{magnitude}.png"
        data = Path(download.path()).read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n" and len(data) > 20_000

    with page.expect_download() as info:
        page.click("#sky-download")
    assert Path(info.value.path()).read_bytes()[:4] == b"\x89PNG"
    assert page.errors == []  # type: ignore[attr-defined]
    page.close()


def test_show_stars_resets_magnitude_type(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    page.select_option("#magtype", "Absolute")
    assert page.input_value("#magtype") == "Absolute"
    page.click("#show-stars")
    # Show Stars puts the choice back to Apparent once the new stars are in (as ass.py does).
    page.wait_for_function("document.getElementById('magtype').value === 'Apparent'", timeout=60_000)
    assert page.text_content("#count") == "1,369"
    page.close()


def test_falls_back_to_saved_copy_when_gavo_is_down(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["with"], gavo_down)
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert page.text_content("#source") == "Gaia DR3 saved copy from ESA (2026-11-02)"
    assert page.text_content("#count") == "1,369"
    assert "saved copy was used" in page.text_content("#status")
    page.close()


def test_saved_copy_not_used_for_a_bigger_field(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["with"], gavo_down)
    enter_m13(page)
    page.fill("#fov", "0.5")
    page.click("#show-stars")
    page.wait_for_selector("#error:not([hidden])")
    assert "no saved copy covers this field" in page.text_content("#error")
    assert page.evaluate("document.getElementById('upload').open") is True
    assert page.is_hidden("#stars-content")
    page.close()


@pytest.mark.parametrize("fixture", ["m13_esa.vot", "m13_esa.csv", "m13_esa.fits"])
def test_upload_esa_file_when_gavo_is_down(browser, sites, fixture: str) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_down)
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_selector("#error:not([hidden])")
    assert "gea.esac.esa.int" not in page.text_content("#error")
    assert "FROM gaiadr3.gaia_source" in page.text_content("#adql")
    page.set_input_files("#file", str(FIXTURES / fixture))
    page.click("#use-file")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert page.text_content("#count") == "1,369"
    assert page.text_content("#source") == f"Gaia DR3 from uploaded ESA file ({fixture})"
    page.close()


def test_lookup_button_fills_in_coordinates(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    button = page.get_by_role("button", name="Look up")
    assert button.is_enabled()
    assert page.get_attribute("#object", "aria-describedby") == "object-hint"
    page.fill("#object", "ngc5466")
    button.click()
    page.wait_for_function("document.getElementById('ra').value !== '+00 : 00 : 00'")
    assert page.input_value("#ra").startswith("+14 : 05")
    assert "Found NGC5466" in (page.text_content("#status") or "")
    page.focus("#object")  # keyboard: Tab from the box reaches the button
    page.keyboard.press("Tab")
    assert page.evaluate("document.activeElement.id") == "lookup"
    page.close()


def test_split_download_has_no_outside_note(browser, sites) -> None:  # type: ignore[no-untyped-def]
    from fakegaia import FakeGaia

    page = open_page(browser, sites["without"], fake_gavo(FakeGaia(timeout_radius=0.06)))
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=90_000)
    assert page.text_content("#notes") == ""  # pieces overhang the field by design
    page.close()


def test_bad_inputs_show_messages(browser, sites, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    page.fill("#object", "andromeda")
    page.press("#object", "Enter")
    page.wait_for_selector("#error:not([hidden])")
    assert "Messier, NGC or IC" in page.text_content("#error")

    enter_m13(page)
    page.fill("#fov", "9")
    page.click("#show-stars")
    page.wait_for_function("document.getElementById('error').textContent.includes('Field of view')")

    page.fill("#fov", "0.2")
    bad = tmp_path / "notes.csv"
    bad.write_text("ra,dec\n1,2\n")
    page.set_input_files("#file", str(bad))
    page.evaluate("document.getElementById('upload').open = true")
    page.click("#use-file")
    page.wait_for_function(
        "document.getElementById('error').textContent.includes('missing these Gaia columns')"
    )
    assert page.errors == []  # type: ignore[attr-defined]
    page.close()


def test_no_requests_leave_the_site_except_gavo(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = browser.new_page()
    hosts: set[str] = set()
    page.on("request", lambda r: hosts.add(r.url.split("/")[2]))
    page.route(f"{GAVO}*", gavo_ok([]))
    page.goto(sites["without"])
    page.wait_for_function(
        "document.getElementById('status').textContent.startsWith('Ready')", timeout=LOAD_TIMEOUT
    )
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert {h for h in hosts if not h.startswith("127.0.0.1")} == {"dc.g-vo.org"}
    page.close()


def test_home_page_links_to_show_clusters(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = browser.new_page()
    page.goto(sites["root"])
    assert page.title() == "StackImages"
    assert page.evaluate("getComputedStyle(document.querySelector('.panel')).borderRadius") == "14px"
    page.click("#open-clusters")
    page.wait_for_url("**/ShowClusters/")
    assert page.title() == "Show Clusters with Gaia Data"
    page.click(".eyebrow a")  # back to the home page
    page.wait_for_function("document.title === 'StackImages'")
    page.close()


def test_all_steps_visible_before_any_data(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    for heading in ("1. Choose a field", "2. The stars", "3. HR diagram"):
        assert page.is_visible(f"h2:text-is('{heading}')"), heading
    assert page.is_visible("#stars-empty") and page.is_hidden("#stars-content")
    assert page.is_disabled("#hr-button")
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert page.is_hidden("#stars-empty") and page.is_enabled("#hr-button")
    page.close()


def fake_gavo(fake):  # type: ignore[no-untyped-def]
    """Route GAVO requests to a FakeGaia, which filters by each query's circle."""

    def handle(route) -> None:  # type: ignore[no-untyped-def]
        query = parse_qs(route.request.post_data or "")["QUERY"][0]
        status, body = fake.answer(query)
        kind = "text/csv" if status == 200 else "text/xml"
        route.fulfill(
            status=status, body=body, headers={"content-type": kind, "access-control-allow-origin": "*"}
        )

    return handle


def stars_shown(page) -> str:  # type: ignore[no-untyped-def]
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=90_000)
    return page.text_content("#count") or ""


@pytest.mark.parametrize("timeout_radius", [0.11, 0.06])  # one level of splitting, then two
def test_split_field_shows_the_same_result_as_one_query(browser, sites, timeout_radius: float) -> None:  # type: ignore[no-untyped-def]
    from fakegaia import FakeGaia

    results = []
    for fake in (FakeGaia(), FakeGaia(timeout_radius=timeout_radius)):
        page = open_page(browser, sites["without"], fake_gavo(fake))
        enter_m13(page)
        page.click("#show-stars")
        results.append((stars_shown(page), page.text_content("#brightest"), len(fake.queries)))
        page.click("#hr-button")
        page.wait_for_function("!document.getElementById('hr').hidden")
        assert page.errors == []  # type: ignore[attr-defined]
        page.close()
    (whole_count, whole_bright, whole_queries), (split_count, split_bright, split_queries) = results
    assert whole_queries == 1 and split_queries > 1
    assert split_count == whole_count == "1,369"
    assert split_bright == whole_bright


def test_too_dense_field_says_so(browser, sites) -> None:  # type: ignore[no-untyped-def]
    from fakegaia import FakeGaia

    fake = FakeGaia(timeout_radius=0.0)  # every query times out, however small
    page = open_page(browser, sites["without"], fake_gavo(fake))
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_selector("#error:not([hidden])", timeout=90_000)
    assert "too many stars" in (page.text_content("#error") or "")
    assert "unavailable" not in (page.text_content("#error") or "")
    assert len(fake.queries) > 5  # it split before giving up
    page.close()


def test_star_data_download_and_reupload(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    count = stars_shown(page)
    assert "ESA/Gaia/DPAC" in (page.text_content("#stars-content .credit") or "")
    with page.expect_download() as info:
        page.click("#data-download")
    download = info.value
    assert download.suggested_filename.startswith(
        "M13_FOV0.2_G16_GAVO_"
    ) and download.suggested_filename.endswith("_stars.csv")
    text = Path(download.path()).read_text()
    assert text.startswith("# Show Clusters star data: M13") and "CC BY-NC 3.0 IGO" in text

    page.evaluate("document.getElementById('upload').open = true")
    page.set_input_files("#file", download.path())
    page.click("#use-file")
    page.wait_for_function(
        "document.getElementById('source').textContent.includes('uploaded')", timeout=60_000
    )
    assert page.text_content("#count") == count == "1,369"
    assert page.text_content("#notes") == ""
    assert page.errors == []  # type: ignore[attr-defined]
    page.close()


def test_hr_points_download_follows_the_diagram_shown(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    assert page.is_hidden("#hr-csv-download") and page.is_hidden("#hr-credit")
    page.select_option("#magtype", "Absolute")
    page.click("#hr-button")
    page.wait_for_function("!document.getElementById('hr').hidden")
    assert "HYG database v3" in (page.text_content("#hr-credit") or "")
    page.select_option("#magtype", "Apparent")  # menu changed, diagram not redrawn
    with page.expect_download() as info:
        page.click("#hr-csv-download")
    download = info.value
    assert download.suggested_filename == "M13_HR_Absolute_points.csv"
    lines = Path(download.path()).read_text().splitlines()
    header = next(line for line in lines if not line.startswith("#"))
    assert header.split(",") == [
        "catalog",
        "ra_deg",
        "dec_deg",
        "parallax_mas",
        "apparent_mag",
        "absolute_mag",
        "color_index",
        "x_plotted",
        "y_plotted",
        "drawn",
    ]
    assert sum(1 for line in lines if not line.startswith("#")) - 1 == 1369
    assert any("CC BY-SA 2.5" in line for line in lines if line.startswith("#"))
    page.close()


def test_footer_credits_gaia_and_hyg(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    footer = page.text_content("footer") or ""
    for text in (
        "European Space Agency (ESA) mission",
        "Gaia Multilateral Agreement",
        "CC BY-NC 3.0 IGO",
        "HYG database",
        "CC BY-SA 2.5",
        "GAVO",
    ):
        assert text in footer, text
    page.close()
