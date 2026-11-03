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
if not (ROOT / "build" / "vendor" / "pyodide.js").is_file():
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


def fulfill_gavo(route) -> None:  # type: ignore[no-untyped-def]
    route.fulfill(
        status=200,
        body=(FIXTURES / "m13_gavo.csv").read_bytes(),
        headers={"content-type": "text/csv", "access-control-allow-origin": "*"},
    )


def test_progress_while_gavo_is_slow(browser, sites) -> None:  # type: ignore[no-untyped-def]
    """A slow answer shows a running clock and asks the student to wait, then clears."""
    held: list = []  # requests answered only when the test says so
    page = open_page(browser, sites["without"], lambda route: held.append(route))
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function(
        "document.getElementById('progress-text').textContent.startsWith('2 s so far')", timeout=10_000
    )
    assert page.is_visible("#progress")
    assert "keep this tab open" in page.text_content("#progress-text")
    assert page.get_attribute("#progress-text", "aria-hidden") == "true"  # not read out every second
    assert page.text_content("#status") == "Downloading star data from Gaia (via GAVO)…"
    assert len(held) == 1
    fulfill_gavo(held.pop())
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert page.is_hidden("#progress")
    assert page.text_content("#count") == "1,369"
    page.close()


def test_progress_counts_pieces_after_a_split(browser, sites) -> None:  # type: ignore[no-untyped-def]
    from fakegaia import TIMEOUT_BODY

    held: list = []

    def gavo(route) -> None:  # type: ignore[no-untyped-def]
        if not held and not getattr(gavo, "timed_out", False):
            gavo.timed_out = True  # type: ignore[attr-defined]
            route.fulfill(
                status=500,
                body=TIMEOUT_BODY,
                headers={"content-type": "text/xml", "access-control-allow-origin": "*"},
            )
            return
        held.append(route)

    page = open_page(browser, sites["without"], gavo)
    enter_m13(page)
    page.click("#show-stars")
    page.wait_for_function(
        "document.getElementById('status').textContent.includes('0 of 4 parts')", timeout=20_000
    )
    assert "too crowded to fetch at once" in page.text_content("#status")
    page.wait_for_function("document.getElementById('progress-bar').max === 4")
    for _ in range(100):  # the four pieces are requested together, not one after another
        if len(held) == 4:
            break
        page.wait_for_timeout(100)
    assert len(held) == 4, f"{len(held)} of the 4 pieces requested at once"
    fulfill_gavo(held.pop(0))
    page.wait_for_function(
        "document.getElementById('status').textContent.includes('1 of 4 parts')", timeout=20_000
    )
    assert page.evaluate("document.getElementById('progress-bar').value") == 1
    for route in held:
        fulfill_gavo(route)
    page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
    assert page.text_content("#count") == "1,369"  # the same stars as one query
    assert page.is_hidden("#progress")
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
    notes = page.evaluate("[...document.querySelectorAll('#notes li')].map((li) => li.textContent)")
    assert notes == [
        "Picture settings restored from the file.",
        "Field (object, RA, Dec, field of view, maximum magnitude) filled in from the file.",
    ]  # and no "only N of M stars inside" warning
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


@pytest.mark.parametrize("width", [1100, 390])  # desktop and phone layouts
def test_sliders_stay_put_while_values_change(browser, sites, width: int) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    page.set_viewport_size({"width": width, "height": 900})
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    ids = ("colsat", "gamma", "mini", "maxi")

    def boxes() -> list[tuple[float, float, float]]:
        return page.evaluate(
            "(ids) => ids.map((id) => { const r = document.getElementById(id).getBoundingClientRect();"
            " return [r.x, r.y, r.width]; })",
            list(ids),
        )

    before = boxes()
    texts = set()
    for slider in ids:
        for end in ("min", "max"):
            page.evaluate(
                "([id, end]) => { const s = document.getElementById(id); s.value = s[end];"
                " s.dispatchEvent(new Event('input')); }",
                [slider, end],
            )
            page.wait_for_timeout(50)
            texts.add(page.text_content(f"#{slider}-value"))
            clipped = page.evaluate(
                "(ids) => ids.filter((id) => { const o = document.getElementById(id + '-value');"
                " return o.scrollWidth > o.clientWidth; })",
                list(ids),
            )
            assert clipped == [], (slider, end, clipped)  # every value fits its column
            after = boxes()
            for (bx, by, bw), (ax, ay, aw) in zip(before, after, strict=True):
                assert (ax, ay, aw) == pytest.approx((bx, by, bw), abs=0.5), (slider, end, before, after)
    assert "59047.00" in texts  # the widest value (nonlinear stretch at its maximum) was shown
    page.close()


@pytest.mark.parametrize(
    ("width", "height"), [(1280, 800), (1440, 900), (1024, 768), (768, 1024), (390, 844)]
)
def test_picture_and_sliders_fit_on_screen_together(browser, sites, width: int, height: int) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    page.set_viewport_size({"width": width, "height": height})
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    span = page.evaluate(
        """() => {
            const ids = ['sky', 'colsat', 'gamma', 'mini', 'maxi', 'reset-sliders'];
            const boxes = ids.map((id) => document.getElementById(id).getBoundingClientRect());
            const min = (key) => Math.min(...boxes.map((b) => b[key]));
            const max = (key) => Math.max(...boxes.map((b) => b[key]));
            return { top: min('top'), bottom: max('bottom'), left: min('left'), right: max('right'),
                     sky: boxes[0].height };
        }"""
    )
    assert span["bottom"] - span["top"] <= height, span  # one scroll position shows them all
    assert span["left"] >= 0 and span["right"] <= width, span
    assert span["sky"] >= min(width, height) * 0.45, span  # the picture stays a useful size
    page.close()


def test_reset_puts_the_sliders_back(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    ids = ["colsat", "gamma", "mini", "maxi"]
    snapshot = "(ids) => ids.map((id) => document.getElementById(id).value)"
    pixels = "() => document.getElementById('sky').toDataURL()"
    initial, initial_picture = page.evaluate(snapshot, ids), page.evaluate(pixels)
    reset = page.get_by_role("button", name="Reset sliders")
    for slider in ids:
        page.evaluate(
            "(id) => { const s = document.getElementById(id); s.value = s.min;"
            " s.dispatchEvent(new Event('input')); }",
            slider,
        )
    page.wait_for_timeout(100)
    assert page.evaluate(snapshot, ids) != initial and page.evaluate(pixels) != initial_picture
    reset.click()
    page.wait_for_timeout(100)
    assert page.evaluate(snapshot, ids) == initial
    assert page.evaluate(pixels) == initial_picture  # the picture is back exactly as first drawn
    page.close()


def test_upload_option_is_part_of_choosing_a_field(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    inside = page.evaluate(
        "document.getElementById('upload').closest('section').querySelector('h2').textContent"
    )
    assert inside == "1. Choose a field"
    assert "nothing to install" not in (page.text_content(".intro") or "")
    page.close()


@pytest.mark.parametrize(
    ("fov", "shown"),
    [
        ("0.3", "degrees (18 arcminutes)"),
        ("0.2", "degrees (12 arcminutes)"),
        ("2", "degrees (120 arcminutes)"),
        ("0.125", "degrees (7.5 arcminutes)"),
        ("0.01", "degrees (0.6 arcminutes)"),
        ("", "degrees"),
        ("-1", "degrees"),
    ],
)
def test_fov_shows_arcminutes(browser, sites, fov: str, shown: str) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    assert page.text_content("#fov-unit") == "degrees (18 arcminutes)"  # the default 0.3
    page.fill("#fov", fov)  # fires input events as the student types
    assert page.text_content("#fov-unit") == shown
    page.close()


def viewer_state(page) -> dict:  # type: ignore[no-untyped-def]
    return page.evaluate(
        """() => {
            const box = (id) => document.getElementById(id).getBoundingClientRect();
            const viewer = document.querySelector('.viewer').getBoundingClientRect();
            const ids = ['sky', 'colsat', 'gamma', 'mini', 'maxi', 'reset-sliders', 'exit-fullscreen'];
            const parts = ids.map(box);
            const onScreen = (b) =>
                b.top >= 0 && b.left >= 0 && b.bottom <= innerHeight && b.right <= innerWidth;
            return {
                viewer: [viewer.left, viewer.top, viewer.width, viewer.height],
                sky: box('sky').width,
                inside: parts.every(onScreen),
                exitShown: !document.getElementById('exit-fullscreen').hidden,
                pressed: document.getElementById('fullscreen').getAttribute('aria-pressed'),
                focus: document.activeElement.id,
                bodyOverflow: getComputedStyle(document.body).overflow,
            };
        }"""
    )


@pytest.mark.parametrize(("width", "height"), [(1280, 800), (390, 844)])
def test_full_screen_picture_and_back(browser, sites, width: int, height: int) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    page.set_viewport_size({"width": width, "height": height})
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    normal = viewer_state(page)
    assert not normal["exitShown"] and normal["pressed"] == "false"

    page.get_by_role("button", name="Full screen").click()
    full = viewer_state(page)
    assert full["viewer"] == pytest.approx([0, 0, width, height], abs=1)  # covers the whole window
    assert full["inside"], full  # picture, every slider, Reset and Exit all on screen
    assert full["sky"] > normal["sky"] * 1.05, (normal["sky"], full["sky"])  # a bigger picture
    assert full["exitShown"] and full["focus"] == "exit-fullscreen" and full["pressed"] == "true"
    assert full["bodyOverflow"] == "hidden"  # the page behind does not scroll

    before = page.evaluate("document.getElementById('sky').toDataURL()")
    page.evaluate(
        "const s = document.getElementById('gamma'); s.value = s.max; s.dispatchEvent(new Event('input'));"
    )
    page.wait_for_timeout(100)
    assert (
        page.evaluate("document.getElementById('sky').toDataURL()") != before
    )  # sliders work in full screen

    page.get_by_role("button", name="Exit full screen").click()
    page.wait_for_function("!document.fullscreenElement")  # the browser leaves its full screen
    back = viewer_state(page)
    assert back["sky"] == pytest.approx(normal["sky"], abs=1) and not back["exitShown"]
    assert back["focus"] == "fullscreen" and back["pressed"] == "false" and back["bodyOverflow"] != "hidden"
    assert page.errors == []  # type: ignore[attr-defined]
    page.close()


def test_escape_leaves_full_screen(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    page.click("#fullscreen")
    assert viewer_state(page)["exitShown"]
    page.keyboard.press("Escape")
    state = viewer_state(page)
    assert not state["exitShown"] and state["pressed"] == "false" and state["focus"] == "fullscreen"
    page.close()


SLIDERS = ["colsat", "gamma", "mini", "maxi"]


def set_sliders(page, fractions: list[float]) -> None:  # type: ignore[no-untyped-def]
    page.evaluate(
        """([ids, fractions]) => ids.forEach((id, i) => {
            const s = document.getElementById(id);
            const lo = Number(s.min), hi = Number(s.max);
            s.value = String(lo + (hi - lo) * fractions[i]);
            s.dispatchEvent(new Event('input'));
        })""",
        [SLIDERS, fractions],
    )
    page.wait_for_timeout(100)


def test_star_data_brings_back_the_slider_positions(browser, sites, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    defaults = page.evaluate("(ids) => ids.map((id) => document.getElementById(id).value)", SLIDERS)
    set_sliders(page, [0.83, 0.37, 0.12, 0.41])
    tuned = page.evaluate("(ids) => ids.map((id) => document.getElementById(id).value)", SLIDERS)
    tuned_picture = page.evaluate("document.getElementById('sky').toDataURL()")
    with page.expect_download() as info:
        page.click("#data-download")
    path = tmp_path / info.value.suggested_filename
    info.value.save_as(path)  # downloads vanish when their page closes
    assert "# Display sliders: colsat=" in path.read_text()
    page.close()

    fresh = open_page(browser, sites["without"], gavo_ok([]))  # a new visit, e.g. the grader
    enter_m13(fresh)
    fresh.evaluate("document.getElementById('upload').open = true")
    fresh.set_input_files("#file", str(path))
    fresh.click("#use-file")
    fresh.wait_for_function(
        "document.getElementById('source').textContent.includes('uploaded')", timeout=60_000
    )
    assert fresh.evaluate("(ids) => ids.map((id) => document.getElementById(id).value)", SLIDERS) == tuned
    assert fresh.evaluate("document.getElementById('sky').toDataURL()") == tuned_picture  # pixel for pixel
    assert "Picture settings restored from the file." in (fresh.text_content("#notes") or "")
    fresh.click("#reset-sliders")
    assert fresh.evaluate("(ids) => ids.map((id) => document.getElementById(id).value)", SLIDERS) == defaults
    assert fresh.errors == []  # type: ignore[attr-defined]
    fresh.close()


def test_picture_png_records_settings_and_credit(browser, sites) -> None:  # type: ignore[no-untyped-def]
    import struct
    import zlib

    from PIL import Image

    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    set_sliders(page, [0.6, 0.25, 0.05, 0.5])
    exact = page.evaluate("(ids) => ids.map((id) => id + '=' + document.getElementById(id).value)", SLIDERS)
    with page.expect_download() as info:
        page.click("#sky-download")
    data = Path(info.value.path()).read_bytes()

    offset, chunks = 8, []  # every chunk's CRC-32 must be right
    while offset < len(data):
        (length,) = struct.unpack(">I", data[offset : offset + 4])
        kind, body = data[offset + 4 : offset + 8], data[offset + 8 : offset + 8 + length]
        (crc,) = struct.unpack(">I", data[offset + 8 + length : offset + 12 + length])
        assert zlib.crc32(kind + body) == crc, kind
        chunks.append(kind)
        offset += 12 + length
    assert chunks[0] == b"IHDR" and chunks[-1] == b"IEND" and chunks.count(b"tEXt") == 10

    image = Image.open(info.value.path())
    image.load()  # Pillow also checks the file
    assert image.size == (1200, 1200)
    text = image.text  # type: ignore[attr-defined]
    assert "; ".join(exact) in text["Description"]
    assert text["Title"].startswith("Show Clusters picture: M13, RA +16 : 41 : 42")
    assert text["Source"] == "Gaia DR3 via GAVO"
    assert "ESA/Gaia/DPAC" in text["Copyright"] and "HYG" in text["Copyright"]
    assert text["Object Catalog Name"] == "M13" and text["Right Ascension (RA)"] == "+16 : 41 : 42"
    assert text["Declination (Dec)"] == "+36 : 27 : 39"
    assert text["Field of View (FOV)"] == "0.2" and text["Maximum Magnitude"] == "16"
    page.close()


def test_uploading_star_data_fills_in_the_field(browser, sites, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    with page.expect_download() as info:
        page.click("#data-download")
    path = tmp_path / info.value.suggested_filename
    info.value.save_as(path)
    page.close()

    fresh = open_page(browser, sites["without"], gavo_ok([]))  # nothing typed in: RA/Dec 0, FOV 0.3, G 16
    fresh.evaluate("document.getElementById('upload').open = true")
    fresh.set_input_files("#file", str(path))
    fresh.click("#use-file")
    fresh.wait_for_function(
        "document.getElementById('source').textContent.includes('uploaded')", timeout=60_000
    )
    form = [fresh.input_value(f"#{i}") for i in ("object", "ra", "dec", "fov", "maxmag")]
    assert form == ["M13", "+16 : 41 : 42", "+36 : 27 : 39", "0.2", "16"]
    assert fresh.text_content("#fov-unit") == "degrees (12 arcminutes)"
    assert fresh.text_content("#count") == "1,369"
    assert "filled in from the file" in (fresh.text_content("#notes") or "")
    assert fresh.errors == []  # type: ignore[attr-defined]
    fresh.close()


CONTRAST = """(selectors) => {
    const rgb = (text) => text.match(/[\\d.]+/g).map(Number);
    const lum = ([r, g, b]) => {
        const f = (c) => { c /= 255; return c <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4; };
        return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b);
    };
    const background = (el) => {  // the first non-transparent background behind an element
        for (let node = el; node; node = node.parentElement) {
            const c = rgb(getComputedStyle(node).backgroundColor);
            if (c.length < 4 || c[3] > 0) return c.slice(0, 3);
        }
        return [255, 255, 255];
    };
    return selectors.flatMap((selector) => [...document.querySelectorAll(selector)].map((el) => {
        const [a, b] = [lum(rgb(getComputedStyle(el).color)), lum(background(el))];
        return [selector, (Math.max(a, b) + 0.05) / (Math.min(a, b) + 0.05)];
    }));
}"""


def test_full_screen_is_black_with_readable_dark_controls(browser, sites) -> None:  # type: ignore[no-untyped-def]
    page = open_page(browser, sites["without"], gavo_ok([]))
    enter_m13(page)
    page.click("#show-stars")
    stars_shown(page)
    page.click("#fullscreen")
    colors = page.evaluate(
        """() => ({
            viewer: getComputedStyle(document.querySelector('.viewer')).backgroundColor,
            panel: getComputedStyle(document.querySelector('.controls')).backgroundColor,
        })"""
    )
    assert colors["viewer"] == "rgb(0, 0, 0)"  # neutral black around the picture
    panel = [int(v) for v in colors["panel"].strip("rgb()").split(",")[:3]]
    assert max(panel) <= 48, colors  # a dark panel, not a bright one beside the stars
    selectors = [
        ".controls label",
        ".controls output",
        ".controls .credit",
        ".controls .credit a",
        "#reset-sliders",
        "#sky-download",
        "#data-download",
        "#exit-fullscreen",
    ]
    ratios = page.evaluate(CONTRAST, selectors)
    assert len(ratios) >= len(selectors)
    low = [(s, round(r, 2)) for s, r in ratios if r < 4.5]
    assert low == [], low  # WCAG AA for normal text
    page.click("#exit-fullscreen")
    page.wait_for_function("!document.fullscreenElement")
    normal = page.evaluate("getComputedStyle(document.querySelector('.controls')).backgroundColor")
    assert normal != colors["panel"]  # the page's normal look is unchanged
    page.close()
