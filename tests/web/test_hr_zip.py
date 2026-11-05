"""The downloadable copy of "HR Diagram of Your Own Image" (scripts/zip_hr.py): the page on
its own, with its own Pyodide, served the way a personal web site serves it."""

from __future__ import annotations

import contextlib
import re
import shutil
import socketserver
import sys
import threading
import zipfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).parent))

from test_hr_browser import LOAD_TIMEOUT, WORK_TIMEOUT  # noqa: E402
from test_site_zip import HOST, handler_for, watch  # noqa: E402

import build_site  # noqa: E402
import hrfixtures  # noqa: E402
import zip_hr  # noqa: E402
import zip_site  # noqa: E402


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not (ROOT / "build" / "vendor" / "pyodide.js").is_file():
        pytest.fail("build/vendor is missing: run scripts/vendor_pyodide.py first")
    site = tmp_path_factory.mktemp("hrzip") / "site"
    build_site.build(site, ROOT / "build" / "vendor", None)
    return site


@pytest.fixture(scope="module")
def archive(site: Path) -> Path:
    zip_hr.write_zip(site, site / zip_hr.NAME)  # where pages.yml puts it
    return site / zip_hr.NAME


def test_zip_holds_the_page_and_its_own_pyodide(site: Path, archive: Path) -> None:
    names = zipfile.ZipFile(archive).namelist()
    assert all(n.startswith("HRDiagram/") for n in names)
    page = {
        f"HRDiagram/{p.relative_to(site / 'HRDiagram').as_posix()}" for p in (site / "HRDiagram").rglob("*")
    }
    vendor = site / "ShowClusters" / "vendor"
    pyodide = {
        f"HRDiagram/vendor/{p.relative_to(vendor).as_posix()}" for p in vendor.rglob("*") if p.is_file()
    }
    assert set(names) == page | pyodide | {"HRDiagram/HOSTING.txt"}
    assert "HRDiagram/vendor/pyodide.asm.wasm" in names
    assert not any("data/" in n or "ShowClusters" in n for n in names)  # none of Show Clusters' data
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert z.read("HRDiagram/hrimage.py") == (site / "HRDiagram" / "hrimage.py").read_bytes()
    assert archive.stat().st_size < 55e6  # Pyodide only, not the 150 MB site


def test_the_unzipped_page_needs_nothing_outside_its_folder(archive: Path) -> None:
    with zipfile.ZipFile(archive) as z:
        worker = z.read("HRDiagram/worker.js").decode()
        page = z.read("HRDiagram/index.html").decode()
        app = z.read("HRDiagram/app.js").decode()
    assert 'new URL("vendor/", self.location.href)' in worker
    assert 'href="https://princetonuniversity.github.io/StackImages/"' in page  # the StackImages link
    for name, text in (("worker.js", worker), ("index.html", page), ("app.js", app)):
        assert not re.search(r"""["'`]\.\./""", text), f"{name} still points outside the folder"


def test_rewrites_fail_loudly_when_the_page_changes(site: Path, tmp_path: Path) -> None:
    copy = tmp_path / "site"
    shutil.copytree(site, copy)
    worker = copy / "HRDiagram" / "worker.js"
    worker.write_text(worker.read_text().replace('"../ShowClusters/vendor/"', '"../elsewhere/"'))
    with pytest.raises(SystemExit, match="worker.js"):
        zip_hr.write_zip(copy, tmp_path / "x.zip")
    assert not (tmp_path / "x.zip").exists() and not (tmp_path / "x.zip.partial").exists()


def test_zip_is_reproducible_and_fingerprinted(site: Path, archive: Path, tmp_path: Path) -> None:
    again = tmp_path / "again.zip"
    zip_hr.write_zip(site, again)
    assert again.read_bytes() == archive.read_bytes()
    assert zip_hr.main(["--fingerprint", str(archive)]) == 0
    copy = tmp_path / "site"
    shutil.copytree(site, copy)
    (copy / "HRDiagram" / "hrimage.py").write_text("# changed\n")
    zip_hr.write_zip(copy, again)
    assert zip_site.fingerprint(again) != zip_site.fingerprint(archive)


def test_site_zip_leaves_out_the_hr_zip(site: Path, archive: Path, tmp_path: Path) -> None:
    """Both zips sit in _site/; the site's zip must not carry the HR one inside it."""
    out = tmp_path / "site.zip"
    zip_site.write_zip(site, out)
    assert f"StackImages/{zip_hr.NAME}" not in zipfile.ZipFile(out).namelist()
    hosting = zipfile.ZipFile(out).read("StackImages/HOSTING.txt").decode()
    assert "StackImages/HRDiagram/" in hosting  # the site's copy names both pages


def test_dry_run_writes_nothing(site: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "x.zip"
    assert zip_hr.main(["--site", str(site), "--out", str(out), "--dry-run"]) == 0
    assert not out.exists()
    printed = capsys.readouterr().out
    assert "HRDiagram/HOSTING.txt" in printed and "HRDiagram/vendor/pyodide.js" in printed


def test_refuses_a_directory_that_is_not_a_site(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="not a built site"):
        zip_hr.write_zip(tmp_path, tmp_path / "x.zip")


# --- the unzipped copy in a browser, on a server like Bob's ------------------------------


@contextlib.contextmanager
def personal_site(archive: Path, root: Path) -> Iterator[str]:
    """Unzip into ~rvdb/ of a fresh docroot and serve it; yields the page's URL."""
    with zipfile.ZipFile(archive) as z:
        z.extractall(root / "~rvdb")
    with socketserver.ThreadingTCPServer(("127.0.0.1", 0), handler_for(root, HOST["types"])) as server:
        threading.Thread(target=server.serve_forever, daemon=True).start()
        try:
            yield f"http://127.0.0.1:{server.server_address[1]}/~rvdb/HRDiagram/"
        finally:
            server.shutdown()


def test_unzipped_copy_draws_a_diagram_under_a_strict_host(
    browser: Any, archive: Path, tmp_path: Path
) -> None:
    """Bob's site's headers (nosniff, its CSP, octet-stream for .wasm/.whl): no server change."""
    image = tmp_path / "M2.fit"
    image.write_bytes(hrfixtures.fits_bytes(hrfixtures.cluster_rgb(), {"OBJECT": "M 2"}))
    with personal_site(archive, tmp_path / "docroot") as url:
        page = browser.new_page()
        problems = watch(page)
        outside: list[str] = []
        page.on("request", lambda r: None if r.url.startswith(url) else outside.append(r.url))
        page.goto(url)
        page.wait_for_function(
            "document.querySelector('#status').textContent.startsWith('Ready')", timeout=LOAD_TIMEOUT
        )
        page.set_input_files("#file", str(image))
        page.wait_for_function("!document.querySelector('#hr-button').disabled", timeout=WORK_TIMEOUT)
        page.click("#hr-button")
        page.wait_for_function(
            "!document.querySelector('#diagram').hidden && document.querySelector('#diagram').complete",
            timeout=WORK_TIMEOUT,
        )
        assert page.evaluate("document.querySelector('#diagram').naturalWidth") > 0  # shown, not blocked
        assert page.is_hidden("#error")
        assert problems == []
        assert outside == []  # everything from the unzipped folder
        page.close()
