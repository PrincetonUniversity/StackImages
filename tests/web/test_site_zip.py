"""The downloadable copy of the site (scripts/zip_site.py), served the way a personal web
site serves it: from a subdirectory, with nosniff and a strict Content-Security-Policy."""

from __future__ import annotations

import contextlib
import functools
import http.server
import json
import os
import re
import shutil
import socketserver
import sys
import threading
import zipfile
from collections.abc import Iterator
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).parent))

from test_browser import GAVO, LOAD_TIMEOUT, enter_m13, gavo_ok  # noqa: E402

import build_site  # noqa: E402
import zip_site  # noqa: E402

# What vanderbei.princeton.edu sends (captured 2026-11-03): nosniff, its CSP, and its
# content types; a type it has no entry for comes back as application/octet-stream.
HOST = json.loads((FIXTURES / "personal_site_headers.json").read_text())


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not (ROOT / "build" / "vendor" / "pyodide.js").is_file():
        pytest.fail("build/vendor is missing: run scripts/vendor_pyodide.py first")
    site = tmp_path_factory.mktemp("zip") / "site"
    build_site.build(site, ROOT / "build" / "vendor", None)
    return site


@pytest.fixture(scope="module")
def archive(site: Path) -> Path:
    zip_site.write_zip(site, site / zip_site.NAME)  # where pages.yml puts it
    return site / zip_site.NAME


def test_zip_holds_the_site_in_one_folder(site: Path, archive: Path) -> None:
    names = zipfile.ZipFile(archive).namelist()
    assert all(n.startswith("StackImages/") for n in names)
    expected = {f"StackImages/{p.relative_to(site).as_posix()}" for p in site.rglob("*") if p.is_file()}
    expected -= {"StackImages/.nojekyll", f"StackImages/{zip_site.NAME}"}
    assert set(names) == expected | {"StackImages/HOSTING.txt"}
    assert "StackImages/ShowClusters/vendor/pyodide.asm.wasm" in names
    assert "StackImages/ShowClusters/data/snapshots/index.json" in names
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert z.read("StackImages/ShowClusters/app.js") == (site / "ShowClusters" / "app.js").read_bytes()


def test_zip_is_reproducible(site: Path, archive: Path, tmp_path: Path) -> None:
    again = tmp_path / "again.zip"
    zip_site.write_zip(site, again)
    assert again.read_bytes() == archive.read_bytes()
    assert not list(tmp_path.glob("*.partial"))


def test_fingerprint_ignores_only_the_build_time(site: Path, archive: Path, tmp_path: Path) -> None:
    """A redeploy of the same site must not make a new release; any real change must."""
    copy = tmp_path / "site"
    shutil.copytree(site, copy)
    manifest = copy / "ShowClusters" / "manifest.json"
    record = json.loads(manifest.read_text())
    record["built"] = "2099-01-01T00:00:00+00:00"
    manifest.write_text(json.dumps(record, indent=1) + "\n")
    rebuilt = tmp_path / "rebuilt.zip"
    zip_site.write_zip(copy, rebuilt)
    assert rebuilt.read_bytes() != archive.read_bytes()
    assert zip_site.fingerprint(rebuilt) == zip_site.fingerprint(archive)

    (copy / "ShowClusters" / "app.js").write_text("// changed\n")
    zip_site.write_zip(copy, rebuilt)
    assert zip_site.fingerprint(rebuilt) != zip_site.fingerprint(archive)


def test_fingerprint_ignores_the_version_stamp(site: Path, archive: Path, tmp_path: Path) -> None:
    """CI fingerprints the site before deciding its version, then stamps it: the stamp must not
    count as a change (or every deploy would be a new version)."""
    import web_version

    copy = tmp_path / "site"
    shutil.copytree(site, copy)
    web_version.stamp(copy, "1.7", "2026-11-07")
    stamped = tmp_path / "stamped.zip"
    zip_site.write_zip(copy, stamped)
    assert stamped.read_bytes() != archive.read_bytes()
    assert zip_site.fingerprint(stamped) == zip_site.fingerprint(archive)
    with zipfile.ZipFile(stamped) as z:
        for name in (
            "StackImages/version.js",
            "StackImages/ShowClusters/version.js",
            "StackImages/HRDiagram/version.js",
        ):
            assert 'version: "1.7"' in z.read(name).decode()


def test_fingerprint_command(archive: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert zip_site.main(["--fingerprint", str(archive)]) == 0
    assert capsys.readouterr().out.strip() == zip_site.fingerprint(archive)
    assert len(zip_site.fingerprint(archive)) == 64


def test_every_source_map_named_is_shipped(site: Path) -> None:
    """Browsers' developer tools fetch these; a missing one shows as a 404 in the console."""
    scripts = list((site / "ShowClusters").rglob("*.js"))
    assert scripts
    for script in scripts:
        for name in re.findall(r"//# sourceMappingURL=(\S+)\s*$", script.read_text(errors="replace"), re.M):
            assert (script.parent / name).is_file(), f"{script.name} names {name}, which is not shipped"


def test_pages_use_only_relative_links(site: Path) -> None:
    for page in ("index.html", "ShowClusters/index.html"):
        text = (site / page).read_text()
        for absolute in ('href="/', 'src="/', 'href="https://princetonuniversity.github.io'):
            assert absolute not in text, (page, absolute)


def test_dry_run_writes_nothing(site: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "x.zip"
    assert zip_site.main(["--site", str(site), "--out", str(out), "--dry-run"]) == 0
    assert not out.exists()
    assert "StackImages/HOSTING.txt" in capsys.readouterr().out


def test_refuses_a_directory_that_is_not_a_site(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="not a built site"):
        zip_site.write_zip(tmp_path, tmp_path / "x.zip")


# --- the unzipped copy in a browser, on a server like Bob's ------------------------------


def handler_for(directory: Path, types: dict[str, str]) -> type[http.server.SimpleHTTPRequestHandler]:
    class Host(http.server.SimpleHTTPRequestHandler):
        def guess_type(self, path: str | os.PathLike[str]) -> str:
            return types.get(Path(path).suffix.lower(), "application/octet-stream")

        def end_headers(self) -> None:
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy", HOST["content-security-policy"])
            super().end_headers()

        def log_message(self, *args: object) -> None:
            pass

    return functools.partial(Host, directory=str(directory))  # type: ignore[return-value]


@contextlib.contextmanager
def personal_site(archive: Path, root: Path, types: dict[str, str]) -> Iterator[str]:
    """Unzip into ~rvdb/ of a fresh docroot and serve it; yields the Show Clusters URL."""
    with zipfile.ZipFile(archive) as z:
        z.extractall(root / "~rvdb")
    with socketserver.ThreadingTCPServer(("127.0.0.1", 0), handler_for(root, types)) as server:
        threading.Thread(target=server.serve_forever, daemon=True).start()
        try:
            yield f"http://127.0.0.1:{server.server_address[1]}/~rvdb/StackImages/ShowClusters/"
        finally:
            server.shutdown()


def watch(page) -> list[str]:  # type: ignore[no-untyped-def]
    problems: list[str] = []
    page.on("pageerror", lambda e: problems.append(f"error: {e}"))
    page.on("console", lambda m: problems.append(m.text) if "Content Security Policy" in m.text else None)
    return problems


def test_unzipped_copy_works_under_a_strict_host(browser, archive: Path, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    """Bob's site as it is, VM or Blob: no entries for .mjs, .wasm, .npy or .whl, so they all
    come back as application/octet-stream, with nosniff. No server change may be needed."""
    with personal_site(archive, tmp_path, HOST["types"]) as url:
        page = browser.new_page()
        problems = watch(page)
        page.route(f"{GAVO}*", gavo_ok([]))
        page.goto(url)
        page.wait_for_function(
            "document.getElementById('status').textContent.startsWith('Ready')", timeout=LOAD_TIMEOUT
        )
        enter_m13(page)
        page.click("#show-stars")
        page.wait_for_function("!document.getElementById('stars-content').hidden", timeout=60_000)
        assert page.text_content("#count") == "1,369"
        page.click("#hr-button")
        page.wait_for_function(
            "!document.getElementById('hr').hidden && document.getElementById('hr').complete", timeout=60_000
        )
        assert page.evaluate("document.getElementById('hr').naturalWidth") == 1000  # shown, not blocked
        page.click("a[href='../']")
        page.wait_for_selector("#open-clusters")  # the home page in the same folder
        assert problems == []
        page.close()
