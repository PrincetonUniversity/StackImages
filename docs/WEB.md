# Show Clusters on the web

**Live page:** https://princetonuniversity.github.io/StackImages/ShowClusters/

The site root, https://princetonuniversity.github.io/StackImages/, is a StackImages home page
(`web/home.html`) that links to Show Clusters and to the desktop app. GitHub serves a project site
at `<org>.github.io/<repo>/`, so the app lives in a `ShowClusters/` folder of this repo's site rather
than at `…github.io/ShowClusters/`, which would need a renamed or separate repo.

The web page does what `ShowClusters` (in `ass.py`) does: look up a Messier, NGC or IC
object, download its Gaia DR3 stars, draw the simulated star field, and plot the HR
diagram. It runs entirely in the student's browser with
[Pyodide](https://pyodide.org/) (Python compiled to WebAssembly), so there is nothing to
install. The desktop app is unchanged and still ships in `PyWin.tar`.

## How it fits together

```
StackImages/
├── ass.py, ShowClusters.py, …        desktop app (Bob's code, unchanged)
├── web/
│   ├── home.html                     StackImages home page (site root)
│   ├── index.html, style.css         the Show Clusters page (site/ShowClusters/)
│   ├── app.js                        form, Gaia requests, picture sliders, downloads
│   ├── worker.js                     starts Pyodide in a background thread
│   └── showclusters.py               the science: a GUI-free port of ShowClusters
├── scripts/
│   ├── vendor_pyodide.py             downloads + verifies Pyodide (build/vendor/)
│   ├── build_site.py                 assembles the site (_site/)
│   ├── build_snapshots.py            saves Gaia results for the fallback
│   ├── run_web.sh                    build and serve locally
│   └── package_pywin.sh              desktop bundle (PyWin.tar)
├── tests/web/                        unit, snapshot and in-browser tests
├── tests/desktop/                    desktop GUI tests + equivalence test vs ass.py
└── .github/workflows/                ci.yml · pages.yml · snapshots.yml
```

In the browser:

```
index.html ─ app.js ──fetch──▶ GAVO (Gaia DR3)         live star data
               │     ──fetch──▶ data/snapshots/…        saved copy, if GAVO fails
               │     ◀─file──── student's ESA download  manual fallback
               │ postMessage (bytes, inputs)
           worker.js ─ Pyodide ─ showclusters.py        select, render, HR diagram
```

`app.js` does every network request and hands the bytes to Python, so the Python code
never touches the network and is tested natively with the same files.

### What `showclusters.py` ports from `ass.py`

| `ass.py` | Port | Notes |
|---|---|---|
| `getRaDec`, `decimal_to_min_sec`, `min_sec_to_decimal` | `lookup_object`, `format_sexagesimal`, `parse_sexagesimal` | Input is validated first |
| Query in `get_star_info` | `gaia_adql(..., "esa" or "gavo")` | Same ADQL; GAVO's table name and colour columns differ |
| Filtering in `get_star_info` | `select_stars` | Same rules, including float32 storage and the HYG bright stars |
| `show_gaia_image` (2D) | `frame_stars`, `render_field` | Star selection on the original 4000 px grid; picture drawn at 1200 px |
| `slidershow` | `display_base` (Python) + `draw()` in `app.js` | Same sliders and formulas |
| `GUI_makeHRdiag0` | `hr_arrays`, `hr_png` | Same points, sizes, colours, axes |

**Equivalence test.** `tests/desktop/test_web_equivalence.py` imports the real `ass.py` in
the desktop container, replaces only its wx widgets and the Gaia network call (with a
recorded ESA result for M13), and requires the port to match: the query text, the selected
stars in order, the star count and brightest star, the full 4000 px picture, and the HR
diagram in both magnitude modes. It runs in `ci.yml` whenever `ass.py` or
`web/showclusters.py` changes.

**When Bob sends a new `ass.py`:** commit it as usual. If CI's equivalence test fails, the
math changed; update `web/showclusters.py` to match and the test tells you when you have.

### Deliberate differences from the desktop app

| Desktop app | Web page | Why |
|---|---|---|
| 4000 × 4000 picture | 1200 × 1200 (same star positions and counts) | ~240 MB arrays are too heavy for a browser tab |
| "Show Stars" silently resets the magnitude type to Apparent, but the dropdown can still say Absolute | The dropdown visibly switches to Apparent | Same rule, without the hidden state |
| Repeated "HR Diagram" clicks in Absolute mode keep editing the stored arrays (removing another star each time) | Every diagram is computed fresh | Repeatable results |
| Plot title accumulates names (`M13, M44`) within one session | Title is the current object | Clarity |
| 3D, "future" and Color Diagram options | Not included | Not used by HW #7 |
| FOV up to any value | 0.01–5°, magnitude 4.5–21 | Browser memory; the desktop changes behaviour above 5° |
| "Save Figure" | "Download HR diagram" / "Download picture" | Browser equivalent |

## Where the star data comes from

Browsers only let a page read data from another site if that site allows it (CORS).
ESA's Gaia archive does not allow it; the German Astrophysical Virtual Observatory's copy
of Gaia DR3 (`dc.g-vo.org`, table `gaia.dr3lite`) does. So the page tries, in order:

1. **GAVO, live.** The same query the desktop app sends to ESA.
2. **A saved copy**, if GAVO fails or times out and a snapshot covers the requested object,
   field and magnitude limit. The page filters it to the requested field. Snapshots are
   built by `snapshots.yml` (from ESA when it answers, otherwise GAVO) and published as
   `snapshots-YYYYMMDD-HHMM` releases; the Pages build downloads the latest one.
3. **A file the student downloads from the ESA Gaia archive** (CSV, VOTable or FITS). The
   page shows the exact query to paste into the archive's ADQL search.

Every result is labelled with its source on the page and in the HR diagram's corner.

## Will star counts match last year's?

**They should, for the same inputs.** We checked:

| Check | Result |
|---|---|
| Port vs original `ass.py` on identical Gaia data (M13, 0.2°, mag 16) | Identical: 1,369 stars, brightest 11.36, same picture and HR points |
| GAVO vs ESA, M13, 0.2°, mag 16 | Identical: 1,498 rows → 1,369 plotted |
| GAVO vs ESA, M13, 0.2°, mag 18 (HW #7 Problem 1) | Identical: 6,771 rows → 4,588 plotted, brightest 11.36 |
| GAVO vs ESA, M13, 0.3°, mag 16 (the page's defaults) | Identical: 1,609 rows → 1,472 plotted, brightest 10.34 |
| GAVO vs ESA, M44, 1.5°, mag 18 (HW #7 Problem 2) | Not yet compared: ESA's archive returned HTTP 500 on every attempt (2026-11-02). GAVO gives 3,893 plotted |
| Gaia DR3 itself | Frozen since its June 2022 release, so ESA answers the same as last year |

Small differences are still possible in these cases:

- **Different inputs.** The count depends on FOV, maximum magnitude and the centre. Both
  versions query with RA/Dec rounded to whole seconds as shown in the boxes; typing in those
  boxes moves the centre. Last year, a student who skipped Enter got RA = Dec = 0.
- **Stars exactly on an edge.** A star within a rounding error of the picture's edge or the
  query radius can fall either way. The desktop app's numpy 1.26 and the browser's numpy
  2.2 round some float32 arithmetic differently, and saved or uploaded files are cut to the
  field by the page rather than by ESA's server. This affects at most a star or two, and
  none in the fields we tested.
- **GAVO's copy.** It is a re-hosted copy of Gaia DR3 and matched ESA exactly in every
  field we compared, but it is a separate service.
- **HR diagrams after repeated clicks.** Not a count difference: last year's app changed
  its stored data each time "HR Diagram" was clicked in Absolute mode, so a student who
  clicked several times could have a slightly different diagram than a first click gives.

## Workflows

| Workflow | Runs on | Does |
|---|---|---|
| `ci.yml` (desktop) | push/PR to `main`, except web-only files; manual | Desktop GUI tests under Xvfb, equivalence test, packaging tests in `ci/desktop-tests.Dockerfile`; builds `PyWin.tar`; on `main`, publishes release `v<ass.py version>.<n>` when the bundle's app/data files changed (`scripts/next_release.py`) |
| `pages.yml` (web) | push/PR touching the web files; manual | eslint + tsc (Node container), ruff + mypy and all web tests incl. Chromium (`ci/web-tools.Dockerfile`), verified Pyodide download, site build with the latest snapshots; deploys on `main` |
| `snapshots.yml` | manual only | Queries Gaia for ~75 clusters (`scripts/build_snapshots.py`), validates, publishes a `snapshots-*` release, redeploys the page. Resumes from the last release; re-run to fill gaps |

Nothing generated is committed: Pyodide (`build/vendor/`), the site (`_site/`) and
snapshots (`build/snapshots/`) are gitignored and rebuilt or downloaded by the workflows.

## Working on it locally

Everything runs in containers (`compose.yaml`), so only Docker is needed:

```bash
docker-compose build web-tools
docker-compose run --rm web-tools python scripts/vendor_pyodide.py   # once; cached in build/vendor
docker-compose run --rm web-tools                                     # all web tests (unit + Chromium)
docker-compose run --rm lint-js sh -c "npm ci && npm run check"       # eslint + tsc
docker-compose run --rm desktop-tests                                 # desktop + equivalence tests
scripts/run_web.sh                                                    # build and serve on :8000
```

Saved copies for a few clusters, for trying the fallback locally:

```bash
docker-compose run --rm web-tools python scripts/build_snapshots.py --only M13 M44 --source gavo
scripts/run_web.sh --snapshots build/snapshots
```

Without Docker, a `.venv` works too: `uv venv --python 3.13 .venv && uv pip install
--python .venv/bin/python -r requirements-web.lock`, then the same `python …` commands.

## Limits and safeguards

- Inputs are validated before any query: object names must look like `m13`, `ngc5466` or
  `ic4665`; RA/Dec, FOV and magnitude must be in range. Queries are built from numbers only.
- Uploads: at most 50 MB and 300,000 stars, with the required Gaia columns; anything else is
  rejected with a message naming the problem. A file centred on a different object is
  rejected or flagged.
- The only outside service the page contacts is GAVO; the browser tests enforce this.
- Pyodide and every wheel are checked against pinned SHA-256 hashes at build time.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| "The live Gaia data source is unavailable…" | GAVO is down or slow. Saved copies cover common clusters; otherwise use an ESA file (bottom of the page). |
| Page stuck on "Loading Python…" | First visit downloads ~35 MB. On failure, "Reload tool" appears with details. |
| "The Python worker stopped" | Browser ran out of memory; use a smaller field or magnitude limit. |
| `pages.yml` warns "No snapshots-* release yet" | Run `snapshots.yml` once (Actions → Snapshots → Run workflow). |
| `ci.yml` fails in `test_web_equivalence.py` | `ass.py` and the port disagree; see "When Bob sends a new `ass.py`". |
