# StackImages on the web

**Live pages:**

- Show Clusters: https://princetonuniversity.github.io/StackImages/ShowClusters/
- HR Diagram of Your Own Image: https://princetonuniversity.github.io/StackImages/HRDiagram/
  (see [HR diagram of your own image](#hr-diagram-of-your-own-image) below)

The site root, https://princetonuniversity.github.io/StackImages/, is a StackImages home page
(`web/home.html`) that links to both pages and to the desktop app. GitHub serves a project site
at `<org>.github.io/<repo>/`, so the app lives in a `ShowClusters/` folder of this repo's site rather
than at `…github.io/ShowClusters/`, which would need a renamed or separate repo.

Show Clusters does what `ShowClusters` (in `ass.py`) does: look up a Messier, NGC or IC
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
| "Show Stars" silently resets the magnitude type to Apparent (to absolute magnitudes above 5°), but the dropdown can still say otherwise | The dropdown visibly switches to Apparent, or to Absolute above 5° | Same rule, without the hidden state. Above 5° `ass.py` sets an internal `'Abs'` that plots absolute magnitudes but skips Absolute's removal of the brightest point; choosing Absolute in either app gives the same diagram |
| Repeated "HR Diagram" clicks in Absolute mode keep editing the stored arrays (removing another star each time) | Every diagram is computed fresh | Repeatable results |
| Plot title accumulates names (`M13, M44`) within one session | Title is the current object | Clarity |
| 3D, "future" and Color Diagram options | Not included | Not used by HW #7 |
| FOV up to any value | 0.01–90°, magnitude 4.5–21; live data for fields wider than ~1.3° only up to magnitude 8.5 (deeper wide fields need an ESA file) | Beyond ~90° the tangent-plane picture breaks down; GAVO can't send wide, deep fields in time (below) |
| "Save Figure" | "Download HR diagram" / "Download picture" | Browser equivalent |
| Picture in its own window, sliders under it | Picture sized to the screen with the sliders beside it (or below on narrow screens), plus a **full screen** icon (⤢) on the picture's corner (black surround, dark controls; Esc or ✕ to exit) and **Reset sliders** set apart below the downloads, which asks before discarding moved sliders | Picture and sliders are always visible together; Reset restores the starting stretch without being hit by accident; full screen gives the picture the whole screen |

## Where the star data comes from

Browsers only let a page read data from another site if that site allows it (CORS).
ESA's Gaia archive does not allow it; the German Astrophysical Virtual Observatory's copy
of Gaia DR3 (`dc.g-vo.org`, table `gaia.dr3lite`) does. So the page tries, in order:

1. **GAVO, live.** The same query the desktop app sends to ESA. Large or crowded fields are
   fetched in pieces (see below).
2. **A saved copy**, if GAVO fails or times out and a snapshot covers the requested object,
   field and magnitude limit. The page filters it to the requested field. Snapshots are
   built by `snapshots.yml` (from GAVO in pieces, or ESA if GAVO fails) and published as
   `snapshots-YYYYMMDD-HHMM` releases; the Pages build downloads the latest one.
3. **A file the student downloads from the ESA Gaia archive** (CSV, VOTable or FITS). The
   page shows the exact query to paste into the archive's ADQL search.

Every result is labelled with its source on the page and in the HR diagram's corner.

### Large and crowded fields

GAVO's quick endpoint stops any query after about 20 seconds ("Query timed out"). Measured on
2026-11-02, a circle up to 1° radius answers in 1–2 s in an ordinary field, but 1.2° times out
(M44 with FOV 2 did), and in the dense Milky Way bulge (M7, M22) only about 0.2° answers in time.
Many of those dense-field timeouts weren't about size. With a plain
`phot_g_mean_mag < 16` condition, GAVO's database can choose its magnitude index, which
in a dense field means scanning every bright star in the sky. Measured on 2026-11-03,
M23 at FOV 0.3 timed out that way, while the same circle without the magnitude
condition returned 37,793 stars in 11 s. The page therefore writes the condition as
`phot_g_mean_mag + 0 < 16`. That selects the same stars but keeps the database on the
position index:

| Field | Before | With `+ 0` |
|---|---|---|
| M23, FOV 0.3, G < 16 | timed out (20 s) | 1,319 stars, ~2–15 s |
| M22, FOV 0.5, G < 18 | timed out | 44,252 stars, ~25 s |
| NGC 6397, FOV 0.5, G < 17 | 9.0 s | 1.9 s |

M13, M44 and M45 were unchanged, and every result both forms returned was identical.
Only GAVO's query changed; the ESA query stays exactly as `ass.py` writes it.

So the page:

1. Sends a field whose query circle is at most 0.9° radius (FOV up to ~1.27°) as **one query**,
   exactly as before.
2. Covers a larger field with **pieces** of at most 0.9° radius (FOV 2 → 9 pieces), fetched four
   at a time. That limit applies to every piece, including the smaller ones made by splitting. Up
   to 2026-11-03, the pieces of a split came one after another, which made crowded fields very slow.
3. **Splits a piece that times out** into four smaller ones, down to 0.05°; if even those time
   out, it says the field has too many stars and suggests a smaller FOV or magnitude limit.
   A field that would need more than 64 pieces (a 5° field needs 36) isn't attempted (see below).
4. **Reassembles** the pieces into exactly the one-query result: stars seen in two overlapping
   pieces are kept once, the result is trimmed to the query's circle and magnitude limit (using
   the centre and radius exactly as written in the query), and stars are put in a fixed order.

`tests/web/test_tiles.py` checks that this gives the **same end result** as a single query,
byte for byte: star count, brightest star, every star's values, the picture's pixels, the HR
diagram's points and the downloaded PNG, for one and for two levels of splitting. It also
checks that the pieces cover the field (including near the poles and across RA 0h).
`test_browser.py` checks the same through the real page, including that the four pieces of a
split are requested together.

While the data downloads, the page shows a progress bar (parts done, for a field in pieces)
and a clock: "*N* s so far. Please keep this tab open; crowded parts of the Milky Way can take a
minute or two." Only the status line (which changes when a part arrives) is read by screen readers.
When pieces have been split, the status line says so. The page gives up after 3 minutes and
falls back to a saved copy or an ESA file. M23 at FOV 2 took 139 s on 2026-11-03.

### Wide fields (constellations)

Pieces don't scale to constellation-sized fields: Orion at FOV 60 would need several thousand.
But with a bright limit, the opposite trick works: GAVO's magnitude index (the plain condition,
without `+ 0`) starts from the brightest stars in the sky, and answers the whole circle at once.
Measured on 2026-11-07, centred on Orion:

| Field | Magnitude index, one query | `+ 0`, one query |
|---|---|---|
| FOV 60, G < 7 (HW Problem 1) | 3,285 stars, 2 s | timed out |
| FOV 60, G < 8 | 9,503 stars, 3 s | |
| FOV 60, G < 8.5 | 16,274 stars, 2 s | |
| FOV 60, G < 9 | timed out (28 s); 22 s near the Galactic centre | |
| FOV 45, G < 10 | timed out | |
| FOV 30, G < 12 | timed out | timed out even at FOV ~7 |

So a field wider than one piece with Maximum Magnitude **8.5 or less** is sent as one
magnitude-index query (not split if it times out). A wider field with a deeper limit that would
need more than 64 pieces is not sent at all: the page says at once that wide fields come from
the live source only up to magnitude 8.5, and opens the ESA download (the archive's
asynchronous jobs have no 20 s limit). ESA counted 44,639 stars for FOV 45, G < 10 and 123,711
for FOV 30, G < 12, centred on Orion; as CSV with the page's query columns those are about 8 and
21 MB, well inside the 50 MB upload limit (which matches the 300,000-star limit: ~50 MB of CSV,
~30 MB VOTable, ~23 MB FITS).

## Downloads

| Button | File | Contents |
|---|---|---|
| Download picture (PNG) | `M13.png` | The star field as currently stretched by the sliders; PNG text fields record the field, the exact slider positions, the data source and the credits |
| Download star data (CSV) | `M13_FOV0.2_G16_GAVO_2026-11-03_stars.csv` | The field's Gaia rows exactly as used, after any reassembly and trimming, in ESA's `gaiadr3.gaia_source` columns (`ra, dec, parallax, phot_g_mean_mag, bp_rp, pmra, pmdec`) |
| Download HR diagram (PNG) | `M13_HR_Absolute.png` | The diagram shown |
| Download HR points (CSV) | `M13_HR_Absolute_points.csv` | One row per star in the picture: catalog (Gaia DR3 or HYG), position, parallax, apparent and absolute magnitude, color index, the plotted x and y, and `drawn` |

- **Re-uploadable star data.** The upload box reads the star-data file back unchanged and
  reproduces the same count, picture and HR points (tested), so a student can submit the exact
  data behind a diagram, or redo the work offline.
- **Picture settings travel with the data.** The star-data file records the sliders' exact positions (`# Display sliders: colsat=…; gamma=…; mini=…; maxi=… (shown as …)`). Uploading it again sets the sliders back, so the picture matches pixel for pixel (tested); Reset still returns to the standard starting positions. Malformed or edited lines are ignored.
- **The field travels with the data.** All four downloads record the five form inputs under the form's own labels (`Field: Object Catalog Name=M13; Right Ascension (RA)=…; Declination (Dec)=…; Field of View (FOV)=0.2; Maximum Magnitude=16`; separate text fields in the PNGs). Uploading a star-data file fills in the form from it, so the file is shown as saved whatever the form held (e.g. a grader who has typed nothing). Invalid records are ignored.
- **Provenance and credit in the file.** Both CSVs start with `#` comment lines: the field, the
  data source and date, and the credits below. The upload reader skips `#` lines.
- **`drawn`** comes from the diagram's own point sizes, so the drawn rows are exactly the diagram's
  points. As in the desktop app, the first star, the faintest star (size 0) and, in Absolute mode,
  the brightest star are not drawn.
- **HYG rows** use V and B−V (Gaia rows use G and BP−RP), as the desktop app mixes them; their
  parallax is blank.
- **HR points follow the diagram on screen**, not the menu: changing Apparent/Absolute without
  redrawing still downloads the points of the diagram shown.

## Data credits

The page footer, the text beside the download buttons, and every CSV carry these credits (also in
the README's "Data and credits"):

- **Gaia DR3:** ESA/Gaia/DPAC, [CC BY-NC 3.0 IGO](https://www.cosmos.esa.int/web/gaia-users/license),
  with ESA's acknowledgement sentence and the two papers ESA asks to be cited (Gaia Collaboration
  2016, A&A 595, A1; 2023, A&A 674, A1).
- **HYG database v3:** David Nash, CC BY-SA 2.5. The site redistributes a subset of it
  (`data/bright_stars.npy`, built from `csv/hygdata_v3.csv`), which remains under CC BY-SA 2.5.
- **GAVO** for serving the live data.

Classroom use fits the non-commercial Gaia licence. The credit text lives in
`web/showclusters.py` (`GAIA_CREDIT`, `HYG_CREDIT`, `GAVO_CREDIT`) and `web/index.html`.

## Will star counts match last year's?

**They should, for the same inputs.** We checked:

| Check | Result |
|---|---|
| Port vs original `ass.py` on identical Gaia data (M13, 0.2°, mag 16) | Identical: 1,369 stars, brightest 11.36, same picture and HR points |
| Port vs original `ass.py` on identical Gaia data (Orion, 60°, mag 7, centred on Alnilam) | Identical stars, picture and HR points |
| GAVO vs ESA, Orion, 60°, mag 7 | Identical: 3,259 rows |
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

## HR diagram of your own image

**Live page:** https://princetonuniversity.github.io/StackImages/HRDiagram/

The page does what these desktop buttons do, in the same order as the StackImages window:

| Desktop (StackImages window) | `ass.py` | Page | Port (`web/hrimage.py`) |
|---|---|---|---|
| Open an RGB fits file: **Select** | `get_RGB_image` → `show_stacked_image` → `slidershow` | 1. choose a FITS file; 2. the picture and its sliders | `read_rgb_fits`, `slidershow_input` + `showclusters.display_base`, `web/stretch.js` |
| Brightest star (mag), Faintest star (mag) | `mag0`, `mag1` (empty → 10 and 20) | the same two boxes | `parse_magnitude` |
| **HR diagram** | `GUI_HRdiag` → `makeHRdiag_new(..., shape='auto', fignum=1, thresh=1e-5)` → `findstars`, `fwhm` | **HR diagram** | `measure` (→ `findstars`, `fwhm`), `plot_points`, `hr_figure` |
| **Color Diagram** | `GUI_ColorDiag` → the same, `fignum=6` | **Color Diagram** | `color_figure` |

Files: `web/hrdiagram/` (the page: `index.html`, `app.js`, `worker.js`), `web/hrimage.py`, and the
shared `web/stretch.js` (slidershow's stretch, also used by Show Clusters), `web/style.css` and
`web/showclusters.py` (`display_base`). The site build puts them in `HRDiagram/`. The page uses
Show Clusters' copy of Pyodide (`../ShowClusters/vendor/`), plus scipy for `convolve2d`.

**Equivalence test.** `tests/desktop/test_hr_equivalence.py` runs Bob's unmodified
`get_RGB_image`, `GUI_HRdiag` and `GUI_ColorDiag` in the desktop container, replacing only the wx
widgets, the file dialog and the plot window, on synthetic star clusters (`tests/web/hrfixtures.py`:
stars of several widths, hot pixels, extreme colors). The port must give the **identical** image
array, name and picture input, and identical plotted points (positions, sizes, face and edge
colors), axis limits, titles and labels, for both diagrams and several magnitude settings. In
the browser, `tests/web/test_hr_browser.py` checks that Pyodide's plotted points match the native
port's.

### What the port does differently (and why the points are still the same)

The arithmetic follows `ass.py` statement by statement, with the same dtypes. Only memory and
skipped work differ. Run unchanged, `makeHRdiag_new` peaks at about 5.5 GB for a 19-megapixel
image, above a browser tab's 4 GB limit; the port peaks at about 1.7 GB:

- `findstars` computes one filter at a time, applies its test, and frees it, instead of holding
  18 full-size arrays. The eight neighbour filters (M1–M8) are only tested for ≥ 0, so they are
  written as the equivalent neighbour comparisons. Filters it never uses (M9–M12, M21, M22) are
  skipped.
- FWHM is measured only at candidate stars, in the same order, instead of looping over every pixel.
- The `shape='auto'` cluster mask is not computed: `ass.py` computes it but never applies it
  (`k = k&mask` is commented out).
- Detection and photometry don't depend on the magnitude boxes, so they run once per image;
  changing the magnitudes, or switching to the Color Diagram, redraws from the stored
  measurements. The desktop app recomputes everything on each click, with the same result.

### Deliberate differences from the desktop app

| Desktop app | Web page | Why |
|---|---|---|
| Picture at full resolution | Averaged down to at most 1600 px on its longest side; the slider ranges and starting positions come from the full image | Browser memory and drawing speed; the diagrams always use every pixel |
| Pixel info, star info, FWHM, Astrometric Analysis, Deconvolve, Unsharp Mask, stacking | Not included | Not needed for the HR and color diagrams |
| Writes `txt/color_mag.txt` | **Download points (CSV)**: every plotted star's position, fluxes, color and magnitude | Browser equivalent, with more detail |
| "Save" in the plot window | **Download diagram (PNG)** | Browser equivalent |
| No checks on the file | Files over 250 MB or 20 megapixels, non-FITS files, images with NaN/infinite pixels and impossible magnitudes are refused with an explanation | `ass.py` would crash or plot nothing |

### Kept exactly as in `ass.py`

These look like slips in `ass.py`. They are reproduced on purpose so the page matches the desktop
app, and are filed as issues for Bob. If one is fixed in `ass.py`, change `web/hrimage.py` and
`tests/desktop/test_hr_equivalence.py` with it.

- The color filter `abs(log(G/B) > 2)` takes the absolute value of the comparison, so it drops
  stars with G/B > e² but not G/B < e⁻² (the R/G test is a true `|log(R/G)| > 2`).
  ([#21](https://github.com/PrincetonUniversity/StackImages/issues/21))
- The `shape='auto'` cluster mask is computed but not applied (see above).
  ([#22](https://github.com/PrincetonUniversity/StackImages/issues/22))
- The Color Diagram's upper axis limits are padded using the already-padded lower limits.
  ([#24](https://github.com/PrincetonUniversity/StackImages/issues/24))

Two more issues found while porting don't affect the page: `fig5axes` is never set when `GUI` is
False ([#23](https://github.com/PrincetonUniversity/StackImages/issues/23); the equivalence test
supplies the axes), and `ass.py`'s star imports only work at the pinned numpy and scipy versions
([#25](https://github.com/PrincetonUniversity/StackImages/issues/25)).

### Limits and timing

Measured in Chromium on an Apple-silicon Mac with a synthetic 5,345 × 3,551 RGB image (217 MB,
like a cropped M2 stack): opening it took 5 s, the first HR diagram 32 s, and each later diagram
0.3 s. Most of the first diagram's time is `fwhm`, called once per candidate star exactly as
`ass.py` does. The page refuses files over 250 MB or 20 megapixels; a smaller crop around the
cluster is faster.

## Workflows

| Workflow | Runs on | Does |
|---|---|---|
| `ci.yml` (desktop) | push/PR to `main`, except web-only files and Markdown; manual | Desktop GUI tests under Xvfb, equivalence test, packaging tests in `ci/desktop-tests.Dockerfile`; builds `PyWin.tar`; on `main`, publishes release `v<ass.py version>.<n>` when the bundle's app/data files changed (`scripts/next_release.py`) |
| `pages.yml` (web) | push/PR touching the web files (not Markdown); manual | eslint + tsc (Node container), ruff + mypy and all web tests incl. Chromium (`ci/web-tools.Dockerfile`), verified Pyodide download, site build with the latest snapshots, `StackImages-web.zip` of the site and `HRDiagram-web.zip` of the HR page; deploys on `main`, then publishes each zip whose contents changed as a release (`web-*`, `hrdiagram-*`) |
| `snapshots.yml` | manual only | Fetches ~75 clusters (`scripts/build_snapshots.py`: GAVO in pieces, ESA as fallback, lower magnitude limit if both fail), validates, publishes a `snapshots-*` release, redeploys the page. Resumes from the last release; re-run to fill gaps |

Documentation-only changes (`README.md`, anything under `docs/`, any other `.md` file) start no
workflow: there is nothing to test, deploy or release. The README inside `PyWin.tar` is refreshed
by the next desktop release, and README/LICENSE edits alone never make one. To redeploy or release
anyway, run the workflow manually (Actions → Run workflow).

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

## Hosting a copy elsewhere

Every deploy also publishes the whole site as one zip:
https://princetonuniversity.github.io/StackImages/StackImages-web.zip
(`scripts/zip_site.py`, run by `pages.yml` after the site build). Each version is also
kept as a `web-YYYYMMDD-HHMM` release, never marked Latest. A new release is made only
when the site's contents change: `zip_site.py --fingerprint` hashes the zip without the
build time. Unzipped, it is a
`StackImages/` folder: the home page, `ShowClusters/` with Pyodide and the saved
snapshots, `HRDiagram/`, and `HOSTING.txt` with these instructions. Copy the folder unchanged into
any directory of a static web site. All links are relative, so it works from a
subdirectory such as `/~rvdb/StackImages/`.

No server settings are needed. Many hosts send `application/octet-stream` for
extensions they don't know. Browsers refuse that type for JavaScript modules
(`.mjs`), so the page loads Pyodide differently:

- **The worker is a classic one.** It loads `pyodide.js` with `importScripts`, and
  every host types `.js` correctly. `pyodide.mjs` isn't shipped.
- **`.wasm` compiles from its bytes** when it doesn't arrive as `application/wasm`.
  Pyodide would otherwise stall on it.
- **The data files are fetched as plain bytes**, so `.whl`, `.npy`, `.csv` and `.json`
  work with any content type.

`test_site_zip.py` serves the unzipped copy the way vanderbei.princeton.edu does,
using its headers saved in `tests/web/fixtures/personal_site_headers.json`:

- `nosniff` and its Content-Security-Policy;
- the content types from its upload map;
- `application/octet-stream` for anything unmapped.

Its CSP permits `data:` images but not `blob:` ones, so the HR diagram is shown
through a `data:` URL.

On vanderbei.princeton.edu this holds both on the current nginx VM and on the coming
Azure Blob storage. The new extensions (`.wasm`, `.whl`, `.npy`) only need the
vanderbei-deploy step for any new file type: `tools/gen-mime-map.sh` records whatever
nginx serves for them, so the autosync to Blob and the cutover parity check include
them.

## Limits and safeguards

- Inputs are validated before any query: object names must look like `m13`, `ngc5466` or
  `ic4665`; RA/Dec, FOV and magnitude must be in range. Queries are built from numbers only.
- Uploads: at most 50 MB and 300,000 stars, with the required Gaia columns; anything else is
  rejected with a message naming the problem. A file centred on a different object is
  rejected or flagged.
- The only outside service Show Clusters contacts is GAVO; the browser tests enforce this.
- HR Diagram of Your Own Image contacts no outside service, and the student's image never
  leaves the browser. It refuses files over 250 MB or 20 megapixels, non-FITS files, images with
  NaN or infinite pixels, and impossible magnitudes, with a message naming the problem.
- Pyodide and every wheel are checked against pinned SHA-256 hashes at build time.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| "The live Gaia data source is unavailable…" | GAVO is down or slow. Saved copies cover common clusters; otherwise use an ESA file (bottom of the page). |
| "This field has too many stars for the live data source…" | A wide or Milky-Way-dense field even GAVO's smallest pieces can't answer in time. Lower the FOV or Maximum Magnitude, or use an ESA file. |
| "Wide fields come from the live data source only up to Maximum Magnitude 8.5…" | A wide field with a deep limit (e.g. FOV 45, magnitude 10). Download it from the ESA Gaia archive with the query shown, and choose the file at the bottom of the page. |
| Page stuck on "Loading Python…" | First visit downloads ~35 MB. On failure, "Reload tool" appears with details. |
| "The Python worker stopped" | Browser ran out of memory; use a smaller field or magnitude limit. |
| `pages.yml` warns "No snapshots-* release yet" | Run `snapshots.yml` once (Actions → Snapshots → Run workflow). |
| `ci.yml` fails in `test_web_equivalence.py` | `ass.py` and the port disagree; see "When Bob sends a new `ass.py`". |
