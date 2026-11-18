# StackImages on the web

- Home: https://princetonuniversity.github.io/StackImages/ (`web/home.html`)
- Show Clusters: https://princetonuniversity.github.io/StackImages/ShowClusters/
- HR Diagram of Your Own Image: https://princetonuniversity.github.io/StackImages/HRDiagram/

Both tools run Bob's desktop code, ported, in the browser with [Pyodide](https://pyodide.org/).
The desktop app is unchanged and still ships in `PyWin.tar`. GitHub serves project sites at
`<org>.github.io/<repo>/`, hence the `StackImages/` prefix.

## Layout

```
StackImages/
├── ass.py, ShowClusters.py, …        desktop app (unchanged)
├── web/
│   ├── home.html                     site root
│   ├── index.html, app.js, worker.js Show Clusters page (site/ShowClusters/)
│   ├── showclusters.py               GUI-free port of ShowClusters
│   ├── hrdiagram/, hrimage.py        HR Diagram page and its port
│   └── stretch.js, style.css         shared by both pages
├── scripts/                          vendor_pyodide.py, build_site.py, build_snapshots.py,
│                                     zip_site.py, web_version.py, run_web.sh, package_pywin.sh
├── tests/web/                        unit, snapshot and browser tests
├── tests/desktop/                    desktop GUI tests + equivalence tests vs ass.py
└── .github/workflows/                ci.yml · pages.yml · snapshots.yml
```

```
index.html ─ app.js ──fetch──▶ GAVO (Gaia DR3)         live star data
               │     ──fetch──▶ data/snapshots/…        saved copy, if GAVO fails
               │     ◀─file──── student's ESA download  manual fallback
               │ postMessage (bytes, inputs)
           worker.js ─ Pyodide ─ showclusters.py        select, render, HR diagram
```

`app.js` makes every network request and hands bytes to Python, so the Python never touches the
network and is tested natively on the same files.

## Show Clusters

### Port of `ass.py`

| `ass.py` | `showclusters.py` | Notes |
|---|---|---|
| `getRaDec`, `decimal_to_min_sec`, `min_sec_to_decimal` | `lookup_object`, `format_sexagesimal`, `parse_sexagesimal` | Input validated first |
| Query in `get_star_info` | `gaia_adql(..., "esa" or "gavo")` | Same ADQL; GAVO's table and colour columns differ |
| Filtering in `get_star_info` | `select_stars` | Same rules, incl. float32 and HYG bright stars |
| `show_gaia_image` (2D) | `frame_stars`, `render_field` | Selection on the 4000 px grid; drawn at 1200 px |
| `slidershow` | `display_base` + `draw()` in `app.js` | Same sliders and formulas |
| `GUI_makeHRdiag0` | `hr_arrays`, `hr_png` | Same points, sizes, colours, axes |

**Equivalence test.** `tests/desktop/test_web_equivalence.py` runs the real `ass.py` in the
desktop container, replacing only wx and the Gaia call (a recorded ESA result for M13), and
requires identical query text, selected stars in order, count, brightest star, 4000 px picture
and HR diagram in both modes. It runs in `ci.yml` when `ass.py` or `web/showclusters.py` changes.

**When Bob sends a new `ass.py`:** commit it. If the equivalence test fails, the math changed;
update `web/showclusters.py` until it passes.

### Differences from the desktop app

| Desktop | Web | Why |
|---|---|---|
| 4000 × 4000 picture | 1200 × 1200, same stars | ~240 MB is too much for a tab |
| "Show Stars" silently resets the magnitude type (Apparent; Absolute above 5°) | The dropdown visibly switches | Same rule, no hidden state. Above 5° `ass.py`'s internal `'Abs'` skips Absolute's removal of the brightest point; choosing Absolute gives the same diagram in both |
| Repeated Absolute "HR Diagram" clicks remove another star each time | Each diagram computed fresh | Repeatable |
| Title accumulates names (`M13, M44`) | Current object only | Clarity |
| 3D, "future", Color Diagram | Not included | Not used by HW #7 |
| Any FOV | 0.01–90°, magnitude 4.5–21 | Tangent-plane picture fails beyond ~90° |
| "Save Figure" | Download buttons | Browser equivalent |
| Picture in its own window | Picture beside the sliders (below on narrow screens), ⤢ full screen (Esc or ✕), **Reset sliders** below the downloads with a confirmation | Picture and sliders always visible together |

### Star data sources

ESA's Gaia archive doesn't allow cross-site reads (CORS); GAVO's copy of Gaia DR3
(`dc.g-vo.org`, `gaia.dr3lite`) does. The page tries, in order:

1. **GAVO, live**: the desktop app's query, in pieces or as an async job when large (below).
2. **A saved copy**, if GAVO fails and a snapshot covers the object, field and magnitude.
   `snapshots.yml` builds them (GAVO, or ESA as fallback) as `snapshots-YYYYMMDD-HHMM` releases;
   the Pages build includes the latest.
3. **An ESA archive download** (CSV, VOTable or FITS); the page shows the ADQL to paste.

The source is labelled on the page and in the HR diagram's corner.

### Large and crowded fields

GAVO's quick endpoint (`tap/sync`) stops queries after ~20 s. A 1° radius answers in 1–2 s in an
ordinary field; in the Milky Way bulge only ~0.2° does.

**`+ 0` trick.** With `phot_g_mean_mag < 16`, GAVO may use its magnitude index and scan every
bright star in the sky. The page writes `phot_g_mean_mag + 0 < 16`: same stars, position index.
Measured 2026-11-03 (results identical in every case; ESA's query is unchanged):

| Field | Plain | With `+ 0` |
|---|---|---|
| M23, FOV 0.3, G < 16 | timed out | 1,319 stars, 2–15 s |
| M22, FOV 0.5, G < 18 | timed out | 44,252 stars, ~25 s |
| NGC 6397, FOV 0.5, G < 17 | 9.0 s | 1.9 s |

**Pieces.** A query circle up to 0.9° radius (FOV ≤ ~1.27°) is one query. Larger fields are
covered by pieces of at most 0.9° radius (FOV 2 → 9), fetched four at a time. A piece that times
out splits into four, down to 0.05°; below that the page says the field has too many stars. The
pieces are reassembled into exactly the one-query result (duplicates dropped, trimmed to the
circle and magnitude, fixed order). `tests/web/test_tiles.py` checks that byte for byte (stars,
picture, HR points, PNG; one and two levels of splitting; poles and RA 0h), and `test_browser.py`
checks it through the page. The page shows progress and a clock, and gives up after 3 minutes
(M23 at FOV 2 took 139 s).

**Wide fields.** Constellation-sized fields would need thousands of pieces. With a bright limit,
the plain magnitude-index query answers the whole circle at once. Measured 2026-11-07 on Orion:

| Field | Magnitude index, one query |
|---|---|
| FOV 60, G < 7 (HW Problem 1) | 3,285 stars, 2 s (`+ 0` times out) |
| FOV 60, G < 8.5 | 16,274 stars, 2 s |
| FOV 60, G < 9 | timed out (22 s near the Galactic centre) |

So a field wider than one piece with Maximum Magnitude ≤ 8.5 is one magnitude-index query.

**Async jobs.** A wider, deeper field (over 64 pieces; a 5° field needs 36) is one job on GAVO's
`tap/async` (no 20 s limit, CORS allowed). The page creates and starts the job, polls every 3 s,
downloads the CSV, and always deletes the job; it gives up after 15 minutes and then offers the
ESA download. Measured 2026-11-07 on Orion:

| Field | Job | Stars | Plotted (Chromium) |
|---|---|---|---|
| FOV 45, G < 10 | 112–224 s | 44,735 (4.3 MB) | 27,503, < 1 s; matches ESA |
| FOV 30, G < 12 | 279 s | 120,593 (11.8 MB) | 70,054, 0.9 s |

**ESA downloads stop at 50,000 stars** unless logged in (marked `OVERFLOW` in VOTable, unmarked in
CSV). Uploads with exactly 50,000 stars or `OVERFLOW` are refused as incomplete.

### Downloads

| Button | Example file | Contents |
|---|---|---|
| Download picture (PNG) | `M13.png` | Field as stretched; PNG text fields hold the field, slider positions, source and credits |
| Download star data (CSV) | `M13_FOV0.2_G16_GAVO_2026-11-03_stars.csv` | The Gaia rows used, in ESA's columns (`ra, dec, parallax, phot_g_mean_mag, bp_rp, pmra, pmdec`) |
| Download HR diagram (PNG) | `M13_HR_Absolute.png` | The diagram shown |
| Download HR points (CSV) | `M13_HR_Absolute_points.csv` | One row per star: catalog, position, parallax, magnitudes, colour, plotted x/y, `drawn` |

- All four record the form inputs (`Field: Object Catalog Name=M13; …`); the CSVs also start with
  `#` lines for the source, date and credits.
- Uploading the star-data file reproduces the count, picture and HR points (tested), restores the
  form and the slider positions (`# Display sliders: …`). Invalid lines are ignored.
- `drawn` follows the diagram's point sizes: as in the desktop app, the first star, the faintest
  (size 0) and, in Absolute mode, the brightest are not drawn.
- HYG rows use V and B−V (Gaia rows G and BP−RP), as the desktop app does; their parallax is blank.
- HR points are those of the diagram shown, even if the Apparent/Absolute menu has changed since.

### Will star counts match last year's?

Yes, for the same inputs:

| Check | Result |
|---|---|
| Port vs `ass.py`, same data: M13 0.2°/16; Orion 60°/7 | Identical stars, picture, HR points |
| GAVO vs ESA: Orion 60°/7; M13 0.2°/16, 0.2°/18 (HW #7 P1), 0.3°/16 | Identical rows |
| GAVO vs ESA: M44 1.5°/18 (HW #7 P2) | Not compared (ESA HTTP 500, 2026-11-02); GAVO gives 3,893 plotted |
| Gaia DR3 | Frozen since June 2022 |

Possible small differences: different inputs (RA/Dec are rounded to whole seconds as shown; last
year, skipping Enter gave RA = Dec = 0); a star or two on an edge (numpy 1.26 vs 2.2 float32
rounding; saved and uploaded files are cut by the page); and repeated Absolute clicks in last
year's app.

## HR diagram of your own image

| Desktop button | `ass.py` | `web/hrimage.py` |
|---|---|---|
| Open an RGB fits file | `get_RGB_image` → `show_stacked_image` → `slidershow` | `read_rgb_fits`, `slidershow_input` + `display_base`, `stretch.js` |
| Brightest / Faintest star (mag) | `mag0`, `mag1` (empty → 10, 20) | `parse_magnitude` |
| HR diagram | `GUI_HRdiag` → `makeHRdiag_new(..., shape='auto', fignum=1, thresh=1e-5)` → `findstars`, `fwhm` | `measure`, `plot_points`, `hr_figure` |
| Color Diagram | `GUI_ColorDiag` (`fignum=6`) | `color_figure` |

The page (`web/hrdiagram/`) uses Show Clusters' Pyodide (`../ShowClusters/vendor/`) plus scipy.

**Equivalence test.** `tests/desktop/test_hr_equivalence.py` runs Bob's unmodified functions,
replacing only wx, the file dialog and the plot window, on synthetic clusters
(`tests/web/hrfixtures.py`). The port must give identical image arrays, plotted points, axes and
labels for both diagrams. `tests/web/test_hr_browser.py` checks Pyodide against the native port.

**Memory.** The arithmetic follows `ass.py` statement by statement. `makeHRdiag_new` peaks at
~5.5 GB on a 19-megapixel image (a tab allows 4 GB); the port peaks at ~1.7 GB because it:

- computes `findstars`' filters one at a time, writes M1–M8 (only tested ≥ 0) as neighbour
  comparisons, and skips the unused M9–M12, M21, M22;
- measures FWHM only at candidate stars;
- skips the `shape='auto'` mask, which `ass.py` computes but never applies;
- detects and measures once per image, so changing magnitudes or diagrams just redraws.

**Differences from the desktop app.** The picture is averaged to ≤ 1600 px (diagrams use every
pixel). Pixel/star info, FWHM, astrometry, deconvolution, unsharp mask and stacking are left out.
**Download points (CSV)** replaces `txt/color_mag.txt`; **Download diagram (PNG)** replaces
"Save". Files over 250 MB or 20 megapixels, non-FITS files, NaN/infinite pixels and impossible
magnitudes are refused with a reason.

**Kept on purpose**, filed for Bob; if one is fixed in `ass.py`, update `web/hrimage.py` and
the equivalence test with it:

- `abs(log(G/B) > 2)` drops G/B > e² but not G/B < e⁻² ([#21](https://github.com/PrincetonUniversity/StackImages/issues/21))
- the `shape='auto'` mask is not applied ([#22](https://github.com/PrincetonUniversity/StackImages/issues/22))
- the Color Diagram pads its upper limits from the padded lower ones ([#24](https://github.com/PrincetonUniversity/StackImages/issues/24))

Not affecting the page: `fig5axes` unset when `GUI` is False
([#23](https://github.com/PrincetonUniversity/StackImages/issues/23)); `ass.py`'s star imports need
the pinned numpy and scipy ([#25](https://github.com/PrincetonUniversity/StackImages/issues/25)).

**Timing** (Chromium, Apple silicon, 5,345 × 3,551 RGB, 217 MB): open 5 s, first diagram 32 s
(mostly `fwhm`), later diagrams 0.3 s.

## Data credits

The page footer, the download area and every CSV carry the Gaia DR3 (ESA/Gaia/DPAC, CC BY-NC 3.0
IGO, with ESA's acknowledgement and two citations), HYG v3 (David Nash, CC BY-SA 2.5; the site
ships a subset as `data/bright_stars.npy`, under the same licence) and GAVO credits. Full text: the
README's "Data and credits"; in code, `GAIA_CREDIT`, `HYG_CREDIT`, `GAVO_CREDIT` in
`web/showclusters.py` and `web/index.html`. Classroom use fits the non-commercial Gaia licence.

## Workflows

| Workflow | Runs on | Does |
|---|---|---|
| `ci.yml` | push/PR to `main` except web-only and Markdown files; manual | Desktop GUI, equivalence and packaging tests (`ci/desktop-tests.Dockerfile`, Xvfb); builds `PyWin.tar`; on `main`, releases `v<ass.py version>.<n>` if the bundle's app/data changed (`scripts/next_release.py`) |
| `pages.yml` | push/PR touching web files; manual | eslint + tsc; ruff + mypy and all web tests incl. Chromium (`ci/web-tools.Dockerfile`); verified Pyodide; site build with latest snapshots; both zips; version stamp; deploy and `web-1.N` release on `main` |
| `snapshots.yml` | manual | Fetches ~75 clusters (GAVO in pieces, ESA fallback, lower limit if both fail), validates, releases `snapshots-*`, redeploys. Resumes; re-run to fill gaps |

`install-commands.yml` runs the README's download and update commands on Windows (cmd.exe) and
macOS (zsh) runners against the live `PyWin.tar` on README changes, weekly, or on demand.
Otherwise Markdown-only changes start no workflow; the README in `PyWin.tar` refreshes with the
next desktop release. To redeploy anyway, run the workflow manually. Nothing generated is committed: Pyodide
(`build/vendor/`), the site (`_site/`) and snapshots (`build/snapshots/`) are gitignored.

## Local development

Only Docker is needed (`compose.yaml`):

```bash
docker-compose build web-tools
docker-compose run --rm web-tools python scripts/vendor_pyodide.py   # once; cached in build/vendor
docker-compose run --rm web-tools                                     # web tests (unit + Chromium)
docker-compose run --rm lint-js sh -c "npm ci && npm run check"       # eslint + tsc
docker-compose run --rm desktop-tests                                 # desktop + equivalence tests
scripts/run_web.sh                                                    # build and serve on :8000

# saved copies, to try the fallback
docker-compose run --rm web-tools python scripts/build_snapshots.py --only M13 M44 --source gavo
scripts/run_web.sh --snapshots build/snapshots
```

Without Docker: `uv venv --python 3.13 .venv && uv pip install --python .venv/bin/python -r requirements-web.lock`.

## Versions

The whole site has one number, **web 1.N**, shown on every page ("Version 1.4 (2026-11-07)")
with a link to [CHANGELOG.md](CHANGELOG.md). In `pages.yml`, `scripts/web_version.py`:

1. hashes the built site (`zip_site.py --fingerprint`, excluding build time and `version.js`);
2. keeps the newest `web-1.N` number and date if the hash matches, otherwise takes the next
   number, dated today (the major number comes from `web/VERSION`; raising it restarts at `.0`);
3. writes `version.js` into each page folder (unstamped builds say "Development build") and
   re-zips.

A PR that changes the site fails until `docs/CHANGELOG.md` has a `## 1.N (date)` entry for the
number CI names. On `main`, the new version is deployed and released as `web-1.N` with both zips
and that entry as notes; a snapshot-only redeploy without an entry gets a note naming the
snapshots. Before 1.0, releases were `web-YYYYMMDD-HHMM` and `hrdiagram-YYYYMMDD-HHMM`.

## Hosting a copy elsewhere

`StackImages-web.zip` (site root; `scripts/zip_site.py`) unzips to a `StackImages/` folder with
the home page, `ShowClusters/` (Pyodide and snapshots), `HRDiagram/` and `HOSTING.txt`.
`HRDiagram-web.zip` holds the HR page alone. Copy the folder unchanged anywhere on a static site;
links are relative, so `/~rvdb/StackImages/` works.

No server settings are needed, even where unknown extensions are served as
`application/octet-stream`:

- the worker is classic and loads `pyodide.js` with `importScripts` (no `.mjs`);
- `.wasm` is compiled from bytes if not served as `application/wasm`;
- data files (`.whl`, `.npy`, `.csv`, `.json`) are fetched as bytes;
- the HR diagram is shown through a `data:` URL, which strict CSPs allow where `blob:` is blocked.

`test_site_zip.py` serves the unzipped copy with vanderbei.princeton.edu's headers
(`tests/web/fixtures/personal_site_headers.json`: `nosniff`, its CSP, its type map,
`octet-stream` otherwise). This holds on its nginx VM and on the coming Azure Blob storage; the
new extensions only need vanderbei-deploy's `tools/gen-mime-map.sh` run so the Blob sync and
cutover parity check include them.

## Limits and safeguards

- Inputs are validated before any query (names like `m13`, `ngc5466`, `ic4665`; RA/Dec, FOV and
  magnitude in range); queries are built from numbers only.
- Uploads: ≤ 50 MB and 300,000 stars with the required Gaia columns; truncated ESA files are
  refused; a file centred elsewhere is rejected or flagged.
- Show Clusters contacts only GAVO (sync and async TAP; enforced by the browser tests) and follows
  only job IDs on GAVO's own job list (letters, digits, `-`, `_`).
- HR Diagram contacts nothing; the image never leaves the browser.
- Pyodide and every wheel are checked against pinned SHA-256 hashes at build time.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| "The live Gaia data source is unavailable…" | GAVO down or slow. Saved copies cover common clusters; otherwise use an ESA file. |
| "This field has too many stars for the live data source…" | Even the smallest pieces time out. Lower FOV or magnitude, or use an ESA file. |
| "Gaia is preparing the stars for this wide field…" | An async job; a few minutes. Keep the tab open. |
| "The live data source couldn't prepare this wide field…" | Job failed or passed 15 minutes. Retry, or use ESA (≤ 50,000 stars unless logged in). |
| "…exactly 50,000 stars…" / "…marked this file as incomplete…" | ESA cut the file off. Use Show Stars, or a smaller field or limit. |
| Stuck on "Loading Python…" | First visit downloads ~35 MB; on failure "Reload tool" shows details. |
| "The Python worker stopped" | Out of memory; smaller field or limit. |
| `pages.yml`: "No snapshots-* release yet" | Run `snapshots.yml` once. |
| `ci.yml` fails `test_web_equivalence.py` | `ass.py` changed; see "When Bob sends a new `ass.py`". |
