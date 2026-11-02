# ShowClusters Web — Implementation Plan

2026-11-02 · Michael Bino

Ship ShowClusters as a static Pyodide page on GitHub Pages so students can do HW #7 from a
link, with no conda, wxPython or StackImages install. HW #7 is due the weekend of Nov 7–8,
so the page must be live and checked by Bob by **Wed Nov 4**.

## Status (updated 2026-11-02)

Live at https://princetonuniversity.github.io/StackImages/. How it works: [WEB.md](WEB.md).

| Item | State |
|---|---|
| Desktop CI fixed; `PyWin.tar` packaging script | Done (PR #1) |
| `web/showclusters.py` port + equivalence test vs `ass.py` | Done (PR #2) |
| Page, worker, GAVO live data, ESA upload | Done (PR #2) |
| Pages workflow and deployment | Done (PR #2) |
| Snapshot builder with ESA → GAVO fallback; `snapshots.yml` | Done (PR #3); first run pending |
| eslint + tsc in CI; web tests in a container | Done (PR #3) |
| Bob reviews the page with HW #7; homework text updated | Open |
| Bob's decisions in §9 (snapshot list, caps, quirks, `Coathanger.png`) | Open |

The rest of this document is the original plan; where it differs, WEB.md describes what was built.

## 1. Scope

In scope:

- Object lookup by Messier / NGC / IC name → RA, Dec (`getRaDec`)
- Gaia star query by RA, Dec, FOV, max magnitude (`get_star_info`)
- Simulated star-field image with display sliders (`show_gaia_image`, 2D only)
- HR diagram, apparent or absolute magnitude, PNG download (`GUI_makeHRdiag0`)
- Three data sources: GAVO live, saved-file fallback, manual ESA upload

Out of scope:

- The rest of `ass.py` (stacking, alignment, astrometry, deconvolution, etc.)
- 3D and "future" image types, and therefore the Animation Rate (proper motion) field
- Color Diagram (not used by the homework; feasible later from ESA uploads, which carry fluxes)
- Changes to the desktop app, except the separate CI fix (§10)

## 2. Architecture

**Decided:** the page lives in the StackImages repo under `web/`, published to
`https://princetonuniversity.github.io/StackImages/`. The desktop app (`ass.py`,
`ShowClusters.py`, `StackImages.py`, `GetFit.py`, `csv/`, `fit/`, …) stays **untouched**, so
Bob can keep sending whole new versions of `ass.py`.

The web version is a port (option 1): `web/showclusters.py` re-implements only the functions
ShowClusters needs, without wx/Qt/cv2/astroquery. An equivalence test (§7) runs the original
`ass.py` functions and the port on the same inputs and fails on any difference, so drift is
caught when Bob changes the math.

```
StackImages/
├── ass.py, ShowClusters.py, StackImages.py, GetFit.py, csv/, fit/, png/, txt/   unchanged
├── web/                 index.html, app.js, worker.js, style.css, showclusters.py
├── scripts/
│   ├── package_pywin.sh     minimal desktop bundle (done, PR #1)
│   ├── lock_conda_env.sh    desktop CI lockfile (done, PR #1)
│   ├── vendor_pyodide.py    download pinned Pyodide + wheels, verify sha256 (not committed)
│   └── build_site.py        assemble _site/: web/ + vendor + generated data
├── tests/desktop/       GUI smoke tests (wx under Xvfb)
├── tests/packaging/     PyWin.tar tests
├── tests/web/           port unit tests, equivalence vs ass.py, JS + browser tests
└── .github/workflows/   ci.yml (desktop, existing)  ·  pages.yml (web, new, independent)
```

Runtime in the browser, same pattern as `pu-shd/hrdiagram`:

```
index.html / app.js (UI, sliders, file drop, downloads, all network fetches)
      │ postMessage
worker.js ── Pyodide 0.29.3 (numpy, matplotlib, astropy)
      │ runs
showclusters.py   lookup · query building · parsing · filtering · render · HR plot
```

- **Pyodide is never committed** (hrdiagram doesn't commit it either: `web/vendor/` is
  gitignored and fetched by a script). `vendor_pyodide.py` pins the release and verifies
  sha256 for every file, including the core runtime files, which hrdiagram's script does not check.
- **Generated data is never committed.** `build_site.py` derives it from files already in the
  repo during the Pages build:
  - `csv/hygdata_v3.csv` (32 MB) → only rows `ass.py` can use (−5 < mag < 5.5) → small JSON
  - `fit/sxv_psf4-RGB.fit` (3 MB) → the 500×500 crop used by `show_gaia_image`, downsampled
  - Messier / NGC / IC CSVs copied as is (~400 KB)
  - Snapshots (§3) queried from ESA with astroquery during the build
- Network fetches happen in JS and bytes are handed to Python, so Python stays pure and
  testable natively. astroquery is not used in the browser (it can't run in Pyodide).
- `PyWin.tar` excludes `web/`, `scripts/`, `tests/`, `ci/`, `docs/` (enforced by
  `package_pywin.sh` and its tests).

## 3. Data sources

The browser can only read responses from servers that send CORS headers. Tested 2026-11-02:

| Server | CORS | Result |
|---|---|---|
| ESA Gaia archive `gea.esac.esa.int` (what astroquery uses) | No — preflight `403` | Not usable from the page |
| GAVO `dc.g-vo.org`, table `gaia.dr3lite` | Yes (echoes Origin) | M13 0.2°/mag 18: 490 KB, 2.5 s · M44 1.5°/mag 18: 480 KB, 4.3 s |

Order on "Show Stars":

1. **GAVO live** (default). POST to `https://dc.g-vo.org/tap/sync`, `FORMAT=csv`, explicit
   `MAXREC`; timeout ~30 s. Columns: `ra, dec, parallax, phot_g_mean_mag, phot_bp_mean_mag,
   phot_rp_mean_mag`; `bp_rp = bp − rp` (Gaia's own definition of `bp_rp`).
   If the row count equals `MAXREC`, warn that results were truncated.
2. **Saved snapshot** (automatic fallback) when GAVO errors or times out *and* a snapshot
   exists whose stored FOV and magnitude limit cover the request; filter it client-side.
   Snapshots are built with astroquery from **ESA itself** during the Pages build (no CORS
   issue outside a browser), so the fallback matches last year exactly. They are deployed, not
   committed; if ESA is unreachable at build time, the build fails rather than shipping without them.
3. **Manual ESA upload** (always available, and offered when 1 and 2 fail):
   - Page shows the ADQL for the current inputs (same query as `get_star_info`) with a
     "Copy query" button and a link to the ESA Gaia archive's ADQL search page.
   - Student runs it, downloads VOTable, CSV or FITS, drops the file on the page.
   - Parsed in Pyodide with astropy (VOTable, FITS) or csv.

Every result carries a **source label** shown in the UI and in the HR diagram title:
"Gaia DR3 via GAVO", "Gaia DR3 saved copy (YYYY-MM-DD)", "Gaia DR3 uploaded ESA file".

Snapshot set (to confirm with Bob): M13, M44 and the other Messier globular and open clusters,
plus a short list of popular NGC/IC clusters. Size budget ≤ 30 MB total, stored as compressed
float32 arrays.

## 4. Porting from `ass.py`

Port function by function into `showclusters.py`, keeping the numerical rules unchanged:

| `ass.py` | Port | Notes |
|---|---|---|
| `getRaDec`, `decimal_to_min_sec`, `min_sec_to_decimal` | As is | Input validated first (§6) |
| `get_star_info` query | Same ADQL; radius `fov/√2`; `G < maxmag` | Built from validated floats only |
| `get_star_info` filters | Keep: skip masked/NaN, keep `mag ≥ 4.5 and bp_rp < 3`; `absmag = mag − 2.5·log10(100²/plx²)` with `abs(parallax)` | |
| Hipparcos merge | Keep: stars within `fov·√2`, `−5 < mag < 5.5` | Uses HYG B−V colour mixed with Gaia BP−RP — existing behaviour, kept |
| `show_gaia_image` render | 2D only; 4000 px → ~1200 px with PSF scaled by the same factor | 4000² float32 RGB with padding is ~240 MB and a per-star Python loop — too heavy in wasm |
| `slidershow` controls | Replicated in JS on the rendered 8-bit image | Confirm which controls it exposes during the port |
| `GUI_makeHRdiag0` | matplotlib Agg → PNG; same colours, sizes, axes, title | Download button replaces "Save Figure" |

Behaviours to decide with Bob (default: keep, for parity with last year):

- "Show Stars" resets magnitude type to Apparent when FOV ≤ 5°; students must pick Absolute *after* it.
- In Absolute mode the brightest star is replaced by the last one before plotting.
- The HR loop starts at index 1, so the first star is never drawn.
- "Number of plotted stars" counts stars inside the square frame, not the query circle.

## 5. UI and homework flow

| HW #7 step (last year) | Web page |
|---|---|
| `conda activate BobsAstro`, `pythonw -m ShowClusters` | Open the link; status shows "Loading Python…" then "Ready" |
| Type `m13`, press Enter → RA/Dec fill in | Same; Enter or Tab triggers lookup |
| FOV (default 0.3), Max Magnitude (default 16) | Same fields and defaults |
| "Show Stars" → picture, sliders; terminal prints star count | Same; star count, brightest star and source label shown on the page |
| "HR Diagram" → Figure 6, "Save Figure" | HR diagram shown inline, "Download PNG" |
| Apparent / Absolute choice | Same dropdown |

Bob should update the homework text: link instead of install steps, "Download PNG" instead of
"Save Figure", and a note that counts may differ slightly from last year when the GAVO source is used.

## 6. Validation and security

- Object name: `^(m|ngc|ic)\s*\d{1,4}$` (case-insensitive); unknown → clear message, no query.
- RA 0–24 h, Dec −90–90°, FOV 0.01–5°, max magnitude 4.5–21; out of range → rejected with the allowed range.
- ADQL is built only from parsed floats, never from raw text.
- Uploads: ≤ 50 MB, ≤ 300k rows, required columns present (`ra, dec, parallax,
  phot_g_mean_mag, bp_rp`), numeric and finite after the same masking rules; anything else
  rejected with a message naming what's missing.
- Upload coverage check: warn if the file's median position is more than FOV from the requested centre.
- No cookies, no analytics, nothing about the student is sent anywhere except the query to GAVO.

## 7. Testing and CI

Every test asserts a specific outcome; empty output, missing figures or zero rows are failures.

- **Python (pytest, native and in Docker):**
  - lookup: M13, M44, NGC5466, IC cases, unknown names, bad input
  - parsing: recorded GAVO CSV for M13; ESA VOTable, CSV and FITS fixtures; malformed and truncated files
  - filtering and absmag: golden values on a small hand-built table
  - equivalence vs `ass.py` (option 1's safety net): in the desktop container (which has wx,
    cv2, astroquery), import the original `ass.py`, replace `astroquery.gaia.Gaia` with a fake
    that returns a recorded ESA result, run `get_star_info`, `getRaDec`, `gnomonic` and the
    HR-diagram arrays, and compare with the port on the same input. Any difference fails, so a
    math change in a new `ass.py` from Bob is caught.
  - render and HR plot: produce a non-empty PNG of the expected size; plotted count matches
- **JS (same runner as hrdiagram):** source selection with mocked `fetch` — GAVO OK; GAVO fails →
  snapshot; snapshot misses → upload prompt; `MAXREC` truncation warning.
- **End-to-end (Playwright in Docker, network mocked):** the HW #7 M13 flow from page load to PNG download.
- **One live smoke test** against GAVO, run manually and on a schedule, not on every push.
- **CI, two independent workflows** (path filters, neither blocks the other):
  - `ci.yml` (desktop, existing, fixed in PR #1): desktop + packaging + equivalence tests in the
    pinned conda container under Xvfb; builds `PyWin.tar`; release on push to `main`.
  - `pages.yml` (web, new): ruff + mypy on `web/showclusters.py`, `scripts/`, `tests/web`;
    eslint + tsc on the JS; pytest and browser tests in a pinned web-tools container; then
    vendor, build `_site/` and deploy on `main`.
  - Dependencies pinned (Python lockfile, `package-lock.json`, conda lock); base images pinned.
  - Lint applies to new code only; Bob's `ass.py` is not linted.

## 8. Deployment

- GitHub Pages via `pages.yml` (`actions/upload-pages-artifact` + `actions/deploy-pages`) on
  push to `main` and manual dispatch; Pages source = "GitHub Actions".
- **Prerequisite:** a repo admin enables Pages for `PrincetonUniversity/StackImages`
  (currently not enabled; the org may restrict Pages).
- `_site/` = `web/` + `vendor/` + generated `data/` + a manifest with sha256 and build date.
- zsh `scripts/run_web.sh` serves `_site/` locally for testing; `--dry-run`, `set -euo pipefail`, re-runnable.

## 9. Risks and open questions

| Risk / question | Mitigation / owner |
|---|---|
| GAVO down or rate-limited during homework week | Snapshot fallback + ESA upload |
| GAVO drops CORS in future | Same fallbacks; only a proxy would restore live ESA access |
| `gaia.dr3lite` differs from `gaiadr3.gaia_source` for our columns | Compare M13/M44 counts against ESA before launch |
| Rendering speed on low-end laptops / Chromebooks | 1200 px render; test on a slow machine |
| Snapshot cluster list | Bob |
| Max FOV and magnitude caps | Bob |
| Keep the four `ass.py` quirks in §4? | Bob |
| Pages enablement in the PrincetonUniversity org | Michael (repo admin / org owner) |
| `png/Coathanger.png` (40 MB) is in `PyWin.tar` but no code references it | Bob: keep or drop |

## 10. Milestones

1. **Done (PR #1):** desktop CI fixed; `PyWin.tar` packaging script with tests.
2. **Nov 2–3:** port `web/showclusters.py`; equivalence test vs `ass.py`; GAVO path; render and
   HR diagram working in the browser; Python tests.
3. **Nov 3:** `vendor_pyodide.py`, `build_site.py` (bright stars, PSF, snapshots); fallback and
   upload paths; JS and browser tests; `pages.yml` green.
4. **Nov 4:** Pages enabled and deployed; Bob runs HW #7 Problems 1–3 on it; fix findings; Bob
   updates the homework text.
5. **Nov 5:** Homework released with the link.
