# ShowClusters Web: Plan and Status

2026-11-02 · Michael Bino · status updated 2026-11-05

Goal: ShowClusters as a static Pyodide page on GitHub Pages, so students do HW #7 from a link
instead of installing conda and wxPython. Due for Bob's review Wed Nov 4; homework weekend Nov 7–8.

Live at https://princetonuniversity.github.io/StackImages/ShowClusters/. What was built, and where
it differs from this plan: [WEB.md](WEB.md). The full original plan is in git history.

## Status

| Item | State |
|---|---|
| Desktop CI fixed; `PyWin.tar` packaging | Done (PR #1) |
| `web/showclusters.py` port, equivalence test vs `ass.py`, page, GAVO live data, ESA upload, Pages deploy | Done (PR #2) |
| Snapshot fallback (`snapshots.yml`); eslint + tsc; web tests in a container | Done (PR #3) |
| Bob reviews the page with HW #7 and updates the homework text | Open |
| Bob's decisions below | Open |

Added since: the site zip for other hosts, HR Diagram of Your Own Image, and web version numbers.

## Decisions made

- The page lives in this repo under `web/`; the desktop app stays untouched so Bob can keep
  sending whole new `ass.py` files.
- `web/showclusters.py` ports only what ShowClusters needs; an equivalence test against the real
  `ass.py` catches any math change.
- JavaScript makes all network requests; Python stays pure and natively testable.
- Pyodide, snapshots and the built site are never committed; Pyodide is SHA-256 verified.
- Live data from GAVO (ESA blocks cross-site reads), with saved snapshots and a manual ESA upload
  as fallbacks; every result is labelled with its source.
- `ass.py`'s quirks are kept for parity with last year (see WEB.md).

## Open questions for Bob

- Which clusters the snapshots should cover.
- Maximum FOV and magnitude limits.
- Whether to keep the `ass.py` quirks.
- Homework text: the link instead of install steps; download buttons instead of "Save Figure".
- `png/Coathanger.png` (40 MB) ships in `PyWin.tar` but no code uses it: keep or drop?
