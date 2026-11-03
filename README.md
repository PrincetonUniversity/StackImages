StackImages
============

StackImages is a Python script and graphical interface for astrophotography image integration.

**Show Clusters in the browser:** https://princetonuniversity.github.io/StackImages/ShowClusters/ — the
ShowClusters tool (Gaia star fields and HR diagrams) with nothing to install. See
[Show Clusters in the browser](#show-clusters-in-the-browser) below.

macOS installation (copy-paste into Terminal)
--------------------------------------------
1. Download and extract StackImages:

```bash
curl -L https://vanderbei.princeton.edu/tmp/PyWin.tar -o $HOME/Downloads/StackImages.tar && \
if [[ ! -d $HOME/Downloads/StackImages ]]; then mkdir -p $HOME/Downloads/StackImages ; fi && \
tar -xvf $HOME/Downloads/StackImages.tar --directory=$HOME/Downloads/StackImages --strip-components=1 PyWin
```

2. Download Miniconda:

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-`arch`.sh -o ~/Downloads/miniconda.sh
```

3. Install Miniconda:

```bash
chmod +x ~/Downloads/miniconda.sh ; ~/Downloads/miniconda.sh -bu
```

4. Initialize conda (then restart Terminal):

```bash
~/miniconda3/bin/conda init bash ; ~/miniconda3/bin/conda init zsh
```

5. Accept the conda Terms of Service:

```bash
conda tos accept
```

6. Configure conda channels (prioritize astropy/apple/conda-forge):

```bash
conda config --add channels conda-forge --add channels apple --add channels astropy
```

7. Create the conda environment `BobsAstro`:

```bash
conda create -n BobsAstro astropy astroquery scipy==1.11.2 matplotlib wxpython opencv pyqt qt=5 imageio -y
```

Usage (macOS)
--------------
Activate the environment and run the helper to download FITS files:

```bash
conda activate BobsAstro
pythonw $HOME/Downloads/StackImages/GetFit.py
```

Run the main GUI:

```bash
conda activate BobsAstro
pythonw $HOME/Downloads/StackImages/StackImages.py
```

Update the software
-------------------
To update StackImages to the latest versions posted on the server, run:

```bash
curl -L https://vanderbei.princeton.edu/tmp/PyWin.tar -o $HOME/Downloads/StackImages.tar && \
tar -xvkf $HOME/Downloads/StackImages.tar --unlink -C $HOME/Downloads/StackImages --strip-components=1 --exclude='._*'
```

Windows installation (copy-paste into Anaconda Prompt)
------------------------------------------------------
1. Download and extract StackImages:

```powershell
curl --noproxy "*" https://vanderbei.princeton.edu/tmp/PyWin.tar -o %USERPROFILE%\Downloads\StackImages.tar ; \
IF NOT EXIST %USERPROFILE%\Downloads\StackImages mkdir %USERPROFILE%\Downloads\StackImages ; \
 tar -xvf %USERPROFILE%\Downloads\StackImages.tar --directory=%USERPROFILE%\Downloads\StackImages --strip-components=1 PyWin
```

2. Configure conda channels:

```powershell
conda config --add channels conda-forge --add channels astropy
```

3. Accept the Terms of Service:

```powershell
conda tos accept
```

4. Create the `BobsAstro` environment:

```powershell
conda create -n BobsAstro astropy astroquery scipy==1.11.2 matplotlib wxpython pyqt qt=5 imageio -y
```

5. Activate the environment and install OpenCV via pip (Windows workaround):

```powershell
conda activate BobsAstro
pip install opencv-contrib-python
```

Usage (Windows)
---------------
Run GetFit and StackImages from the activated environment:

```powershell
conda activate BobsAstro
pythonw %USERPROFILE%\Downloads\StackImages\GetFit.py
pythonw %USERPROFILE%\Downloads\StackImages\StackImages.py
```

Show Clusters in the browser
----------------------------
https://princetonuniversity.github.io/StackImages/ShowClusters/ runs ShowClusters in any modern browser
(Chrome, Edge, Firefox, Safari) on macOS, Windows, Linux or a Chromebook:

1. Type an object (`m13`, `NGC5466`, `IC4665`) and click **Look up** (or press Enter); RA and Dec fill in.
2. Set the field of view and maximum magnitude, then click **Show Stars**.
3. Adjust the picture with the sliders (**Reset sliders** puts them back; **Full screen** shows the
   picture as large as the screen allows, Esc or ✕ to return); choose Apparent or Absolute and
   click **HR Diagram**.
4. Use **Download HR diagram (PNG)** to save it. **Download star data (CSV)** saves the exact Gaia
   stars used (upload it again to reproduce the result); **Download HR points (CSV)** saves the
   diagram's points.

How it works: the page runs Python in the browser with Pyodide. `web/showclusters.py` is a
GUI-free port of the ShowClusters functions in `ass.py` (which is unchanged), checked against
the original by an automated equivalence test. Star data comes live from GAVO's copy of Gaia
DR3, with saved copies and a manual ESA-archive upload as fallbacks.

Details for maintainers (architecture, data sources, whether star counts match the desktop
app, workflows, local development, troubleshooting): [docs/WEB.md](docs/WEB.md).

The site root, https://princetonuniversity.github.io/StackImages/, is a short StackImages home page
linking to Show Clusters and to the desktop app.

Quick local run (Docker only):

```bash
scripts/run_web.sh            # downloads Pyodide once, builds _site/, serves http://127.0.0.1:8000/
docker-compose run --rm web-tools                         # web unit + browser tests
docker-compose run --rm lint-js sh -c "npm ci && npm run check"   # eslint + tsc
```

Data and credits
----------------
Show Clusters (desktop and web) uses star data from these sources; please keep the credits when
sharing results or the downloaded CSV files (which carry them in their `#` header lines).

- **Gaia DR3** — ESA/Gaia/DPAC, licensed
  [CC BY-NC 3.0 IGO](https://www.cosmos.esa.int/web/gaia-users/license). ESA asks for this
  acknowledgement: *"This work has made use of data from the European Space Agency (ESA) mission
  Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis
  Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has
  been provided by national institutions, in particular the institutions participating in the Gaia
  Multilateral Agreement."* Cite Gaia Collaboration, Prusti et al. 2016, A&A 595, A1 and Gaia
  Collaboration, Vallenari et al. 2023, A&A 674, A1
  ([credit and citation instructions](https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_credit_and_citation_instructions/)).
- **HYG database v3** (`csv/hygdata_v3.csv`, the bright stars) — David Nash,
  [codeberg.org/astronexus/hyg](https://codeberg.org/astronexus/hyg), licensed
  [CC BY-SA 2.5](https://creativecommons.org/licenses/by-sa/2.5/).
- **GAVO** — the German Astrophysical Virtual Observatory serves the web page's live Gaia data
  ([dc.g-vo.org](https://dc.g-vo.org/), table `gaia.dr3lite`).

Notes
-----
- The distributed release artifact produced by CI is named `PyWin.tar` to match legacy expectations.
  It is built by `scripts/package_pywin.sh` and contains only the desktop app and its data under `PyWin/`.
- Releases are versioned `v<ass.py version>.<n>` (e.g. `v21.0`, `v21.1`; a new `version = N` in `ass.py` starts
  `vN.0`). CI publishes one only when the app or data files in `PyWin.tar` changed since the last release;
  README/LICENSE edits alone don't trigger one. Run CI (desktop) manually on `main` to release on demand.
- Maintainers: run the desktop tests locally with `docker-compose run --rm desktop-tests`.
  `environment.yml` and `ci/conda-linux-64.lock` define the CI environment; regenerate the lock with `scripts/lock_conda_env.sh`.
- CI: `ci.yml` (desktop app, `PyWin.tar`), `pages.yml` (web page) and `snapshots.yml` (web fallback data) are
  independent workflows; see [docs/WEB.md](docs/WEB.md#workflows).
