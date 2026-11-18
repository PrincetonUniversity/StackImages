StackImages
===========

StackImages is a Python script and graphical interface for astrophotography image integration.

Two of its tools also run in the browser, with nothing to install: [Show Clusters](#show-clusters-in-the-browser)
and [HR diagram of your own image](#hr-diagram-of-your-own-image-in-the-browser).

macOS
-----
Paste each command into Terminal (Applications → Utilities).

1. Download and extract StackImages:

   ```bash
   curl -fL https://github.com/PrincetonUniversity/StackImages/releases/latest/download/PyWin.tar -o $HOME/Downloads/StackImages.tar && mkdir -p $HOME/Downloads/StackImages && tar -xvf $HOME/Downloads/StackImages.tar --directory=$HOME/Downloads/StackImages --strip-components=1 PyWin
   ```

2. Download and install Miniconda:

   ```bash
   curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-`arch`.sh -o ~/Downloads/miniconda.sh
   chmod +x ~/Downloads/miniconda.sh ; ~/Downloads/miniconda.sh -bu
   ```

3. Initialize conda, then quit and reopen Terminal:

   ```bash
   ~/miniconda3/bin/conda init bash ; ~/miniconda3/bin/conda init zsh
   ```

4. Accept the Terms of Service, set the channels (astropy, then apple, then conda-forge) and create
   the `BobsAstro` environment:

   ```bash
   conda tos accept
   conda config --add channels conda-forge --add channels apple --add channels astropy
   conda create -n BobsAstro astropy astroquery scipy==1.11.2 matplotlib wxpython opencv pyqt qt=5 imageio -y
   ```

Run GetFit (downloads FITS images) or StackImages:

```bash
conda activate BobsAstro
pythonw $HOME/Downloads/StackImages/GetFit.py
pythonw $HOME/Downloads/StackImages/StackImages.py
```

Update to the latest release:

```bash
curl -fL https://github.com/PrincetonUniversity/StackImages/releases/latest/download/PyWin.tar -o $HOME/Downloads/StackImages.tar && tar -xvkf $HOME/Downloads/StackImages.tar --unlink -C $HOME/Downloads/StackImages --strip-components=1 --exclude='._*'
```

Windows
-------
1. [Download and install Miniconda](https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe),
   then open **Anaconda Prompt (Miniconda3)** from the Start menu. Paste each command below into it.
2. Download and extract StackImages:

   ```bat
   curl -fL --noproxy "*" https://github.com/PrincetonUniversity/StackImages/releases/latest/download/PyWin.tar -o %USERPROFILE%\Downloads\StackImages.tar && IF NOT EXIST %USERPROFILE%\Downloads\StackImages mkdir %USERPROFILE%\Downloads\StackImages && tar -xvf %USERPROFILE%\Downloads\StackImages.tar --directory=%USERPROFILE%\Downloads\StackImages --strip-components=1 PyWin
   ```

3. Accept the Terms of Service, set the channels (astropy, then conda-forge) and create the
   `BobsAstro` environment. OpenCV comes from pip, because conda's `opencv` fails to import on Windows.

   ```bat
   conda tos accept
   conda config --add channels conda-forge --add channels astropy
   conda create -n BobsAstro astropy astroquery scipy==1.11.2 matplotlib wxpython pyqt qt=5 imageio -y
   conda activate BobsAstro
   pip install opencv-contrib-python
   ```

Run GetFit or StackImages:

```bat
conda activate BobsAstro
pythonw %USERPROFILE%\Downloads\StackImages\GetFit.py
pythonw %USERPROFILE%\Downloads\StackImages\StackImages.py
```

Update to the latest release:

```bat
curl -fL --noproxy "*" https://github.com/PrincetonUniversity/StackImages/releases/latest/download/PyWin.tar -o %USERPROFILE%\Downloads\StackImages.tar && tar -xvkf %USERPROFILE%\Downloads\StackImages.tar --unlink -C %USERPROFILE%\Downloads\StackImages --strip-components=1 --exclude=._*
```

Start over
----------
To remove the environment, run `conda env remove -n BobsAstro`, then recreate it with the
`conda create` step above.

Show Clusters in the browser
----------------------------
https://princetonuniversity.github.io/StackImages/ShowClusters/ (Chrome, Edge, Firefox or Safari):

1. Type an object (`m13`, `NGC5466`, `IC4665`) and press Enter or click **Look up**.
2. Set the field of view and maximum magnitude, then click **Show Stars**.
3. Adjust the picture with the sliders (⤢ shows it full screen). Choose Apparent or Absolute and
   click **HR Diagram**.
4. Download the diagram (PNG), its points (CSV) or the star data (CSV). Uploading the star-data
   file again reproduces the same picture and diagram.

The page runs a port of `ass.py`'s ShowClusters (`web/showclusters.py`) in the browser with
Pyodide, using live Gaia DR3 data from GAVO.

HR diagram of your own image in the browser
-------------------------------------------
https://princetonuniversity.github.io/StackImages/HRDiagram/ does what the desktop **Open an RGB
fits file**, **HR diagram** and **Color Diagram** buttons do:

1. Choose an RGB (or mono) FITS image, up to 250 MB and 20 megapixels. It never leaves your computer.
2. Adjust the picture with the sliders.
3. Set **Brightest star (mag)** (default 10) and **Faintest star (mag)** (default 20), then click
   **HR diagram** or **Color Diagram**. The first diagram of a large image can take a minute.
4. Download the diagram (PNG) or its points (CSV).

Web pages for maintainers
-------------------------
- How they work, data sources, tests and workflows: [docs/WEB.md](docs/WEB.md).
- Versions (web 1.0, 1.1, …) are shown on every page, described in [docs/CHANGELOG.md](docs/CHANGELOG.md)
  and released as `web-1.N` with the site zips.
- To host a copy elsewhere, unzip
  [StackImages-web.zip](https://princetonuniversity.github.io/StackImages/StackImages-web.zip) (or
  [HRDiagram-web.zip](https://princetonuniversity.github.io/StackImages/HRDiagram-web.zip) for the
  HR page alone) and upload the folder unchanged.
- Local run and tests (Docker only):

  ```bash
  scripts/run_web.sh                                              # build and serve http://127.0.0.1:8000/
  docker-compose run --rm web-tools                               # web unit + browser tests
  docker-compose run --rm lint-js sh -c "npm ci && npm run check" # eslint + tsc
  docker-compose run --rm desktop-tests                           # desktop + equivalence tests
  ```

Desktop releases
----------------
- CI (`ci.yml`) builds `PyWin.tar` (`scripts/package_pywin.sh`) with only the desktop app and data.
- Releases are `v<ass.py version>.<n>` (`v21.0`, `v21.1`, …), published only when the app or data
  in `PyWin.tar` changed. Run `ci.yml` manually on `main` to release on demand.
- `environment.yml` and `ci/conda-linux-64.lock` define the CI environment; regenerate the lock
  with `scripts/lock_conda_env.sh`.
- `install-commands.yml` runs this README's download and update commands on Windows and macOS
  runners (`scripts/readme_install.py`) when the README changes and weekly. Other Markdown-only
  changes start no workflow.

Data and credits
----------------
Keep these credits when sharing results; the downloaded CSV files carry them in their `#` lines.

- **Gaia DR3**: ESA/Gaia/DPAC, [CC BY-NC 3.0 IGO](https://www.cosmos.esa.int/web/gaia-users/license).
  ESA asks for this acknowledgement: *"This work has made use of data from the European Space
  Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing
  and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for
  the DPAC has been provided by national institutions, in particular the institutions participating
  in the Gaia Multilateral Agreement."* Cite Gaia Collaboration, Prusti et al. 2016, A&A 595, A1 and
  Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1
  ([instructions](https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_credit_and_citation_instructions/)).
- **HYG database v3** (`csv/hygdata_v3.csv`, the bright stars): David Nash,
  [codeberg.org/astronexus/hyg](https://codeberg.org/astronexus/hyg),
  [CC BY-SA 2.5](https://creativecommons.org/licenses/by-sa/2.5/).
- **GAVO** ([dc.g-vo.org](https://dc.g-vo.org/), table `gaia.dr3lite`) serves the web page's live Gaia data.
