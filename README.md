StackImages
============

StackImages is a Python script and graphical interface for astrophotography image integration.

This README includes installation and usage instructions adapted from the project's documentation:
https://orfe.io/stackimg.html

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
A browser version of ShowClusters (no installation) is published at
https://princetonuniversity.github.io/StackImages/. It is built from `web/` by
`.github/workflows/pages.yml`; `ass.py` is unchanged. To build and test it locally:

```bash
uv venv --python 3.13 .venv && uv pip install --python .venv/bin/python -r requirements-web.lock
.venv/bin/python scripts/vendor_pyodide.py          # pinned, verified Pyodide into build/vendor
.venv/bin/python -m playwright install chromium
.venv/bin/python -m pytest tests/web               # unit + browser tests
.venv/bin/python scripts/build_site.py --no-snapshots && .venv/bin/python -m http.server -d _site 8000
```

Notes
-----
- The distributed release artifact produced by CI is named `PyWin.tar` to match legacy expectations.
  It is built by `scripts/package_pywin.sh` and contains only the desktop app and its data under `PyWin/`.
- Maintainers: run the desktop tests locally with `docker-compose run --rm desktop-tests`.
  `environment.yml` and `ci/conda-linux-64.lock` define the CI environment; regenerate the lock with `scripts/lock_conda_env.sh`.
- This README text is adapted from https://orfe.io/stackimg.html. For the latest instructions or additional details, see that page.
