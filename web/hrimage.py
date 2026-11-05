"""HR diagram of your own image: a GUI-free port of the parts of ass.py it needs.

Desktop workflow being reproduced (the StackImages window):

1. "Open an RGB fits file: Select" -> get_RGB_image: read the image, its name and the
   picture shown with slidershow's sliders.
2. "Brightest star (mag)" / "Faintest star (mag)", then "HR diagram" -> GUI_HRdiag ->
   makeHRdiag_new(RGB, objname=name, minmag0, shape='auto', fignum=1, thresh=1e-5, maxmag0).
3. "Color Diagram" -> GUI_ColorDiag -> the same call with fignum=6.

Runs in Pyodide (numpy, scipy, matplotlib, astropy) and natively for tests;
tests/desktop/test_hr_equivalence.py checks it against the original functions.

The arithmetic follows ass.py statement by statement, with the same dtypes, so the
plotted points are the same. What changed is only how much is in memory at once and
which work is skipped because nothing uses it:

* findstars computes each filter, turns it into its yes/no test and frees it, instead
  of holding 18 full-size arrays. The eight neighbour filters (M1-M8) only have their
  sign tested, so they are written as the equivalent neighbour comparisons. Filters whose
  results findstars never uses (M9-M12, M21, M22) are not computed.
* findstars measures FWHM only at candidate stars, in the same order, instead of looping
  over every pixel and skipping the rest.
* makeHRdiag_new's shape='auto' cluster mask is computed in ass.py but never applied
  (``k = k&mask`` is commented out), so it is not computed here.
* Detection and photometry do not depend on the magnitude boxes, so they run once per
  image; changing the magnitudes redraws from the stored measurements.
"""

from __future__ import annotations

import io
import math
import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np

# --- Limits for the browser (see docs/WEB.md, "HR diagram of your own image") ---------
MAX_FILE_BYTES = 250 * 1024 * 1024
MAX_PIXELS = 20_000_000  # width x height
MIN_SIDE = 64  # findstars ignores 20 pixels at every edge
DISPLAY_SIZE = 1600  # longest side of the picture shown in the page
MAG_RANGE = (-30.0, 60.0)  # sanity bounds for the magnitude boxes

# ass.py globals that makeHRdiag_new and fwhm read; the desktop app never changes them
# for this path (setExptime is never called; the FWHM filter is relative, so the pixel
# scale cancels out).
EXPTIME = 1.0
ARCSEC_PER_PIXEL = 1
THRESH = 0.00001  # GUI_HRdiag's thresh


class InputError(ValueError):
    """A problem with the student's file or settings, explained in plain words."""


# --- get_RGB_image --------------------------------------------------------------------


@dataclass
class RGBImage:
    rgb: np.ndarray  # (m, n, 3) float32, exactly as get_RGB_image builds RGB
    name: str  # get_RGB_image's name: the OBJECT header, else from the file name
    ra_hours: float | None = None  # CRVAL1 / 15, when the header has it
    dec_deg: float | None = None  # CRVAL2


def _name_from_filename(filename: str) -> str:
    # ass.py: (((file.split(slash)[-1]).split('.fit')[0]).split('-')[0]).split('_')[0]
    base = re.split(r"[\\/]", filename)[-1]
    return base.split(".fit")[0].split("-")[0].split("_")[0]


def read_rgb_fits(data: bytes, filename: str) -> RGBImage:
    """get_RGB_image: the first HDU holding a 3-D (RGB) or 2-D (mono) image."""
    from astropy.io import fits

    if len(data) > MAX_FILE_BYTES:
        raise InputError(f"This file is larger than {MAX_FILE_BYTES // 2**20} MB. Use a smaller crop.")
    if not data.startswith(b"SIMPLE"):
        raise InputError("This is not a FITS file. Choose an RGB .fit or .fits image.")
    try:
        hdul = fits.open(io.BytesIO(data), memmap=False)
    except Exception as error:  # astropy raises many types for damaged files
        raise InputError("This FITS file could not be read. Is it complete?") from error
    with hdul:
        image = None
        for hdu in hdul:  # ass.py: k = 0, 1, ... until a 3-D or 2-D array turns up
            shape: tuple[int, ...] = tuple(getattr(hdu, "shape", ()))
            if hdu.is_image and len(shape) in (2, 3):
                m, n = shape[-2], shape[-1]
                if m * n > MAX_PIXELS:
                    raise InputError(
                        f"This image is {n:,} x {m:,} pixels ({m * n / 1e6:.1f} megapixels). "
                        f"The limit is {MAX_PIXELS / 1e6:.0f} megapixels; use a smaller crop."
                    )
                image = hdu.data
                break
        if image is None:
            raise InputError("This FITS file has no image in it.")
        image = np.array(image)
        if image.ndim == 3:
            if image.shape[0] < 3:
                raise InputError("This image has fewer than three color planes.")
            m, n = image.shape[1], image.shape[2]
            rgb = np.zeros((m, n, 3), dtype=np.float32)
            rgb[:, :, 0] = image[0, :, :]
            rgb[:, :, 1] = image[1, :, :]
            rgb[:, :, 2] = image[2, :, :]
        else:
            (m, n) = image.shape
            rgb = np.zeros((m, n, 3), dtype=np.float32)
            rgb[:, :, 0] = image
            rgb[:, :, 1] = image
            rgb[:, :, 2] = image
        del image
        if min(m, n) < MIN_SIDE:
            raise InputError(f"This image is too small ({n} x {m} pixels) to find stars in.")
        if not np.isfinite(rgb).all():
            raise InputError("This image has missing (NaN) or infinite pixels. Save it again without them.")
        header = hdul[0].header  # ass.py reads the primary header even for an extension image
        name, ra, dec = "", None, None
        for key in header:
            if key.upper() == "OBJECT":
                name = str(header["OBJECT"]).split("-")[0].split("_")[0]
            elif key.upper() == "CRVAL1":
                ra = float(header["CRVAL1"]) / 15
            elif key.upper() == "CRVAL2":
                dec = float(header["CRVAL2"])
    if name == "":
        name = _name_from_filename(filename)
    return RGBImage(rgb=rgb, name=name, ra_hours=ra, dec_deg=dec)


# --- findstars and fwhm -----------------------------------------------------------------


def _convolve(image: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    from scipy.signal import convolve2d

    return convolve2d(image, kernel, mode="same")


def fwhm(L: np.ndarray, i: int, j: int) -> Any:  # noqa: N803 - ass.py's names
    """ass.py's fwhm(L, i, j), statement for statement."""
    from numpy import arange, maximum, minimum, percentile

    (m, n) = L.shape
    k0 = 0
    k1 = 0
    k2 = 0
    k3 = 0
    width = np.min([i, j, m - 1 - i, n - 1 - j, 20])
    bkg = percentile(L[i - width : i + width, j - width : j + width], 50)
    v0 = L[i, j] - bkg
    for k in arange(i + 1, minimum(i + width, m - 1)):
        if L[k, j] - bkg < v0 / 2:
            v1 = L[k - 1, j] - bkg
            v2 = L[k, j] - bkg
            if v2 != v1:
                k0 = k - (i + 1) + (v0 / 2 - v1) / (v2 - v1)
            else:
                k0 = 100
            break
    for k in arange(i - 1, maximum(i - width, 0), -1):
        if L[k, j] - bkg < v0 / 2:
            v1 = L[k + 1, j] - bkg
            v2 = L[k, j] - bkg
            if v2 != v1:
                k1 = (i - 1) - k + (v0 / 2 - v1) / (v2 - v1)
            else:
                k1 = 100
            break
    for k in arange(j + 1, minimum(j + width, n)):
        if L[i, k] - bkg < v0 / 2:
            v1 = L[i, k - 1] - bkg
            v2 = L[i, k] - bkg
            if v2 != v1:
                k2 = k - (j + 1) + (v0 / 2 - v1) / (v2 - v1)
            else:
                k2 = 100
            break
    for k in arange(j - 1, maximum(j - width, 0), -1):
        if L[i, k] - bkg < v0 / 2:
            v1 = L[i, k + 1] - bkg
            v2 = L[i, k] - bkg
            if v2 != v1:
                k3 = (j - 1) - k + (v0 / 2 - v1) / (v2 - v1)
            else:
                k3 = 100
            break
    return ARCSEC_PER_PIXEL * (k0 + k1 + k2 + k3) / 2


_LAPLACIAN = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])
_LAPLACIAN2 = np.array(
    [
        [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        [-1.0, -1.0, 0.0, 0.0, 0.0, -1.0, -1.0],
        [-1.0, 0.0, 0.0, 5.0, 0.0, 0.0, -1.0],
        [-1.0, 0.0, 5.0, 8.0, 5.0, 0.0, -1.0],
        [-1.0, 0.0, 0.0, 5.0, 0.0, 0.0, -1.0],
        [-1.0, -1.0, 0.0, 0.0, 0.0, -1.0, -1.0],
        [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
    ]
)
_LAPLACIAN3 = np.array(
    [
        [0.0, -1.0, -1.0, -1.0, 0.0],
        [-1.0, 1.0, 2.0, 1.0, -1.0],
        [-1.0, 2.0, 0.0, 2.0, -1.0],
        [-1.0, 1.0, 2.0, 1.0, -1.0],
        [0.0, -1.0, -1.0, -1.0, 0.0],
    ]
)
_LAPLACIAN4 = np.array([[0.0, 3.0, 0.0], [3.0, -4.0, 3.0], [0.0, 3.0, 0.0]])  # "avoid hot pixels"


def findstars(L: np.ndarray, thresh: float) -> np.ndarray:  # noqa: N803
    """ass.py's findstars(L, thresh, fignum=0, bkg=False)[0]: the star mask k."""
    (m, n) = L.shape
    # M1-M8 (convolutions with a +1 centre and one -1 neighbour) are tested only for
    # >= 0, i.e. L[i, j] >= that neighbour. Edges are cleared below, so only interior
    # pixels matter, where the comparison and the convolution agree exactly.
    k = np.zeros((m, n), dtype=bool)
    centre = L[1:-1, 1:-1]
    inner = np.ones((m - 2, n - 2), dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                inner &= centre >= L[1 + di : m - 1 + di, 1 + dj : n - 1 + dj]
    k[1:-1, 1:-1] = inner
    del inner
    # One filter at a time: compute it, apply its test, free it.
    for kernel, relative in (
        (_LAPLACIAN, True),
        (_LAPLACIAN2, True),
        (_LAPLACIAN3, False),
        (_LAPLACIAN4, False),
    ):
        conv = _convolve(L, kernel)
        k &= conv >= thresh * np.max(conv) if relative else conv >= 0
        del conv

    k[:, 0:20] = 0
    k[:, n - 19 : n] = 0
    k[0:20, :] = 0
    k[m - 19 : m, :] = 0

    rows, cols = np.nonzero(k)  # row-major, the order ass.py's double loop visits them
    if rows.size == 0:
        raise InputError("No stars were found in this image.")
    fwhms = np.zeros(rows.size, dtype=np.float32)  # ass.py stores them in a float32 array
    for index, (i, j) in enumerate(zip(rows.tolist(), cols.tolist(), strict=True)):
        fwhms[index] = fwhm(L, i, j)
    med_fwhm = np.percentile(fwhms, 50)
    drop = (fwhms < 0.7 * med_fwhm) | (fwhms > 1.3 * med_fwhm)
    k[rows[drop], cols[drop]] = 0
    return k


# --- makeHRdiag_new ---------------------------------------------------------------------


@dataclass
class Measurements:
    """Detection and photometry: everything makeHRdiag_new computes before the magnitudes."""

    fluxsort: np.ndarray  # L flux of each detected star, faintest first (argsort order)
    Rfluxsort: np.ndarray  # noqa: N815 - ass.py's names
    Gfluxsort: np.ndarray  # noqa: N815
    Bfluxsort: np.ndarray  # noqa: N815
    xsort: np.ndarray  # pixel column of each star
    ysort: np.ndarray  # pixel row


def measure(rgb: np.ndarray, thresh: float = THRESH, exptime: float = EXPTIME) -> Measurements:
    """makeHRdiag_new up to the sorted fluxes (its result does not depend on the magnitudes)."""
    from numpy import arange, array, float32, log10, ones, sqrt  # noqa: F401 - ass.py's spelling

    R = rgb[:, :, 0] / exptime  # noqa: N806
    G = rgb[:, :, 1] / exptime  # noqa: N806
    B = rgb[:, :, 2] / exptime  # noqa: N806
    L = (R + B) / 2  # noqa: N806 - ass.py first computes (R+G+B)/3, then replaces it with this
    k = findstars(L, thresh)

    r0 = 2.5
    r1 = 2.5
    r2 = 4
    xs = arange(-r2, r2 + 1)
    ys = array([arange(-r2, r2 + 1)]).T
    radii = sqrt((ones(((2 * r2 + 1), 1)) * xs) ** 2 + (ys * ones((1, (2 * r2 + 1)))) ** 2)
    disk = float32(radii <= r0)
    annulus = float32((radii <= r2) & (radii > r1))
    annulus = np.sum(disk) * annulus / np.sum(annulus)
    averager = np.asarray(disk - annulus)  # float32, as in ass.py

    # ass.py: k = k & (R*>thresh*max) & (B*...) & (G*...) & (L*...); flux = LconvL2[k], ...
    values = {}
    for key, channel in (("R", R), ("B", B), ("G", G), ("L", L)):
        conv = _convolve(channel, averager)
        k &= conv > thresh * np.max(conv)
        values[key] = conv
    del R, G, B, L
    flux = values["L"][k]
    Rflux = values["R"][k]  # noqa: N806
    Gflux = values["G"][k]  # noqa: N806
    Bflux = values["B"][k]  # noqa: N806
    del values
    if flux.size == 0:
        raise InputError("No stars were found in this image.")
    y, x = np.nonzero(k)
    x = x.astype(np.int16)  # ass.py keeps positions as int16
    y = y.astype(np.int16)

    I = np.argsort(flux)  # noqa: E741, N806 - ass.py's name
    return Measurements(
        fluxsort=flux[I],
        Rfluxsort=Rflux[I],
        Gfluxsort=Gflux[I],
        Bfluxsort=Bflux[I],
        xsort=x[I],
        ysort=y[I],
    )


@dataclass
class Points:
    """The plotted stars, in makeHRdiag_new's order, and the values it plots them with."""

    clr: np.ndarray
    fluxmag: np.ndarray
    rgb: np.ndarray  # (N, 4) RGBA face colours
    marksize: np.ndarray
    Rflux: np.ndarray  # noqa: N815
    Gflux: np.ndarray  # noqa: N815
    Bflux: np.ndarray  # noqa: N815
    x: np.ndarray
    y: np.ndarray
    detected: int  # n in ass.py: stars measured before the magnitude and colour cuts


def plot_points(meas: Measurements, minmag0: float = 10, maxmag0: float = 20) -> Points:
    """makeHRdiag_new from the sorted fluxes to the scatter arrays, statement for statement."""
    from numpy import arange, ceil, float32, log, log10, median, percentile, sqrt, zeros

    fluxsort, Rfluxsort, Gfluxsort, Bfluxsort = (  # noqa: N806
        meas.fluxsort,
        meas.Rfluxsort,
        meas.Gfluxsort,
        meas.Bfluxsort,
    )
    Rflux = Rfluxsort.copy()  # noqa: N806 - ass.py overwrites its (equally long) Rflux array
    Gflux = Gfluxsort.copy()  # noqa: N806
    Bflux = Bfluxsort.copy()  # noqa: N806
    xs = np.zeros(len(fluxsort), dtype=meas.xsort.dtype)
    ys = np.zeros(len(fluxsort), dtype=meas.ysort.dtype)

    n = len(fluxsort)
    color = 2.5 * log10(Rfluxsort / Bfluxsort)
    color2 = zeros(n, dtype=float32)
    color2[:] = 1.0 * color[:]

    fluxmag = 2.5 * log10(fluxsort)
    maxfluxmag = np.max(fluxmag)
    fluxmag = fluxmag - maxfluxmag - minmag0
    fluxmag = -fluxmag
    maxmag = np.max(fluxmag)
    minmag = np.min(fluxmag)

    clr = zeros(n, dtype=float32)
    tmp = zeros(n, dtype=float32)
    fluxmag = zeros(n, dtype=float32)
    marksize = zeros(n, dtype=float32)
    rgb = zeros((n, 4), dtype=float32)
    for k in arange(0, n):
        flxmg = 2.5 * log10(fluxsort[k])
        flxmg = flxmg - maxfluxmag - minmag0
        flxmg = -flxmg
        dk = int(max([ceil((flxmg - 18.0) * 0), 0]))  # builtins max/min throughout, as in ass.py
        tmp[k] = median(color2[max([k - dk, 0]) : min([k + dk + 1, n + 1])])  # type: ignore[type-var]
    median_color = percentile(tmp, 70)
    msize = sqrt(640000.0 / n)
    kk = 0
    maxmag1 = max([maxmag, 20])  # type: ignore[type-var]  # builtin max, as in ass.py
    midmag = (maxmag + minmag) / 2
    for k in arange(0, n):
        flxmg = 2.5 * log10(fluxsort[k])
        flxmg = flxmg - maxfluxmag - minmag0
        flxmg = -flxmg
        if flxmg > maxmag0:
            continue
        # As written in ass.py: the second test is abs(log(G/B) > 2), the absolute value of a
        # comparison, so it drops only G/B > e^2 (not G/B < e^-2). Kept exactly.
        if (abs(log(Rfluxsort[k] / Gfluxsort[k])) > 2) | (abs(log(Gfluxsort[k] / Bfluxsort[k]) > 2)):
            continue
        bright = sqrt(0.00001 + (flxmg - maxmag) / (minmag - maxmag))
        diam = (flxmg - maxmag1) / (midmag - maxmag1)
        dk = int(max(ceil((flxmg - 18.0) * 0), 0))
        clr[kk] = median(color2[max(k - dk, 0) : min(k + dk + 1, n + 1)]) - median_color
        clr_kk = clr[kk]
        r = max([min([1 + 0.7 * clr_kk, 1]), 0])
        b = max([min([1 - 0.7 * clr_kk, 1]), 0])
        g = 1 - 0.5 * (1 - min([r, b]))
        r = bright * r
        r = min([r, 1])
        g = bright * g
        g = min([g, 1])
        b = bright * b
        b = min([b, 1])
        marksize[kk] = 0.001 + diam * msize
        fluxmag[kk] = flxmg
        rgb[kk, 0] = r
        rgb[kk, 1] = g
        rgb[kk, 2] = b
        rgb[kk, 3] = 0.8
        Rflux[kk] = Rfluxsort[k]
        Gflux[kk] = Gfluxsort[k]
        Bflux[kk] = Bfluxsort[k]
        xs[kk] = meas.xsort[k]
        ys[kk] = meas.ysort[k]
        kk += 1
    N = kk  # noqa: N806
    if N == 0:
        raise InputError(
            "No stars are between the brightest and faintest magnitudes. Raise Faintest star (mag)."
        )
    return Points(
        clr=clr[0:N],
        fluxmag=fluxmag[0:N],
        rgb=rgb[0:N, 0:4],
        marksize=marksize[0:N],
        Rflux=Rflux[0:N],
        Gflux=Gflux[0:N],
        Bflux=Bflux[0:N],
        x=xs[0:N],
        y=ys[0:N],
        detected=n,
    )


def hr_axes(points: Points, minmag0: float, maxmag0: float) -> tuple[float, float, float, float]:
    """makeHRdiag_new's fig5axes.axis([...]) (before invert_yaxis)."""
    clr = points.clr
    colorrange = np.minimum(np.maximum(max(clr), -min(clr)) + 0.4, 2.0)
    return (float(-colorrange), float(colorrange), minmag0 - 0.6, maxmag0)


def color_axes(points: Points) -> tuple[float, float, float, float]:
    """makeHRdiag_new's fig6axes.axis([...]), including its use of the already-padded xmin/ymin."""
    from numpy import log, percentile

    gr = log(points.Gflux / points.Rflux)
    bg = log(points.Bflux / points.Gflux)
    xmin = percentile(gr, 0)
    xmax = percentile(gr, 100)
    ymin = percentile(bg, 0)
    ymax = percentile(bg, 100)
    xmin = xmin - 0.1 * (xmax - xmin)
    ymin = ymin - 0.1 * (ymax - ymin)
    xmax = xmax + 0.1 * (xmax - xmin)
    ymax = ymax + 0.1 * (ymax - ymin)
    return (float(xmin), float(xmax), float(ymin), float(ymax))


# --- Figures, as makeHRdiag_new draws them -----------------------------------------------

PAGE_URL = "https://princetonuniversity.github.io/StackImages/HRDiagram/"


def _figure() -> tuple[Any, Any]:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    fig = Figure(figsize=(10, 8.8), dpi=100)  # the desktop "Images" window is 1000 x 880
    ax = fig.add_subplot(1, 1, 1, facecolor="black")
    fig.set_facecolor("black")
    for side in ("bottom", "left", "top", "right"):
        ax.spines[side].set_color("white")
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")
    ax.tick_params(axis="x", colors="white")
    ax.tick_params(axis="y", colors="white")
    return fig, ax


def _png(fig: Any, metadata: dict[str, str]) -> bytes:
    out = io.BytesIO()
    fig.savefig(out, format="png", facecolor="black", metadata=metadata)
    return out.getvalue()


def hr_figure(points: Points, name: str, minmag0: float, maxmag0: float) -> Any:
    fig, ax = _figure()
    ax.scatter(points.clr, points.fluxmag, facecolor=points.rgb, s=points.marksize, edgecolor=points.rgb / 3)
    ax.axis(hr_axes(points, minmag0, maxmag0))
    ax.set_xlabel("Color index (B-R)")
    ax.set_ylabel("Brightness (apparent magnitude)")
    ax.set_title("Hertzsprung-Russell Diagram for " + name.upper(), color="white")
    ax.invert_yaxis()
    return fig


def color_figure(points: Points, name: str) -> Any:
    from numpy import log

    fig, ax = _figure()
    ax.scatter(
        log(points.Gflux / points.Rflux),
        log(points.Bflux / points.Gflux),
        facecolor=points.rgb,
        s=points.marksize,
        edgecolor=points.rgb / 3,
    )
    ax.axis(color_axes(points))
    ax.set_xlabel("Log(G/R)")
    ax.set_ylabel("Log(B/G)")
    ax.set_title("Stellar Locus for " + name, color="white")
    return fig


# --- The picture: slidershow's starting stretch, shown smaller ----------------------------


def slidershow_input(rgb: np.ndarray) -> np.ndarray:
    """What show_stacked_image hands to slidershow: maximum(65535*RGB/np.max(RGB), 0)."""
    return np.maximum(65535 * rgb / np.max(rgb), 0)


def display_image(rgb: np.ndarray, size: int = DISPLAY_SIZE) -> dict[str, Any]:
    """show_stacked_image -> slidershow's preprocessing (showclusters.display_base) at full
    resolution, then averaged down so the longest side is at most `size` pixels for the page.
    The slider ranges and starting positions therefore match the desktop window."""
    import showclusters as sc

    base, pct_lo, pct_hi = sc.display_base(slidershow_input(rgb))
    (m, n, _) = base.shape
    factor = max(1, math.ceil(max(m, n) / size))
    h, w = m // factor, n // factor
    small = base[: h * factor, : w * factor].reshape(h, factor, w, factor, 3).mean(axis=(1, 3))
    return {
        "width": w,
        "height": h,
        "factor": factor,
        "pct_lo": pct_lo,
        "pct_hi": pct_hi,
        "image": small.astype(np.float32).tobytes(),
    }


# --- Downloads ---------------------------------------------------------------------------


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_") or "image"


def points_csv(points: Points, name: str, filename: str, minmag0: float, maxmag0: float) -> bytes:
    """One row per plotted star: where it is in the image, its fluxes, and its plotted values."""
    lines = [
        f"# HR diagram points for {name or filename}, from {filename}",
        f"# Brightest star (mag)={minmag0:g}; Faintest star (mag)={maxmag0:g}",
        f"# Made with StackImages, HR diagram of your own image: {PAGE_URL}",
        "# x and y are zero-based pixel positions (column, row) in the FITS image as stored.",
        "x,y,R_flux,G_flux,B_flux,color_index_B_minus_R,apparent_magnitude",
    ]
    for row in zip(
        points.x.tolist(),
        points.y.tolist(),
        points.Rflux.tolist(),
        points.Gflux.tolist(),
        points.Bflux.tolist(),
        points.clr.tolist(),
        points.fluxmag.tolist(),
        strict=True,
    ):
        x, y, rf, gf, bf, c, mag = row
        lines.append(f"{x},{y},{rf:.6g},{gf:.6g},{bf:.6g},{c:.4f},{mag:.4f}")
    return ("\n".join(lines) + "\n").encode()


# --- The page's session --------------------------------------------------------------------


def parse_magnitude(text: str, default: float, label: str) -> float:
    """The GUI's rule: an empty box means the default; otherwise float(value)."""
    text = str(text).strip()
    if text == "":
        return default
    try:
        value = float(text)
    except ValueError:
        raise InputError(f"{label} must be a number, such as {default:g}.") from None
    if not (math.isfinite(value) and MAG_RANGE[0] <= value <= MAG_RANGE[1]):
        raise InputError(f"{label} must be between {MAG_RANGE[0]:g} and {MAG_RANGE[1]:g}.")
    return value


@dataclass
class Session:
    image: RGBImage | None = None
    filename: str = ""
    _meas: Measurements | None = field(default=None, repr=False)

    def open_image(self, data: Any, filename: str) -> dict[str, Any]:
        """Open an RGB fits file: read it and return the picture to show."""
        self.image, self._meas = None, None
        raw = data.to_bytes() if hasattr(data, "to_bytes") else bytes(data)
        image = read_rgb_fits(raw, filename)
        del raw
        shown = display_image(image.rgb)
        self.image, self.filename = image, filename
        (m, n, _) = image.rgb.shape
        return {
            "name": image.name,
            "pixels": [n, m],
            "ra_hours": image.ra_hours,
            "dec_deg": image.dec_deg,
            **shown,
        }

    def _points(self, brightest: str, faintest: str) -> tuple[Points, float, float]:
        if self.image is None:
            raise InputError("Open an RGB FITS image first.")
        minmag0 = parse_magnitude(brightest, 10, "Brightest star (mag)")
        maxmag0 = parse_magnitude(faintest, 20, "Faintest star (mag)")
        if maxmag0 <= minmag0:
            raise InputError("Faintest star (mag) must be larger than Brightest star (mag).")
        if self._meas is None:
            self._meas = measure(self.image.rgb)
        return plot_points(self._meas, minmag0, maxmag0), minmag0, maxmag0

    def _metadata(self, kind: str, minmag0: float, maxmag0: float) -> dict[str, str]:
        assert self.image is not None
        return {
            "Title": f"{kind} for {self.image.name or self.filename}",
            "Source": self.filename,
            "Software": f"StackImages, HR diagram of your own image, {PAGE_URL}",
            "Comment": f"Brightest star (mag)={minmag0:g}; Faintest star (mag)={maxmag0:g}",
        }

    def hr_diagram(self, brightest: str, faintest: str) -> dict[str, Any]:
        points, minmag0, maxmag0 = self._points(brightest, faintest)
        assert self.image is not None
        fig = hr_figure(points, self.image.name, minmag0, maxmag0)
        png = _png(fig, self._metadata("Hertzsprung-Russell Diagram", minmag0, maxmag0))
        return {"png": png, "plotted": len(points.clr), "detected": points.detected}

    def color_diagram(self, brightest: str, faintest: str) -> dict[str, Any]:
        points, minmag0, maxmag0 = self._points(brightest, faintest)
        assert self.image is not None
        fig = color_figure(points, self.image.name)
        png = _png(fig, self._metadata("Stellar Locus", minmag0, maxmag0))
        return {"png": png, "plotted": len(points.clr), "detected": points.detected}

    def points(self, brightest: str, faintest: str) -> dict[str, Any]:
        points, minmag0, maxmag0 = self._points(brightest, faintest)
        assert self.image is not None
        data = points_csv(points, self.image.name, self.filename, minmag0, maxmag0)
        return {"name": f"{_slug(self.image.name or self.filename)}_HR_points.csv", "data": data}
