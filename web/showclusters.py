"""ShowClusters for the browser: a GUI-free port of the parts of ass.py it needs.

Runs in Pyodide (numpy, matplotlib, astropy) and natively for tests. The numerical
rules follow ass.py (getRaDec, get_star_info, show_gaia_image, slidershow,
GUI_makeHRdiag0); tests/desktop/test_web_equivalence.py checks them against the
original. Network access is done by the page's JavaScript, never here.
"""

from __future__ import annotations

import csv
import io
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

# --- Limits for user input (validated before anything else runs) -------------------

FOV_RANGE = (0.01, 5.0)  # degrees; ass.py switches to a different mode above 5
MAXMAG_RANGE = (4.5, 21.0)  # Gaia G; ass.py drops stars brighter than 4.5
MAX_UPLOAD_BYTES = 50 * 1024 * 1024
MAX_ROWS = 300_000
RENDER_SIZE = 1200  # ass.py renders 4000 px; scaled down for browsers
ORIGINAL_SIZE = 4000
ORIGINAL_PSF_HALF_WIDTH = 250

_NAME = re.compile(r"^\s*(m|ngc|ic)\s*(\d{1,4})\s*$", re.IGNORECASE)
_CATALOGS = {
    "m": ("NI2021_MessierObjects.csv", 1),
    "ngc": ("NI2021_NGCObjects.csv", 3),
    "ic": ("ICObj.csv", 2),
}
# Columns ass.py reads from the ESA query, in its order. The first seven must be present.
REQUIRED = ("ra", "dec", "parallax", "phot_g_mean_mag", "bp_rp", "pmra", "pmdec")


class InputError(ValueError):
    """A problem with what the student entered or uploaded; message is shown as is."""


# --- Coordinates ---------------------------------------------------------------------


@dataclass(frozen=True)
class ObjectPosition:
    name: str
    ra_hours: float
    dec_deg: float


def lookup_object(name: str, catalog_dir: Path) -> ObjectPosition:
    """ass.getRaDec: Messier, NGC or IC number → RA (hours), Dec (degrees)."""
    match = _NAME.match(name)
    if not match:
        raise InputError("Enter a Messier, NGC or IC name, for example m13, NGC5466 or IC4665.")
    prefix, number = match.group(1).lower(), int(match.group(2))
    filename, j0 = _CATALOGS[prefix]
    with open(catalog_dir / filename, encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if not row or int(row[0][j0:]) != number:
                continue
            ra = float(row[1]) + (float(row[2]) + float(row[3]) / 60) / 60
            deg, mins, secs = float(row[4]), float(row[5]), float(row[6])
            if row[4].strip()[0] != "-":
                dec = deg + (mins + secs / 60) / 60
            else:
                dec = deg - (mins + secs / 60) / 60
            return ObjectPosition(f"{prefix.upper()}{number}", ra, dec)
    raise InputError(f"{prefix.upper()}{number} is not in the catalog.")


def format_sexagesimal(angle: float) -> str:
    """ass.decimal_to_min_sec: 16.6948 → '+16 : 41 : 41'."""
    sign = "-" if angle < 0 else "+"
    angle = abs(angle) + 1e-6
    whole = math.floor(angle)
    minutes = math.floor((angle - whole) * 60.0)
    seconds = round(100 * (angle - whole - minutes / 60.0) * 3600.0) / 100
    return f"{sign}{whole:02.0f} : {minutes:02.0f} : {seconds:02.0f}"


def parse_sexagesimal(text: str) -> float:
    """ass.min_sec_to_decimal ('+hh : mm : ss'), with validation."""
    parts = text.split(":")
    try:
        if len(parts) != 3:
            raise ValueError
        head = parts[0].strip()
        sign = -1.0 if head.startswith("-") else 1.0
        whole = int(head.lstrip("+-"))
        minutes = int(parts[1])
        seconds = float(parts[2])
    except ValueError:
        raise InputError(f"'{text}' is not in the form +hh : mm : ss.") from None
    if not (0 <= minutes < 60 and 0 <= seconds < 60):
        raise InputError(f"'{text}' has minutes or seconds out of range.")
    return sign * (whole + (minutes + seconds / 60) / 60)


def validate_query(ra_hours: float, dec_deg: float, fov: float, maxmag: float) -> None:
    values = (ra_hours, dec_deg, fov, maxmag)
    if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in values):
        raise InputError("Right ascension, declination, field of view and magnitude must be numbers.")
    if not 0 <= ra_hours < 24:
        raise InputError("Right ascension must be between 0 and 24 hours.")
    if not -90 <= dec_deg <= 90:
        raise InputError("Declination must be between −90 and +90 degrees.")
    if not FOV_RANGE[0] <= fov <= FOV_RANGE[1]:
        raise InputError(f"Field of view must be {FOV_RANGE[0]}–{FOV_RANGE[1]} degrees.")
    if not MAXMAG_RANGE[0] <= maxmag <= MAXMAG_RANGE[1]:
        raise InputError(f"Maximum magnitude must be {MAXMAG_RANGE[0]}–{MAXMAG_RANGE[1]}.")


def gaia_adql(ra_hours: float, dec_deg: float, fov: float, maxmag: float, service: str) -> str:
    """The query ass.get_star_info sends, for ESA ('esa') or GAVO's copy ('gavo')."""
    validate_query(ra_hours, dec_deg, fov, maxmag)
    where = (
        f"WHERE DISTANCE({ra_hours * 15:14f}, {dec_deg:14f}, gaia.ra, gaia.dec) < {fov / math.sqrt(2):14f} "
        f"AND gaia.phot_g_mean_mag < {maxmag:5f} "
    )
    if service == "esa":
        return (
            "SELECT  gaia.ra, gaia.dec, gaia.parallax, gaia.phot_g_mean_mag, gaia.bp_rp, gaia.pmra, "
            "gaia.pmdec, gaia.phot_g_mean_flux, gaia.phot_rp_mean_flux, gaia.phot_bp_mean_flux "
            "FROM gaiadr3.gaia_source AS gaia " + where
        )
    if service == "gavo":
        return (
            "SELECT gaia.ra, gaia.dec, gaia.parallax, gaia.phot_g_mean_mag, "
            "gaia.phot_bp_mean_mag, gaia.phot_rp_mean_mag, gaia.pmra, gaia.pmdec "
            "FROM gaia.dr3lite AS gaia " + where
        )
    raise ValueError(f"unknown service {service!r}")


# --- Reading Gaia results ------------------------------------------------------------


def parse_gaia_table(data: bytes, filename: str) -> dict[str, np.ndarray]:
    """Read a Gaia result (GAVO CSV, or an ESA download as CSV, VOTable or FITS).

    Returns float64 columns REQUIRED with NaN where a value is missing, in file order.
    """
    if len(data) > MAX_UPLOAD_BYTES:
        raise InputError("The file is larger than 50 MB. Use a smaller field or magnitude limit.")
    if not data.strip():
        raise InputError("The file is empty.")
    lower = filename.lower()
    if lower.endswith((".vot", ".votable", ".xml")) or data.lstrip()[:5] == b"<?xml":
        columns = _read_astropy(data, "votable")
    elif lower.endswith((".fits", ".fit", ".fts")) or data[:6] == b"SIMPLE":
        columns = _read_astropy(data, "fits")
    else:
        columns = _read_csv(data)
    if "bp_rp" not in columns and {"phot_bp_mean_mag", "phot_rp_mean_mag"} <= columns.keys():
        columns["bp_rp"] = columns["phot_bp_mean_mag"] - columns["phot_rp_mean_mag"]
    missing = [name for name in REQUIRED if name not in columns]
    if missing:
        raise InputError("The file is missing these Gaia columns: " + ", ".join(missing) + ".")
    rows = len(columns["ra"])
    if rows > MAX_ROWS:
        raise InputError(f"The file has more than {MAX_ROWS:,} stars.")
    return {name: columns[name] for name in REQUIRED}


def _read_csv(data: bytes) -> dict[str, np.ndarray]:
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise InputError("The file is not a CSV, VOTable or FITS table.") from None
    reader = csv.reader(io.StringIO(text))
    header = [h.strip().lower().removeprefix("gaia.") for h in next(reader, [])]
    wanted = {name: i for i, name in enumerate(header)}
    values: dict[str, list[float]] = {name: [] for name in wanted}
    for line, row in enumerate(reader, start=2):
        if not row:
            continue
        if len(row) != len(header):
            raise InputError(f"Line {line} of the CSV has {len(row)} fields, expected {len(header)}.")
        for name, i in wanted.items():
            values[name].append(_number(row[i]))
    return {name: np.asarray(v, dtype=np.float64) for name, v in values.items()}


def _number(text: str) -> float:
    text = text.strip()
    if text in ("", "--", "null", "NULL", "nan", "NaN"):
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def _read_astropy(data: bytes, fmt: str) -> dict[str, np.ndarray]:
    from astropy.table import Table

    try:
        with warnings.catch_warnings():  # ESA files carry units astropy can't parse; harmless
            warnings.simplefilter("ignore")
            table = Table.read(io.BytesIO(data), format=fmt)
    except Exception:
        raise InputError(f"The file could not be read as {fmt.upper()}.") from None
    columns: dict[str, np.ndarray] = {}
    for name in table.colnames:
        column = table[name]
        key = name.lower().removeprefix("gaia.")
        try:
            values = np.ma.asarray(column, dtype=np.float64).filled(np.nan)
        except (TypeError, ValueError):
            continue
        columns[key] = np.asarray(values, dtype=np.float64)
    return columns


# --- Bright stars (Hipparcos/HYG) ----------------------------------------------------


def bright_star_rows(hyg_csv: Path) -> list[list[float]]:
    """Rows of hygdata_v3.csv that get_star_info can ever use, in file order.

    Columns: ra (hours), dec, mag, absmag, ci, pmra, pmdec. Used by the site build.
    """
    rows: list[list[float]] = []
    with open(hyg_csv, newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if row[7] == "" or row[8] == "" or row[13] == "" or row[16] == "":
                continue
            mag = float(row[13])
            if not -5 < mag < 5.5:
                continue
            rows.append([float(row[i]) for i in (7, 8, 13, 14, 16, 10, 11)])
    return rows


# --- get_star_info -------------------------------------------------------------------


@dataclass
class Stars:
    """Selected stars, as in ass.get_star_info's return (float32, RA in hours)."""

    ra: np.ndarray
    dec: np.ndarray
    mag: np.ndarray
    absmag: np.ndarray
    coloridx: np.ndarray
    parallax: np.ndarray
    pm_ra: np.ndarray
    pm_dec: np.ndarray

    def subset(self, keep: np.ndarray) -> Stars:
        return Stars(**{k: v[keep] for k, v in self.__dict__.items()})


def _radec_to_xyz(ra: Any, dec: Any) -> tuple[Any, Any, Any]:
    # ass.radec_to_xyz; note it treats RA (hours) as degrees, kept for parity.
    ra, dec = np.asarray(ra) * np.pi / 180, np.asarray(dec) * np.pi / 180
    return np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)


def select_stars(
    gaia: dict[str, np.ndarray],
    bright: np.ndarray,
    ra_mid: float,
    dec_mid: float,
    fov: float,
) -> Stars:
    """Everything get_star_info does after the query returns."""
    f32 = np.float32
    present = np.ones(len(gaia["ra"]), dtype=bool)
    for name in REQUIRED:
        present &= ~np.isnan(gaia[name])
    # Divide in float64, then store as float32, exactly as ass.py does (float(r)/15).
    ra = (gaia["ra"][present] / 15).astype(f32)
    dec = gaia["dec"][present].astype(f32)
    parallax = np.abs(gaia["parallax"][present]).astype(f32)
    mag = gaia["phot_g_mean_mag"][present].astype(f32)
    coloridx = gaia["bp_rp"][present].astype(f32)
    pm_ra = (gaia["pmra"][present] / 15).astype(f32)
    pm_dec = gaia["pmdec"][present].astype(f32)
    with np.errstate(divide="ignore", invalid="ignore"):
        absmag = (
            mag.astype(np.float64) - 2.5 * np.log10(100.0**2 / parallax.astype(np.float64) ** 2)
        ).astype(f32)
    keep = (mag >= 4.5) & (coloridx < 3.0)  # NaNs already removed above

    # Hipparcos stars near the field (bright: ra_h, dec, mag, absmag, ci, pmra, pmdec).
    b = np.asarray(bright, dtype=np.float64).reshape(-1, 7)
    b_ra, b_dec = b[:, 0].astype(f32), b[:, 1].astype(f32)
    xc, yc, zc = _radec_to_xyz(ra_mid, dec_mid)
    x1, y1, z1 = _radec_to_xyz(b_ra.astype(np.float64), b_dec.astype(np.float64))
    dot = xc * x1 + yc * y1 + zc * z1
    cross = np.sqrt((yc * z1 - y1 * zc) ** 2 + (zc * x1 - z1 * xc) ** 2 + (xc * y1 - x1 * yc) ** 2)
    with np.errstate(invalid="ignore"):
        angle = 180 * np.arcsin(np.minimum(cross, 1.0)) / np.pi
    near = (dot > 0) & (angle <= fov * math.sqrt(2))

    def join(gaia_values: np.ndarray, bright_values: np.ndarray) -> np.ndarray:
        return np.concatenate([gaia_values[keep], bright_values[near].astype(f32)])

    return Stars(
        ra=join(ra, b_ra),
        dec=join(dec, b_dec),
        mag=join(mag, b[:, 2]),
        absmag=join(absmag, b[:, 3]),
        coloridx=join(coloridx, b[:, 4]),
        parallax=np.concatenate([parallax[keep], np.zeros(int(near.sum()), dtype=f32)]),
        pm_ra=join(pm_ra, b[:, 5] / 15),
        pm_dec=join(pm_dec, b[:, 6]),
    )


# --- show_gaia_image -----------------------------------------------------------------


def _gnomonic(ra: Any, dec: Any, ra_mid: float, dec_mid: float) -> tuple[Any, Any, Any]:
    theta = (np.asarray(ra, dtype=np.float64) - ra_mid) * np.pi / 12
    phi = (90 - np.asarray(dec, dtype=np.float64)) * np.pi / 180
    x0, y0, z0 = np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)
    t = dec_mid * np.pi / 180
    x1 = x0 * np.cos(t) + z0 * np.sin(t)
    z1 = -x0 * np.sin(t) + z0 * np.cos(t)
    with np.errstate(divide="ignore", invalid="ignore"):
        return x1 / np.abs(x1), y0 / np.abs(x1), z1 / np.abs(x1)


@dataclass
class Frame:
    stars: Stars  # stars inside the picture (ass.py's keepit)
    count: int  # "number of plotted stars"
    brightest: float  # "Brightest star" (minmag0)
    pixels: tuple[np.ndarray, np.ndarray]  # (row, col) at ORIGINAL_SIZE, for rendering
    flux_maxmag: float


def frame_stars(stars: Stars, ra_mid: float, dec_mid: float, fov: float) -> Frame:
    """Which stars show_gaia_image draws, computed on its 4000-pixel grid."""
    n = ORIGINAL_SIZE
    _, _, zmax = _gnomonic(ra_mid, dec_mid + fov / 2, ra_mid, dec_mid)
    _, _, zmin = _gnomonic(ra_mid, dec_mid - fov / 2, ra_mid, dec_mid)
    zmid = (zmax + zmin) / 2
    size = zmax - zmid
    ymin, ymax, zmin, zmax = -size, size, zmid - size, zmid + size
    x, y, z = _gnomonic(stars.ra, stars.dec, ra_mid, dec_mid)
    with np.errstate(invalid="ignore"):
        row = n - np.round(n * (z - zmin) / (zmax - zmin))  # j - w
        col = n - np.round(n * (y - ymin) / (ymax - ymin))  # i - w
        keep = (stars.mag < 22) & (x > 0) & (col >= 0) & (col < n) & (row >= 0) & (row < n)
    kept = stars.subset(keep)
    brightest = float(kept.mag.min()) if keep.any() else math.nan
    flux_maxmag = max(float(np.max(stars.mag)), 18.0) if len(stars.mag) else 18.0
    return Frame(kept, int(keep.sum()), brightest, (row[keep], col[keep]), flux_maxmag)


def prepare_psf(raw: np.ndarray) -> np.ndarray:
    """show_gaia_image's PSF from fit/sxv_psf4-RGB.fit data (3, H, W) → (500, 500, 3), peak 1."""
    w = ORIGINAL_PSF_HALF_WIDTH
    raw = np.asarray(raw, dtype=np.float32)
    psf = np.zeros((2 * w, 2 * w, 3), dtype=np.float32)
    for c in range(3):
        psf[:, :, c] = 1 + raw[c, 250 - w : 250 + w, 250 - w : 250 + w]
    clipped = np.maximum(psf - np.percentile(psf, 80), 0)
    return np.asarray(clipped / clipped.max(), dtype=np.float32)


def scale_psf(psf: np.ndarray, size: int) -> np.ndarray:
    """Resample the (2w, 2w, 3) PSF for a render of `size` pixels (bilinear, peak kept)."""
    width = max(2, 2 * round(ORIGINAL_PSF_HALF_WIDTH * size / ORIGINAL_SIZE))
    if width == psf.shape[0]:
        return psf.astype(np.float32)
    src = np.linspace(0, psf.shape[0] - 1, width)
    lo = np.floor(src).astype(int)
    hi = np.minimum(lo + 1, psf.shape[0] - 1)
    frac = (src - lo)[:, None, None]
    rows = psf[lo] * (1 - frac) + psf[hi] * frac
    frac_c = (src - lo)[None, :, None]
    out = rows[:, lo] * (1 - frac_c) + rows[:, hi] * frac_c
    return np.asarray(out / out.max(), dtype=np.float32)


def render_field(frame: Frame, psf: np.ndarray, size: int = RENDER_SIZE) -> np.ndarray:
    """show_gaia_image's 2D picture (A1s), at `size` pixels: log(A + 20) − log(20)."""
    kernel = scale_psf(psf, size)
    w = kernel.shape[0] // 2
    scale = size / ORIGINAL_SIZE
    row = (
        np.round(frame.pixels[0] * scale).astype(int)
        if size != ORIGINAL_SIZE
        else frame.pixels[0].astype(int)
    )
    col = (
        np.round(frame.pixels[1] * scale).astype(int)
        if size != ORIGINAL_SIZE
        else frame.pixels[1].astype(int)
    )
    inside = (row >= 0) & (row < size) & (col >= 0) & (col < size)
    stars = frame.stars.subset(inside)
    row, col = row[inside], col[inside]
    color = stars.coloridx.astype(np.float64) - 0.7
    red = 10 ** (0.5 * color / 2.5)
    flux = 10 ** ((frame.flux_maxmag - stars.mag.astype(np.float64) - 8) / 2.5)
    weights = (red * flux, flux, flux / red)

    padded = size + 2 * w
    fft_size = padded + 2 * w
    image = np.empty((size, size, 3), dtype=np.float32)
    for channel in range(3):
        grid = np.zeros((padded, padded))
        np.add.at(grid, (row + w, col + w), weights[channel])
        spectrum = np.fft.rfft2(grid, s=(fft_size, fft_size))
        spectrum *= np.fft.rfft2(kernel[:, :, channel].astype(np.float64), s=(fft_size, fft_size))
        full = np.fft.irfft2(spectrum, s=(fft_size, fft_size))
        # Patch for a star at (j, i) covers [j - w, j + w): linear conv shifted by w.
        accumulated = full[w : w + padded, w : w + padded][w : size + w, w : size + w]
        image[:, :, channel] = np.log(np.maximum(accumulated, 0) + 20) - np.log(20)
    return image


def display_base(image: np.ndarray) -> tuple[np.ndarray, float, float]:
    """slidershow's preprocessing for an RGB image: (B, pct_lo, pct_hi)."""
    a = image + 0.001 * (np.max(image) - np.min(image))
    b = a / np.max(a)
    dims = [np.percentile(b[:, :, c], 0.01) for c in range(3)]
    low = min(dims)
    for c in range(3):
        b[:, :, c] = b[:, :, c] - dims[c] + low
    positive = b[b > 0]
    pct_lo = float(np.percentile(positive, 0.1)) if positive.size else 0.0
    pct_hi = float(np.percentile(b, 99.9))
    return b.astype(np.float32), pct_lo, pct_hi


# --- GUI_makeHRdiag0 -----------------------------------------------------------------


def hr_arrays(stars: Stars, magnitude_type: str) -> dict[str, np.ndarray]:
    """The scatter data GUI_makeHRdiag0 plots (x, y, sizes, RGBA colours)."""
    if magnitude_type not in ("Apparent", "Absolute"):
        raise InputError("Magnitude type must be Apparent or Absolute.")
    mag0 = (stars.mag if magnitude_type == "Apparent" else stars.absmag).copy()
    coloridx = stars.coloridx.copy()
    n = len(mag0)
    if n == 0:
        raise InputError("No stars to plot. Click Show Stars first, or widen the field.")
    if magnitude_type == "Absolute":
        k = int(np.argmin(mag0))
        mag0[k], coloridx[k] = mag0[n - 1], coloridx[n - 1]
    finite = np.isfinite(mag0)
    maxmag, minmag = float(np.max(mag0[finite])), float(np.min(mag0[finite]))
    clridx = 1.25 * (coloridx - 1.0)
    span = minmag - maxmag
    with np.errstate(divide="ignore", invalid="ignore"):
        diam = 12 * np.sqrt(np.minimum(mag0 - maxmag, 0) / span) if span else np.zeros(n)
        alpha = (mag0 - maxmag) / span if span else np.zeros(n)
    col = 0.6 * clridx
    red = np.clip(1 + col, 0, 1)
    blue = np.clip(1 - col, 0, 1)
    colors = np.stack([red, (red + blue) / 2, blue, np.minimum(alpha, 1)], axis=1).astype(np.float32)
    diam = diam.astype(np.float32)
    diam[0], colors[0] = 0, 0  # ass.py's loop starts at index 1
    keep = finite & np.isfinite(diam)
    return {
        "x": clridx[keep],
        "y": mag0[keep],
        "s": diam[keep],
        "c": colors[keep],
        "limits": np.array([minmag, maxmag], dtype=np.float64),
    }


def hr_png(stars: Stars, magnitude_type: str, title: str, source: str) -> bytes:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    data = hr_arrays(stars, magnitude_type)
    minmag, maxmag = data["limits"]
    fig = Figure(figsize=(10, 8), dpi=100, facecolor="black")
    ax = fig.add_subplot(1, 1, 1, facecolor="black")
    for side in ("bottom", "left", "top", "right"):
        ax.spines[side].set_color("white")
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")
    ax.tick_params(axis="x", colors="white")
    ax.tick_params(axis="y", colors="white")
    ax.scatter(data["x"], data["y"], alpha=0.5, facecolor=data["c"], s=data["s"], edgecolors="none")
    colorrange = min(max(float(np.max(data["x"])), -float(np.min(data["x"]))), 2.0)
    ax.axis((-colorrange, colorrange, minmag - 0.2, maxmag))
    ax.set_xlabel("Color index (B-R)")
    ax.set_ylabel(f"Brightness ({magnitude_type.lower()} magnitude)")
    ax.set_title(f"{title}   Hertzsprung-Russell Diagram using Gaia Data", color="white")
    ax.invert_yaxis()
    fig.text(0.99, 0.01, source, color="0.6", fontsize=8, ha="right", va="bottom")
    out = io.BytesIO()
    fig.savefig(out, format="png", facecolor="black")
    return out.getvalue()
