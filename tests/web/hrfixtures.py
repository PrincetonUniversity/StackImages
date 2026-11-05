"""Synthetic star-cluster FITS images for the HR-diagram tests (web and desktop).

Deterministic (seeded) and generated on the fly, so no image is committed. They are made
to exercise every rule in makeHRdiag_new: star widths that findstars' FWHM test keeps and
rejects, hot pixels, colours its red/green and green/blue tests drop, and a spread of
brightness wide enough for the magnitude boxes to matter.
"""

from __future__ import annotations

import io

import numpy as np


def cluster_rgb(height: int = 600, width: int = 800, stars: int = 900, seed: int = 2) -> np.ndarray:
    """A (3, height, width) float32 RGB image of a centrally concentrated cluster."""
    rng = np.random.default_rng(seed)
    image = rng.normal(100, 2, (3, height, width)).astype(np.float32)
    yy, xx = np.mgrid[-7:8, -7:8]
    radius = np.abs(rng.normal(0, min(height, width) / 6, stars))
    angle = rng.uniform(0, 2 * np.pi, stars)
    ys = np.clip((height / 2 + radius * np.sin(angle)).astype(int), 30, height - 31)
    xs = np.clip((width / 2 + radius * np.cos(angle)).astype(int), 30, width - 31)
    for y, x in zip(ys, xs, strict=True):
        sigma = rng.choice([1.6, 1.6, 1.6, 1.6, 2.6, 1.0])  # mostly one width; some too wide/narrow
        psf = np.exp(-(xx**2 + yy**2) / (2 * sigma**2)).astype(np.float32)
        colour = rng.uniform(0.4, 1.6, 3)
        if rng.uniform() < 0.05:
            colour[2] *= 0.05  # very little blue: G/B above e^2
        elif rng.uniform() < 0.05:
            colour[0] *= 0.05  # very little red: |log(R/G)| above 2
        peak = rng.uniform(20, 6000)
        image[:, y - 7 : y + 8, x - 7 : x + 8] += (peak * colour)[:, None, None].astype(np.float32) * psf
    hot = rng.integers(30, [height - 30, width - 30], (40, 2))
    image[:, hot[:, 0], hot[:, 1]] += 5000  # single hot pixels
    return image


def fits_bytes(data: np.ndarray, header: dict[str, object] | None = None, extension: bool = False) -> bytes:
    """`data` as a FITS file, in the primary HDU or (extension=True) the first extension."""
    from astropy.io import fits

    primary = fits.PrimaryHDU(None if extension else data)
    for key, value in (header or {}).items():
        primary.header[key] = value
    hdus = [primary] + ([fits.ImageHDU(data)] if extension else [])
    out = io.BytesIO()
    fits.HDUList(hdus).writeto(out)
    return out.getvalue()
