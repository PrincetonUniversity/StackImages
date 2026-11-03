"""Shared fixtures for the web tests."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "web"))
sys.path.insert(0, str(Path(__file__).parent))

import showclusters as sc  # noqa: E402


@pytest.fixture(scope="session")
def session() -> sc.Session:
    """A browser session as the worker builds it (catalogs, bright stars, scaled PSF)."""
    from astropy.io import fits

    bright = np.array(sc.bright_star_rows(ROOT / "csv" / "hygdata_v3.csv"))
    psf = sc.scale_psf(sc.prepare_psf(fits.getdata(ROOT / "fit" / "sxv_psf4-RGB.fit")), sc.RENDER_SIZE)
    return sc.Session(ROOT / "csv", bright, psf)
