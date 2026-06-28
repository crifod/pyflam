"""Landscape loading for the GUI: .lcp / GeoTIFF bundle + a synthetic fallback.

Thin wrappers over :meth:`pyflam.Landscape.from_lcp` and
:meth:`pyflam.Landscape.from_geotiffs`, cached so a re-run of a page does not
re-read the raster. ``synthetic_landscape`` builds a small in-memory landscape so
the whole GUI can be exercised offline with no data files.
"""

from __future__ import annotations

import streamlit as st


@st.cache_resource(show_spinner="Reading landscape...")
def load_lcp(path: str):
    """Read a FlamMap ``.lcp`` into a :class:`pyflam.Landscape` (cached by path)."""
    import pyflam
    return pyflam.Landscape.from_lcp(path)


@st.cache_resource(show_spinner="Reading GeoTIFF bands...")
def load_geotiffs(paths: tuple, slope_units: str = "degrees"):
    """Read a band->path GeoTIFF bundle (``paths`` is a tuple of (name, path) pairs)."""
    import pyflam
    return pyflam.Landscape.from_geotiffs(dict(paths), slope_units=slope_units)


def synthetic_landscape(n: int = 120, cellsize: float = 30.0, seed: int = 0):
    """A small conifer-hill landscape with the full canopy stack, for offline use.

    Mirrors ``tests/make_synthetic_canopy_lcp.make_synthetic_canopy_landscape`` but
    is self-contained (no dependency on the tests package).
    """
    import numpy as np
    import pyflam

    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    rng = np.random.default_rng(seed)
    cx, cy = n * 0.45, n * 0.55
    elev = 800.0 + 600.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * (n * 0.3) ** 2))
    gy, gx = np.gradient(elev, cellsize)
    slope = np.clip(np.degrees(np.arctan(np.hypot(gx, gy))), 0.0, 45.0)
    aspect = np.degrees(np.arctan2(gx, -gy)) % 360.0

    fuel = np.full((n, n), 10, dtype=int)
    fuel[:, n // 2 - 1:n // 2 + 1] = 91
    fuel[: n // 4, : n // 4] = 1
    forest = fuel == 10
    cover = np.where(forest, rng.integers(55, 95, (n, n)), 0)
    cbh = np.clip(15 + 0.9 * yy + rng.normal(0, 8, (n, n)), 5, 140)
    cbd = np.clip(8 + 30 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * (n * 0.35) ** 2))
                  + rng.normal(0, 4, (n, n)), 2, 40)
    ch = np.clip(130 + 0.5 * yy + rng.normal(0, 15, (n, n)), 90, 250)
    z = np.zeros((n, n))
    return pyflam.Landscape(
        fuel_model=fuel, slope=slope.astype(int), aspect=aspect.astype(int),
        elevation=elev.astype(int),
        canopy_cover=np.where(forest, cover, 0).astype(int),
        canopy_base_height=np.where(forest, cbh, z).astype(int),
        canopy_bulk_density=np.where(forest, cbd, z).astype(int),
        canopy_height=np.where(forest, ch, z).astype(int),
        cellsize_x=cellsize, cellsize_y=cellsize, west=0.0, north=n * cellsize,
        slope_units="degrees")


def describe(ls) -> dict:
    """Human-readable summary of a landscape for the status panel."""
    import numpy as np
    nrows, ncols = ls.shape
    fuels = np.unique(np.asarray(ls.fuel_model))
    return {
        "shape": f"{nrows} x {ncols}",
        "cellsize": f"{ls.cellsize_x:g} x {ls.cellsize_y:g}",
        "crs": str(ls.crs) if ls.crs is not None else "(none)",
        "fuel models": ", ".join(str(int(f)) for f in fuels[:12])
                       + (" ..." if fuels.size > 12 else ""),
        "has canopy": ls.canopy_bulk_density is not None,
    }
