"""Unit tests for the pyflam_gui core helpers that carry no Streamlit dependency.

Covers AOI bbox geometry, the output-folder convention + zip bundler, and the
refactored pyroconvection classifier (shared with tests/pyroconv_daily.py),
checked against the reference cases in test_pyroconvection_type.py.
"""

from __future__ import annotations

import os
import sys
import zipfile

import numpy as np
import pytest

# Make pyflam_gui importable when running with PYTHONPATH=src.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyflam_gui.core import aoi as aoi_mod
from pyflam_gui.core import outputs
from pyflam_gui.core import pyroconv as pc
from pyflam_gui.core.state import AreaOfInterest


# --- AOI geometry -------------------------------------------------------------

def test_area_of_interest_bbox_and_center():
    a = AreaOfInterest(north=44.6, west=9.6, south=42.2, east=12.5, label="t")
    assert a.bbox == (44.6, 9.6, 42.2, 12.5)
    assert a.lonlat_bbox == (9.6, 12.5, 42.2, 44.6)
    clat, clon = a.center
    assert clat == pytest.approx(43.4)
    assert clon == pytest.approx(11.05)


def test_preset_tuscany_matches_daily_box():
    a = aoi_mod.preset_aoi("Tuscany")
    # tests/pyroconv_daily.py uses LON0,LON1,LAT0,LAT1 = 9.6,12.5,42.2,44.6
    assert a.lonlat_bbox == (9.6, 12.5, 42.2, 44.6)


def test_valid_bbox():
    assert aoi_mod.valid_bbox(44.6, 9.6, 42.2, 12.5)
    assert not aoi_mod.valid_bbox(42.2, 9.6, 44.6, 12.5)   # north <= south
    assert not aoi_mod.valid_bbox(44.6, 12.5, 42.2, 9.6)   # east <= west


def test_bbox_from_geojson_polygon():
    geom = {"type": "Polygon", "coordinates": [[
        [9.6, 42.2], [12.5, 42.2], [12.5, 44.6], [9.6, 44.6], [9.6, 42.2]]]}
    n, w, s, e = aoi_mod.bbox_from_geojson(geom)
    assert (n, w, s, e) == (44.6, 9.6, 42.2, 12.5)


def test_clip_bbox():
    inner = (50.0, 5.0, 40.0, 15.0)
    outer = (44.6, 9.6, 42.2, 12.5)
    assert aoi_mod.clip_bbox(inner, outer) == (44.6, 9.6, 42.2, 12.5)


# --- output folder ------------------------------------------------------------

def test_make_run_dir_and_zip(tmp_path):
    run_dir = outputs.make_run_dir(str(tmp_path), "preview")
    assert os.path.isdir(run_dir)
    assert os.path.basename(run_dir).startswith("preview_")
    with open(os.path.join(run_dir, "a.txt"), "w") as fh:
        fh.write("hello")
    data = outputs.zip_dir(run_dir)
    import io
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        assert "a.txt" in zf.namelist()


# --- pyroconvection classifier ------------------------------------------------

def _levels_from_gamma(gamma_700_500):
    """Build a 700/500 hPa temperature pair giving a target gamma_theta (K/m)."""
    from pyflam.atmosphere import theta_kelvin
    # theta at 700 fixed; solve theta500 = theta700 + gamma*(Z500 - Z700), then
    # back out T from theta. theta = (T+273.15)*(1000/p)^0.286.
    Z = pc.STD_LEVEL_HEIGHT_M
    th700 = theta_kelvin(5.0, 700.0)
    th500 = th700 + gamma_700_500 * (Z[500] - Z[700])
    t500 = th500 * (500.0 / 1000.0) ** 0.286 - 273.15
    return 5.0, t500


def test_classify_fuel_gate_forces_surface_plume():
    # Any atmosphere, but fireline intensity below the 10 MW/m gate -> class 0.
    T2m = np.array([[35.0]])
    RH = np.array([[10.0]])
    t700, t500 = _levels_from_gamma(3.0e-3)
    Tl = {850: np.array([[24.0]]), 700: np.array([[t700]]), 500: np.array([[t500]])}
    cats = pc.classify(T2m, RH, Tl, fli=np.array([[5.0e3]]))
    assert cats[0, 0] == 0   # surface_plume


def test_classify_runs_and_returns_levels_in_range():
    T2m = np.array([[35.0, 20.0]])
    RH = np.array([[10.0, 60.0]])
    t700, t500 = _levels_from_gamma(3.0e-3)
    Tl = {850: np.array([[24.0, 18.0]]),
          700: np.array([[t700, t700]]),
          500: np.array([[t500, t500]])}
    cats = pc.classify(T2m, RH, Tl)
    assert cats.dtype == np.int16
    assert cats.min() >= 0 and cats.max() <= 4


# --- maps / plotting helpers (no Streamlit runtime needed) --------------------

def test_array_to_png_datauri_and_transparency():
    from pyflam_gui.core import maps
    a = np.array([[0.0, 1.0], [np.inf, 0.5]])      # inf -> transparent
    uri = maps._array_to_png_datauri(a, cmap="magma")
    assert uri.startswith("data:image/png;base64,")
    assert len(uri) > 50


def test_wgs84_bounds_none_without_crs():
    from pyflam_gui.core import maps, landscape as lsc
    ls = lsc.synthetic_landscape(n=20, cellsize=30)   # no CRS
    assert maps._wgs84_bounds(ls) is None


def test_fire_size_hist_returns_figure():
    from pyflam_gui.core import plotting
    fig = plotting.fire_size_hist(np.array([1e4, 5e4, 0.0, np.nan, 2e5]), to_ha=1e-4)
    assert fig is not None
