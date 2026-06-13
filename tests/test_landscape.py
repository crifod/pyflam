"""Tests for the Landscape abstraction and whole-landscape fire behavior."""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyflam
from pyflam import basic_fire_behavior
from pyflam.landscape import Landscape
from pyflam.units import mph_to_ft_per_min

MOIST = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
             m_live_herb=0.60, m_live_woody=0.90)
WIND = mph_to_ft_per_min(5)


def _landscape(fuel, slope, **kw):
    return Landscape(
        fuel_model=np.asarray(fuel, dtype=np.int16),
        slope=np.asarray(slope, dtype=np.int16),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=0.0,
        slope_units=kw.pop("slope_units", "degrees"), **kw,
    )


def test_vectorized_matches_scalar():
    """Every cell must equal the scalar Rothermel result for its fuel + slope."""
    fuel = [[1, 102, 165], [10, 145, 183]]
    slope = [[0, 15, 30], [45, 5, 20]]
    ls = _landscape(fuel, slope)
    res = basic_fire_behavior(ls, wind_midflame=WIND, **MOIST)

    for r in range(2):
        for c in range(3):
            fb = pyflam.spread(
                pyflam.get_fuel_model(fuel[r][c]),
                wind_midflame=WIND,
                slope=math.tan(math.radians(slope[r][c])),
                **MOIST,
            )
            assert res["rate_of_spread"][r, c] == pytest.approx(fb.rate_of_spread)
            assert res["flame_length"][r, c] == pytest.approx(fb.flame_length)
            assert res["fireline_intensity"][r, c] == pytest.approx(
                fb.fireline_intensity
            )


def test_nonburnable_cells_are_zero():
    ls = _landscape([[91, 92, 93, 98, 99]], [[10, 10, 10, 10, 10]])
    res = basic_fire_behavior(ls, wind_midflame=WIND, **MOIST)
    assert np.all(res["rate_of_spread"] == 0.0)
    assert np.all(res["flame_length"] == 0.0)


def test_unknown_fuel_code_is_nodata():
    ls = _landscape([[1, -9999, 300]], [[0, 0, 0]])
    res = basic_fire_behavior(ls, wind_midflame=WIND, nodata=np.nan, **MOIST)
    assert math.isfinite(res["rate_of_spread"][0, 0])     # FM1
    assert math.isnan(res["rate_of_spread"][0, 1])        # -9999 fill
    assert math.isnan(res["rate_of_spread"][0, 2])        # undefined code


def test_slope_increases_spread_across_grid():
    ls = _landscape([[1, 1, 1, 1]], [[0, 15, 30, 45]])
    ros = basic_fire_behavior(ls, wind_midflame=WIND, **MOIST)["rate_of_spread"]
    assert list(ros[0]) == sorted(ros[0])
    assert ros[0, 0] < ros[0, -1]


def test_slope_units_percent_vs_degrees():
    deg = _landscape([[1]], [[45]], slope_units="degrees")
    pct = _landscape([[1]], [[100]], slope_units="percent")  # 100% == 45 deg
    a = basic_fire_behavior(deg, wind_midflame=WIND, **MOIST)["rate_of_spread"][0, 0]
    b = basic_fire_behavior(pct, wind_midflame=WIND, **MOIST)["rate_of_spread"][0, 0]
    assert a == pytest.approx(b)


def test_geotransform():
    ls = _landscape([[1, 1]], [[0, 0]])
    ls.west, ls.north = 600000.0, 4900000.0
    assert ls.transform_gdal == (600000.0, 30.0, 0.0, 4900000.0, 0.0, -30.0)


# --- LCP integration ----------------------------------------------------------

def test_landscape_lcp_roundtrip(tmp_path):
    fuel = np.array([[1, 102, 91], [165, 145, 1]], dtype=np.int16)
    slope = np.array([[0, 10, 20], [30, 40, 5]], dtype=np.int16)
    ls = _landscape(fuel, slope, elevation=np.full((2, 3), 1500, np.int16),
                    aspect=np.full((2, 3), 180, np.int16),
                    canopy_cover=np.full((2, 3), 40, np.int16))
    path = tmp_path / "ls.lcp"
    ls.to_lcp(str(path))
    back = Landscape.from_lcp(str(path))
    assert np.array_equal(back.fuel_model, fuel)
    assert np.array_equal(back.slope, slope)
    assert back.cellsize_x == 30.0

    # Behavior is identical whether computed before or after the LCP round-trip.
    a = basic_fire_behavior(ls, wind_midflame=WIND, **MOIST)["rate_of_spread"]
    b = basic_fire_behavior(back, wind_midflame=WIND, **MOIST)["rate_of_spread"]
    assert np.array_equal(a, b)


# --- GeoTIFF integration ------------------------------------------------------

def test_geotiff_roundtrip_and_output(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import Affine

    # Use a CRS if this environment's PROJ can build one; otherwise go CRS-less.
    # (Some installs have a mismatched proj.db that can't resolve EPSG codes.)
    try:
        crs = rasterio.crs.CRS.from_epsg(5070)
    except Exception:
        crs = None

    fuel = np.array([[1, 102], [165, 1]], dtype=np.int16)
    slope = np.array([[0, 20], [40, 10]], dtype=np.float32)
    transform = Affine(30.0, 0.0, 500000.0, 0.0, -30.0, 5000000.0)

    fuel_tif = tmp_path / "fuel.tif"
    slope_tif = tmp_path / "slope.tif"
    for path, arr, dt in ((fuel_tif, fuel, "int16"), (slope_tif, slope, "float32")):
        with rasterio.open(
            str(path), "w", driver="GTiff", height=2, width=2, count=1,
            dtype=dt, transform=transform, crs=crs,
        ) as ds:
            ds.write(arr, 1)

    ls = Landscape.from_geotiffs(
        {"fuel_model": str(fuel_tif), "slope": str(slope_tif)}
    )
    assert np.array_equal(ls.fuel_model, fuel)
    assert ls.cellsize_x == 30.0 and ls.west == 500000.0
    if crs is not None:
        assert ls.crs.to_epsg() == 5070

    res = basic_fire_behavior(ls, wind_midflame=WIND, **MOIST)
    out = tmp_path / "ros.tif"
    ls.to_geotiff(str(out), res["rate_of_spread"])
    with rasterio.open(str(out)) as ds:
        assert ds.shape == ls.shape
        assert np.allclose(ds.read(1), res["rate_of_spread"], equal_nan=True)
