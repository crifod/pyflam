"""Tests for the shared WindField type and the ESRI ASCII reader."""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyflam
from pyflam import wind
from pyflam.landscape import Landscape
from pyflam.units import mph_to_ft_per_min

MOIST = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
             m_live_herb=0.60, m_live_woody=0.90)


def _write_asc(path, arr, *, cellsize=60.0, xll=0.0, yll=0.0, nodata=-9999.0):
    arr = np.asarray(arr, dtype=float)
    with open(path, "w") as fh:
        fh.write(f"ncols {arr.shape[1]}\nnrows {arr.shape[0]}\n")
        fh.write(f"xllcorner {xll}\nyllcorner {yll}\ncellsize {cellsize}\n")
        fh.write(f"NODATA_value {nodata}\n")
        for row in arr:
            fh.write(" ".join(repr(float(v)) for v in row) + "\n")


# --- ESRI ASCII reader --------------------------------------------------------

def test_read_esri_ascii(tmp_path):
    p = tmp_path / "g.asc"
    _write_asc(p, [[1, 2, 3], [4, 5, 6]], cellsize=30.0, xll=500000.0, yll=4000000.0)
    g = wind.read_esri_ascii(str(p))
    assert g.data.shape == (2, 3)
    assert g.data[0, 0] == 1 and g.data[1, 2] == 6  # row 0 is north
    assert g.cellsize == 30.0
    assert g.west == 500000.0
    assert g.north == 4000000.0 + 2 * 30.0


def test_read_esri_ascii_center_origin(tmp_path):
    p = tmp_path / "c.asc"
    with open(p, "w") as fh:
        fh.write("ncols 2\nnrows 2\nxllcenter 15.0\nyllcenter 15.0\n")
        fh.write("cellsize 30.0\nNODATA_value -9999\n1 2\n3 4\n")
    g = wind.read_esri_ascii(str(p))
    assert g.xllcorner == 0.0 and g.yllcorner == 0.0  # center -> corner


# --- WindField ----------------------------------------------------------------

def test_direction_toward():
    wf = wind.WindField(speed=np.array([[4.0]]), direction=np.array([[270.0]]),
                        cellsize=30.0, west=0.0, north=0.0)
    assert wf.direction_toward[0, 0] == 90.0  # from west -> toward east


def test_speed_unit_conversion():
    wf = wind.WindField(speed=np.array([[1.0]]), direction=np.array([[0.0]]),
                        cellsize=30.0, west=0.0, north=0.0, speed_units="mph")
    assert wf.speed_ft_per_min()[0, 0] == pytest.approx(88.0)
    wf.speed_units = "m/s"
    assert wf.speed_ft_per_min()[0, 0] == pytest.approx(60.0 / 0.3048)
    wf.speed_units = "bogus"
    with pytest.raises(ValueError):
        wf.speed_ft_per_min()


def test_resample_nearest_onto_landscape():
    # 2x2 wind at 60 m over the same extent as a 4x4 landscape at 30 m.
    wf = wind.WindField(
        speed=np.array([[10.0, 20.0], [30.0, 40.0]]),
        direction=np.array([[1.0, 2.0], [3.0, 4.0]]),
        cellsize=60.0, west=0.0, north=120.0, speed_units="mph")
    ls = Landscape(fuel_model=np.ones((4, 4), np.int16),
                   slope=np.zeros((4, 4), np.int16),
                   cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=120.0)
    r = wf.to_landscape(ls)
    expected = np.array([[10, 10, 20, 20], [10, 10, 20, 20],
                         [30, 30, 40, 40], [30, 30, 40, 40]], dtype=float)
    assert np.array_equal(r.speed, expected)
    assert r.direction[2, 3] == 4.0


def test_to_midflame_applies_reduction():
    wf = wind.WindField(speed=np.array([[10.0]]), direction=np.array([[0.0]]),
                        cellsize=30.0, west=0.0, north=0.0, speed_units="mph")
    mid = wf.to_midflame(wind_reduction_factor=0.4)
    assert mid[0, 0] == pytest.approx(mph_to_ft_per_min(10) * 0.4)


def test_windfield_drives_basic_fire_behavior():
    fuel = np.array([[1, 102], [165, 1]], dtype=np.int16)
    slope = np.array([[0, 10], [20, 30]], dtype=np.int16)
    ls = Landscape(fuel_model=fuel, slope=slope, cellsize_x=30.0,
                   cellsize_y=30.0, west=0.0, north=60.0)
    wf = wind.WindField(speed=np.array([[3.0, 8.0], [5.0, 12.0]]),
                        direction=np.zeros((2, 2)), cellsize=30.0,
                        west=0.0, north=60.0, speed_units="mph")
    mid = wf.to_midflame(ls, wind_reduction_factor=1.0)
    res = pyflam.basic_fire_behavior(ls, wind_midflame=mid, **MOIST)
    for r in range(2):
        for c in range(2):
            fb = pyflam.spread(
                pyflam.get_fuel_model(int(fuel[r, c])),
                wind_midflame=mph_to_ft_per_min(wf.speed[r, c]),
                slope=math.tan(math.radians(slope[r, c])), **MOIST)
            assert res["rate_of_spread"][r, c] == pytest.approx(fb.rate_of_spread)
