"""Tests for the wind adjustment factor (Albini & Baughman 1979 / Andrews 2012)."""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import wind_reduction as wr


# --- Unsheltered WAF ----------------------------------------------------------

def test_unsheltered_known_values():
    # Hand-computed from WAF = 1.83 / ln((20+0.36H)/(0.13H)).
    assert float(wr.unsheltered_waf(1.0)) == pytest.approx(0.362, abs=0.002)
    assert float(wr.unsheltered_waf(0.4)) == pytest.approx(0.307, abs=0.002)
    assert float(wr.unsheltered_waf(6.0)) == pytest.approx(0.547, abs=0.002)


def test_unsheltered_increases_with_depth():
    depths = [0.2, 0.5, 1.0, 2.5, 6.0]
    wafs = [float(wr.unsheltered_waf(d)) for d in depths]
    assert wafs == sorted(wafs)              # deeper fuel -> less reduction
    assert all(0.0 < w <= 1.0 for w in wafs)


def test_unsheltered_array_and_zero_depth():
    out = wr.unsheltered_waf(np.array([0.0, 1.0, 6.0]))
    assert out[0] == 1.0                     # zero depth -> no reduction (fallback)
    assert out[1] == pytest.approx(0.362, abs=0.002)


# --- Sheltered WAF ------------------------------------------------------------

def test_sheltered_below_unsheltered_under_dense_canopy():
    depth = 1.0
    unshel = float(wr.unsheltered_waf(depth))
    shel = float(wr.sheltered_waf(30.0, 0.6))   # 30 ft canopy, 60% cover
    assert shel < unshel
    assert 0.0 < shel <= 1.0


def test_sheltered_zero_height_is_one():
    assert float(wr.sheltered_waf(0.0, 0.5)) == 1.0


# --- Decision logic -----------------------------------------------------------

def test_chooses_unsheltered_without_canopy_height():
    # No canopy height -> always the per-fuel-depth unsheltered value.
    waf = wr.wind_adjustment_factor(
        np.array([0.4, 6.0]), canopy_height_ft=0.0, canopy_cover_fraction=0.8)
    assert np.allclose(waf, [wr.unsheltered_waf(0.4), wr.unsheltered_waf(6.0)])


def test_chooses_sheltered_when_cover_and_height_present():
    waf = wr.wind_adjustment_factor(
        1.0, canopy_height_ft=40.0, canopy_cover_fraction=0.5)
    assert float(waf) == pytest.approx(float(wr.sheltered_waf(40.0, 0.5)))


def test_low_cover_stays_unsheltered():
    waf = wr.wind_adjustment_factor(
        1.0, canopy_height_ft=40.0, canopy_cover_fraction=0.02)  # below 5%
    assert float(waf) == pytest.approx(float(wr.unsheltered_waf(1.0)))


# --- Landscape fields ---------------------------------------------------------

def _grid(fuel):
    n = 6
    return pyflam.Landscape(
        fuel_model=np.full((n, n), fuel, dtype=int), slope=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        canopy_cover=np.zeros((n, n)), slope_units="degrees")


def test_waf_field_matches_fuel_depth():
    # GR1 (depth 0.4) vs SH5 (depth 6.0): WAF differs and tracks depth.
    gr1 = pyflam.waf_field(_grid(101))
    sh5 = pyflam.waf_field(_grid(145))
    assert np.allclose(gr1, wr.unsheltered_waf(0.4))
    assert np.allclose(sh5, wr.unsheltered_waf(6.0))
    assert sh5.mean() > gr1.mean()


def test_midflame_field_scales_wind():
    ls = _grid(102)                          # GR2, depth 1.0
    mid = pyflam.midflame_field(ls, 1000.0)  # 1000 ft/min 20-ft wind
    assert np.allclose(mid, 1000.0 * wr.unsheltered_waf(1.0))
