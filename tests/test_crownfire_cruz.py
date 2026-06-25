"""Tests for the Cruz et al. (2004/2005) crown-fire models.

The Cruz, Alexander & Wakimoto active crown-fire rate of spread (2005) and the
logistic crown-fire-initiation model (2004) — the empirical successors to the
Rothermel (1991) / Van Wagner stack, which Cruz & Alexander (2010) showed
under-predicts. Coefficients are checked against the published equations.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyflam
from pyflam import crownfire
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)


def _surface(n, wind_mph, slope=0.0):
    return pyflam.spread(pyflam.get_fuel_model(n),
                         wind_midflame=mph_to_ft_per_min(wind_mph), slope=slope,
                         **SCENARIO)


# --- Cruz (2005) active crown ROS ---------------------------------------------

def _u20ft_for_u10(u10_kmh):
    """20-ft wind (ft/min) that maps to a given 10-m open wind (km/h)."""
    return u10_kmh * 0.87 / 0.018288


def test_cruz_ros_matches_published_formula():
    u10 = 30.0
    r = crownfire.active_crown_ros_cruz(
        _u20ft_for_u10(u10), canopy_bulk_density=0.15, m_1h=0.06)
    expect = 11.02 * u10 ** 0.90 * 0.15 ** 0.19 * math.exp(-0.17 * 6.0)
    assert r == pytest.approx(expect, rel=1e-3)


def test_cruz_ros_zero_without_canopy():
    assert crownfire.active_crown_ros_cruz(
        _u20ft_for_u10(30), canopy_bulk_density=0.0, m_1h=0.06) == 0.0


def test_cruz_ros_monotonic():
    base = crownfire.active_crown_ros_cruz(_u20ft_for_u10(20), canopy_bulk_density=0.15, m_1h=0.08)
    assert crownfire.active_crown_ros_cruz(_u20ft_for_u10(35), canopy_bulk_density=0.15, m_1h=0.08) > base
    assert crownfire.active_crown_ros_cruz(_u20ft_for_u10(20), canopy_bulk_density=0.30, m_1h=0.08) > base
    assert crownfire.active_crown_ros_cruz(_u20ft_for_u10(20), canopy_bulk_density=0.15, m_1h=0.04) > base


# --- Cruz (2004) logistic crown-fire initiation -------------------------------

def test_logistic_matches_published_coefficients():
    u10, fsg, sfc, effm = 25.0, 4.0, 2.5, 6.0      # SFC>=2 -> reference category (0)
    g = 4.236 + 0.357 * u10 - 0.710 * fsg + 0.0 - 0.331 * effm
    expect = 1.0 / (1.0 + math.exp(-g))
    p = crownfire.crown_fire_probability(
        wind_10m_kmh=u10, fuel_strata_gap=fsg, surface_fuel_consumption=sfc,
        fine_fuel_moisture=effm)
    assert p == pytest.approx(expect, rel=1e-6)


def test_logistic_bounds_and_monotonicity():
    def P(**kw):
        base = dict(wind_10m_kmh=20.0, fuel_strata_gap=4.0,
                    surface_fuel_consumption=2.5, fine_fuel_moisture=8.0)
        base.update(kw)
        return crownfire.crown_fire_probability(**base)
    assert 0.0 <= P() <= 1.0
    assert P(wind_10m_kmh=40) > P(wind_10m_kmh=5)        # more wind -> more likely
    assert P(fuel_strata_gap=2) > P(fuel_strata_gap=12)  # smaller gap -> more likely
    assert P(fine_fuel_moisture=4) > P(fine_fuel_moisture=20)


def test_logistic_sfc_category_effect():
    # Lower surface fuel consumption suppresses crowning (negative dummies).
    def P(sfc):
        return crownfire.crown_fire_probability(
            wind_10m_kmh=20.0, fuel_strata_gap=4.0,
            surface_fuel_consumption=sfc, fine_fuel_moisture=8.0)
    assert P(0.5) < P(1.5) < P(2.5)


def test_logistic_vectorized():
    u = np.array([5.0, 40.0])
    p = crownfire.crown_fire_probability(
        wind_10m_kmh=u, fuel_strata_gap=4.0, surface_fuel_consumption=2.5,
        fine_fuel_moisture=8.0)
    assert p.shape == (2,) and p[1] > p[0]


# --- crown_spread selection (point + landscape) -------------------------------

def test_crown_spread_default_unchanged():
    """Default must still be the Rothermel-1991 path (backward compatible)."""
    fb = _surface(10, 25)
    cf = pyflam.crown_fire_behavior(
        fb, canopy_base_height=1.5, canopy_bulk_density=0.18, foliar_moisture=100,
        wind_20ft_ft_per_min=mph_to_ft_per_min(30), **SCENARIO)
    expect = crownfire.active_crown_ros(mph_to_ft_per_min(30), **SCENARIO)
    assert cf.active_rate_of_spread == pytest.approx(expect)


def test_cruz_path_uses_full_active_rate_no_reduction():
    fb = _surface(10, 30)
    cf = pyflam.crown_fire_behavior(
        fb, canopy_base_height=1.0, canopy_bulk_density=0.20, foliar_moisture=100,
        wind_20ft_ft_per_min=mph_to_ft_per_min(35), crown_spread="cruz2005",
        **SCENARIO)
    if cf.fire_type == "active":
        assert cf.rate_of_spread == pytest.approx(cf.active_rate_of_spread)
        assert cf.crown_fraction_burned == 1.0


def test_invalid_crown_spread_raises():
    fb = _surface(10, 20)
    with pytest.raises(ValueError):
        pyflam.crown_fire_behavior(
            fb, canopy_base_height=1.5, canopy_bulk_density=0.15, foliar_moisture=100,
            wind_20ft_ft_per_min=mph_to_ft_per_min(20), crown_spread="nope", **SCENARIO)


def _canopy_landscape(n=20):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.full((n, n), 20.0, dtype=float), aspect=np.full((n, n), 180.0),
        canopy_base_height=np.full((n, n), 10.0, dtype=float),    # *10 -> 1 m
        canopy_bulk_density=np.full((n, n), 20.0, dtype=float),   # *100 -> 0.20
        canopy_height=np.full((n, n), 150.0, dtype=float),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def test_landscape_cruz_vs_rothermel_active_cells():
    ls = _canopy_landscape()
    kw = dict(foliar_moisture=100.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
              wind_midflame=mph_to_ft_per_min(12), **SCENARIO)
    roth = pyflam.crownfire.crown_fire_potential(ls, crown_spread="rothermel1991", **kw)
    cruz = pyflam.crownfire.crown_fire_potential(ls, crown_spread="cruz2005", **kw)
    assert set(cruz) == set(roth)
    active = (roth["fire_type"] == 2) & (cruz["fire_type"] == 2)
    if active.any():
        # Cruz removes the under-predicting CFB reduction -> faster on active cells.
        assert cruz["rate_of_spread"][active].mean() >= roth["rate_of_spread"][active].mean()
