"""Tests for the crown-fire module (Van Wagner / Rothermel / Scott & Reinhardt).

Two layers, as elsewhere in pyflam: physics/property tests asserting the
relationships the models must obey, and a small golden-master lock on current
numeric outputs. The true acceptance test is a diff against FlamMap's crown-fire
outputs / NEXUS; see tests/REFERENCE.md.
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
    return pyflam.spread(
        pyflam.get_fuel_model(n),
        wind_midflame=mph_to_ft_per_min(wind_mph), slope=slope, **SCENARIO,
    )


# --- Van Wagner critical thresholds -------------------------------------------

def test_critical_intensity_increases_with_cbh_and_fmc():
    base = crownfire.critical_fireline_intensity(2.0, 100.0)
    assert crownfire.critical_fireline_intensity(4.0, 100.0) > base
    assert crownfire.critical_fireline_intensity(2.0, 120.0) > base


def test_critical_intensity_zero_without_canopy():
    assert crownfire.critical_fireline_intensity(0.0, 100.0) == 0.0


def test_critical_active_ros_inverse_in_cbd():
    assert crownfire.critical_active_ros(0.10) == pytest.approx(30.0)
    assert crownfire.critical_active_ros(0.20) == pytest.approx(15.0)
    assert math.isinf(crownfire.critical_active_ros(0.0))


def test_active_ros_increases_with_wind():
    lo = crownfire.active_crown_ros(mph_to_ft_per_min(10), **SCENARIO)
    hi = crownfire.active_crown_ros(mph_to_ft_per_min(30), **SCENARIO)
    assert hi > lo > 0.0


# --- Fire-type classification (Scott & Reinhardt) -----------------------------

def test_low_wind_stays_surface():
    fb = _surface(8, wind_mph=1)  # closed timber litter, low intensity
    cf = pyflam.crown_fire_behavior(
        fb, canopy_base_height=3.0, canopy_bulk_density=0.10,
        foliar_moisture=100, wind_20ft_ft_per_min=mph_to_ft_per_min(1),
        **SCENARIO,
    )
    assert cf.fire_type == "surface"
    assert cf.crown_fraction_burned == 0.0
    assert cf.rate_of_spread == pytest.approx(cf.surface_rate_of_spread)


def test_high_wind_initiates_crown_fire():
    # Drive surface and crown with the same 20-ft wind (midflame = 0.4 * U20),
    # so the two spread rates are consistent.
    u20 = 25
    fb = _surface(10, wind_mph=0.4 * u20)
    cf = pyflam.crown_fire_behavior(
        fb, canopy_base_height=1.5, canopy_bulk_density=0.20,
        foliar_moisture=90, wind_20ft_ft_per_min=mph_to_ft_per_min(u20),
        **SCENARIO,
    )
    assert cf.initiates
    assert cf.fire_type in ("passive", "active")
    assert 0.0 < cf.crown_fraction_burned <= 1.0
    # Final spread is a convex blend of the surface and active crown rates.
    lo = min(cf.surface_rate_of_spread, cf.active_rate_of_spread)
    hi = max(cf.surface_rate_of_spread, cf.active_rate_of_spread)
    assert lo - 1e-9 <= cf.rate_of_spread <= hi + 1e-9


def test_denser_canopy_more_likely_active():
    fb = _surface(10, wind_mph=30)
    sparse = pyflam.crown_fire_behavior(
        fb, canopy_base_height=1.0, canopy_bulk_density=0.05,
        foliar_moisture=90, wind_20ft_ft_per_min=mph_to_ft_per_min(30), **SCENARIO)
    dense = pyflam.crown_fire_behavior(
        fb, canopy_base_height=1.0, canopy_bulk_density=0.30,
        foliar_moisture=90, wind_20ft_ft_per_min=mph_to_ft_per_min(30), **SCENARIO)
    # Denser canopy has a lower critical active ROS, so crowns more readily.
    order = {"surface": 0, "passive": 1, "active": 2}
    assert order[dense.fire_type] >= order[sparse.fire_type]


def test_cfb_monotonic_in_surface_intensity():
    cfbs = []
    for wind in (5, 10, 15, 20, 30):
        fb = _surface(10, wind_mph=wind)
        cf = pyflam.crown_fire_behavior(
            fb, canopy_base_height=1.5, canopy_bulk_density=0.15,
            foliar_moisture=100, wind_20ft_ft_per_min=mph_to_ft_per_min(wind),
            **SCENARIO)
        cfbs.append(cf.crown_fraction_burned)
    assert cfbs == sorted(cfbs)
    assert all(0.0 <= x <= 1.0 for x in cfbs)


# --- Torching and crowning indices --------------------------------------------

def test_indices_positive_and_ordered():
    fm = pyflam.get_fuel_model(10)
    ti = pyflam.torching_index(fm, canopy_base_height=1.5, foliar_moisture=100,
                               **SCENARIO)
    ci = pyflam.crowning_index(0.15, **SCENARIO)
    assert ti > 0.0 and math.isfinite(ti)
    assert ci > 0.0 and math.isfinite(ci)


def test_higher_canopy_base_raises_torching_index():
    fm = pyflam.get_fuel_model(10)
    low = pyflam.torching_index(fm, canopy_base_height=1.0, foliar_moisture=100,
                                **SCENARIO)
    high = pyflam.torching_index(fm, canopy_base_height=4.0, foliar_moisture=100,
                                 **SCENARIO)
    assert high > low  # a taller canopy base is harder to torch


def test_denser_canopy_lowers_crowning_index():
    sparse = pyflam.crowning_index(0.05, **SCENARIO)
    dense = pyflam.crowning_index(0.25, **SCENARIO)
    assert dense < sparse  # denser canopy crowns at a lower wind


# --- Landscape raster ---------------------------------------------------------

def test_crown_fire_potential_raster():
    n = 12
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        aspect=np.zeros((n, n)),
        canopy_height=np.full((n, n), 150, dtype=int),       # 15 m (m*10)
        canopy_base_height=np.full((n, n), 15, dtype=int),   # 1.5 m (m*10)
        canopy_bulk_density=np.full((n, n), 15, dtype=int),  # 0.15 kg/m^3 (*100)
    )
    out = pyflam.crown_fire_potential(
        ls, foliar_moisture=100,
        wind_20ft_ft_per_min=mph_to_ft_per_min(30),
        wind_midflame=mph_to_ft_per_min(0.4 * 30), **SCENARIO)
    assert out["fire_type"].shape == (n, n)
    assert set(np.unique(out["fire_type"])).issubset({0, 1, 2})
    assert np.all(out["crown_fraction_burned"] >= 0.0)
    assert np.all(out["crown_fraction_burned"] <= 1.0)
    assert np.all(out["rate_of_spread"] >= 0.0)


def test_crown_fire_potential_requires_canopy_bands():
    n = 5
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int), slope=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0)
    with pytest.raises(ValueError):
        pyflam.crown_fire_potential(
            ls, foliar_moisture=100, wind_20ft_ft_per_min=0.0, **SCENARIO)


# --- Golden-master ------------------------------------------------------------
# Provisional values from this implementation; replace with FlamMap/NEXUS-
# verified numbers once cross-checked (see tests/REFERENCE.md).

def test_golden_master_indices():
    fm = pyflam.get_fuel_model(10)
    ti = pyflam.torching_index(fm, canopy_base_height=1.5, foliar_moisture=100,
                               **SCENARIO)
    ci = pyflam.crowning_index(0.15, **SCENARIO)
    assert ti == pytest.approx(6.04, abs=0.05)
    assert ci == pytest.approx(20.66, abs=0.05)
