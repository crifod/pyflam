"""Tests for the crown / plume-coupling building blocks (plan steps 1-2):
the wind-reconciliation helper and the crown-aware spread field."""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import crownfire
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07,
                m_live_herb=0.60, m_live_woody=0.90)


def _mixed_canopy(n=24):
    """Left half: low crown base (will crown); right half: very tall base (won't)."""
    cbh = np.full((n, n), 10.0)          # *10 -> 1 m
    cbh[:, n // 2:] = 2000.0             # *10 -> 200 m, unreachable -> stays surface
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.full((n, n), 20.0, dtype=float), aspect=np.full((n, n), 180.0),
        canopy_base_height=cbh, canopy_bulk_density=np.full((n, n), 20.0),
        canopy_height=np.full((n, n), 150.0, dtype=float),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def _build(ls, crown_spread="cruz2005"):
    return pyflam.crown_spread_field(
        ls, wind_midflame=mph_to_ft_per_min(12), wind_direction=270.0,
        wind_20ft_ft_per_min=mph_to_ft_per_min(30), foliar_moisture=100.0,
        crown_spread=crown_spread, **SCENARIO)


# --- wind reconciliation ------------------------------------------------------

def test_wind_20ft_to_u10_value_and_array():
    # 20 mph 20-ft wind -> ~37 km/h 10-m wind (U10 > U20ft with height).
    assert pyflam.wind_20ft_to_u10_kmh(mph_to_ft_per_min(20)) == pytest.approx(37.0, abs=0.6)
    arr = pyflam.wind_20ft_to_u10_kmh(np.array([mph_to_ft_per_min(10),
                                                mph_to_ft_per_min(20)]))
    assert arr.shape == (2,) and arr[1] > arr[0]


def test_wind_helper_is_the_single_conversion():
    """active_crown_ros_cruz must use exactly the helper's U10."""
    w = mph_to_ft_per_min(25)
    u10 = pyflam.wind_20ft_to_u10_kmh(w)
    direct = 11.02 * u10 ** 0.90 * 0.18 ** 0.19 * np.exp(-0.17 * 6.0)
    assert crownfire.active_crown_ros_cruz(
        w, canopy_bulk_density=0.18, m_1h=0.06) == pytest.approx(direct, rel=1e-6)


# --- crown-aware spread field -------------------------------------------------

def test_crown_aware_field_mixes_surface_and_crown():
    caf = _build(_mixed_canopy())
    surface = caf.fire_type == 0
    crowning = caf.fire_type >= 1
    assert surface.any() and crowning.any()           # genuinely mixed
    # surface cells are untouched; crowning cells are faster + more intense
    assert np.allclose(caf.field.ros_max[surface], caf.surface.ros_max[surface])
    assert (caf.field.ros_max[crowning] >= caf.surface.ros_max[crowning] - 1e-6).all()
    assert (caf.field.fireline_intensity[crowning]
            >= caf.surface.fireline_intensity[crowning] - 1e-6).all()


def test_crown_aware_keeps_ellipse_geometry():
    caf = _build(_mixed_canopy())
    np.testing.assert_array_equal(caf.field.eccentricity, caf.surface.eccentricity)
    np.testing.assert_array_equal(caf.field.heading, caf.surface.heading)


def test_fire_type_matches_crown_fire_potential():
    ls = _mixed_canopy()
    caf = _build(ls)
    ref = pyflam.crownfire.crown_fire_potential(
        ls, foliar_moisture=100.0, wind_20ft_ft_per_min=mph_to_ft_per_min(30),
        wind_midflame=mph_to_ft_per_min(12), crown_spread="cruz2005", **SCENARIO)
    np.testing.assert_array_equal(caf.fire_type, ref["fire_type"])


def test_crown_aware_field_drives_mtt_and_burn_probability():
    caf = _build(_mixed_canopy())
    arr = pyflam.minimum_travel_time(caf.field, [(12, 12)], max_time=30.0)
    assert np.isfinite(arr[12, 12])
    prob, nf = pyflam.burn_probability(caf.field, [(12, 12)] * 3, max_time=20.0)
    assert nf == 3 and prob[12, 12] == 1.0


def test_cruz_field_faster_than_rothermel_on_crown_cells():
    ls = _mixed_canopy()
    cruz = _build(ls, "cruz2005")
    roth = _build(ls, "rothermel1991")
    active = (cruz.fire_type == 2) & (roth.fire_type == 2)
    if active.any():
        assert (cruz.field.ros_max[active].mean()
                >= roth.field.ros_max[active].mean())
