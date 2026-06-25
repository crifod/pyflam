"""Tests for crown-driven ember spotting (coupling-plan step 5).

A crown fire's much higher fireline intensity lofts firebrands higher, so they
travel farther. The spotting models already read ``field.fireline_intensity``,
which is the crown intensity on a crown-aware field, so feeding them a crown field
must produce stronger / farther spotting than the surface field.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import pyflam
from pyflam.units import mph_to_ft_per_min

SC = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07, m_live_herb=0.6, m_live_woody=0.9)


def _canopy_ls(n=41):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        canopy_base_height=np.full((n, n), 10.0),     # *0.1 -> 1 m
        canopy_bulk_density=np.full((n, n), 20.0),    # *0.01 -> 0.20 kg/m^3
        canopy_height=np.full((n, n), 150.0),         # *0.1 -> 15 m
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def _crown_aware(ls):
    return pyflam.crown_spread_field(
        ls, wind_midflame=mph_to_ft_per_min(12), wind_direction=270.0,
        wind_20ft_ft_per_min=mph_to_ft_per_min(30), foliar_moisture=100.0,
        crown_spread="cruz2005", **SC)


# --- deterministic loft physics -----------------------------------------------

def test_max_spot_distance_farther_under_crown_intensity():
    caf = _crown_aware(_canopy_ls())
    crowning = caf.fire_type >= 1
    assert crowning.any()
    model = pyflam.SpottingModel()
    w20 = mph_to_ft_per_min(30)
    d_crown = model.max_spot_distance(caf.field.fireline_intensity[crowning], w20)
    d_surf = model.max_spot_distance(caf.surface.fireline_intensity[crowning], w20)
    assert (d_crown >= d_surf - 1e-9).all()
    assert d_crown.max() > d_surf.max()              # crown reaches strictly farther


def test_firebrand_loft_increases_with_intensity():
    fp = pyflam.FirebrandPhysics()
    wt = 2.0                                          # m/s terminal velocity
    assert fp.loft_height(20000.0, wt) > fp.loft_height(2000.0, wt)   # kW/m


# --- end-to-end: crown field lands embers farther -----------------------------

def _burned_patch(shape, center, r=1):
    arrival = np.full(shape, math.inf)
    cr, cc = center
    arrival[cr - r:cr + r + 1, cc - r:cc + r + 1] = 0.0
    return arrival


def _max_reach(spots, center):
    if not spots:
        return 0.0
    cr, cc = center
    return max(math.hypot(r - cr, c - cc) for r, c, _ in spots)


def test_generate_spots_reach_farther_from_crown_field():
    # Big grid + a modest *spotting* wind so both reaches land in-domain (the crown
    # fire's ~40x intensity would otherwise fling embers right off a small grid).
    ls = _canopy_ls(120)
    caf = _crown_aware(ls)
    center = (60, 10)                                 # upwind (west) side; wind FROM 270 -> east
    arrival = _burned_patch(ls.shape, center)
    model = pyflam.SpottingModel(spot_probability=1.0, launch_fraction=1.0,
                                 min_intensity=0.0, spot_delay=0.0)
    kw = dict(wind_20ft=mph_to_ft_per_min(8), wind_direction=270.0, max_time=240.0)

    # average reach over several seeds to beat the lognormal landing noise
    reach_crown, reach_surf = [], []
    for seed in range(8):
        rng = np.random.default_rng(seed)
        reach_crown.append(_max_reach(
            model.generate_spots(caf.field, arrival, rng=rng, **kw), center))
        rng = np.random.default_rng(seed)
        reach_surf.append(_max_reach(
            model.generate_spots(caf.surface, arrival, rng=rng, **kw), center))
    assert np.mean(reach_crown) > np.mean(reach_surf)


def test_surface_field_spotting_unchanged_by_crown_path():
    """The crown-aware field's *surface* attribute is the plain surface field, so
    spotting from it matches spotting from a directly built surface field."""
    ls = _canopy_ls()
    caf = _crown_aware(ls)
    direct = pyflam.spread_field(ls, wind_midflame=mph_to_ft_per_min(12),
                                 wind_direction=270.0, **SC)
    np.testing.assert_allclose(caf.surface.fireline_intensity,
                               direct.fireline_intensity)
