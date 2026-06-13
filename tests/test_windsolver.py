"""Tests for the native mass-consistent wind solver.

These validate the *physics properties* the solver must obey — mass conservation
(divergence-free output), recovery of the input wind over flat ground, ridge
speed-up, and direction recovery — rather than cell-exact agreement with
WindNinja (a terrain-following FEM code). Matching WindNinja's grids is a
separate validation step against that reference tool.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import windsolver as ws
from pyflam.landscape import Landscape


def _gaussian_hill(n=21, relief=100.0, sigma=5.0, base=1000.0):
    yy, xx = np.mgrid[0:n, 0:n]
    c = n // 2
    return base + relief * np.exp(-(((xx - c) ** 2 + (yy - c) ** 2)
                                    / (2.0 * sigma ** 2)))


# --- Flat ground: must reproduce the input wind exactly -----------------------

def test_flat_recovers_input_wind():
    flat = np.full((16, 16), 1500.0)
    wf, diag = ws.solve_mass_consistent(
        flat, 30.0, speed=5.0, direction=225.0, z0=0.1,
        reference_height=6.1, output_height=6.1, return_diagnostics=True)
    assert wf.speed == pytest.approx(5.0, abs=1e-6)
    assert wf.direction == pytest.approx(225.0, abs=1e-6)
    assert diag["max_divergence"] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("direction", [0.0, 90.0, 180.0, 270.0])
def test_flat_direction_convention(direction):
    flat = np.full((12, 12), 1000.0)
    wf = ws.solve_mass_consistent(flat, 30.0, speed=4.0, direction=direction,
                                  z0=0.05)
    assert wf.direction == pytest.approx(direction, abs=1e-6)


# --- Mass conservation --------------------------------------------------------

def test_output_is_divergence_free():
    """The corrected field must conserve mass over real terrain."""
    hill = _gaussian_hill(20, relief=80.0, sigma=5.0)
    _, diag = ws.solve_mass_consistent(
        hill, 30.0, speed=6.0, direction=270.0, z0=0.1,
        top_margin=120.0, return_diagnostics=True)
    # Tiny relative to fluxes ~ speed * cell-face area (order 1e3+).
    assert diag["max_divergence"] < 1e-2


def test_divergence_free_with_rough_terrain():
    rng = np.random.default_rng(1)
    terrain = 1000.0 + np.cumsum(rng.normal(0, 5, (18, 18)), axis=1)
    _, diag = ws.solve_mass_consistent(
        terrain, 30.0, speed=5.0, direction=200.0, z0=0.2,
        top_margin=120.0, return_diagnostics=True)
    assert diag["max_divergence"] < 1e-2


# --- Terrain effects ----------------------------------------------------------

def test_ridge_speedup():
    """Wind accelerates over a hill crest relative to the upwind base."""
    hill = _gaussian_hill(25, relief=120.0, sigma=6.0)
    wf = ws.solve_mass_consistent(hill, 30.0, speed=5.0, direction=270.0,
                                  z0=0.1, top_margin=150.0)
    c = 25 // 2
    crest = wf.speed[c, c]
    upwind = wf.speed[c, 4]          # toward the west edge (wind from west)
    assert crest > upwind
    assert wf.direction.std() > 0.5  # flow direction veers around the hill


def test_more_stable_pushes_flow_over_not_around():
    hill = _gaussian_hill(25, relief=120.0, sigma=6.0)
    c = 25 // 2
    neutral = ws.solve_mass_consistent(hill, 30.0, speed=5.0, direction=270.0,
                                       z0=0.1, top_margin=150.0,
                                       stability_ratio=1.0)
    stable = ws.solve_mass_consistent(hill, 30.0, speed=5.0, direction=270.0,
                                      z0=0.1, top_margin=150.0,
                                      stability_ratio=0.25)
    assert stable.speed[c, c] >= neutral.speed[c, c]


# --- Roughness ----------------------------------------------------------------

def test_roughness_from_fuel():
    fuel = np.array([[1, 5, 10], [101, 142, 188], [91, 999, 2]])
    z0 = ws.roughness_from_fuel(fuel)
    assert z0[0, 0] == pytest.approx(0.03)   # FM1 grass
    assert z0[0, 1] == pytest.approx(0.30)   # FM5 brush
    assert z0[0, 2] == pytest.approx(1.0)    # FM10 timber
    assert z0[1, 0] == pytest.approx(0.03)   # GR1
    assert z0[1, 2] == pytest.approx(1.0)    # TL8
    assert z0[2, 1] == pytest.approx(0.1)    # unknown -> default


# --- Landscape integration ----------------------------------------------------

def test_wind_field_from_landscape_feeds_fire_behavior():
    hill = _gaussian_hill(20, relief=90.0, sigma=5.0).astype(np.int16)
    fuel = np.full((20, 20), 1, dtype=np.int16)     # FM1 grass everywhere
    slope = np.zeros((20, 20), dtype=np.int16)
    ls = Landscape(fuel_model=fuel, slope=slope, elevation=hill,
                   cellsize_x=30.0, cellsize_y=30.0, west=500000.0, north=5000000.0)

    wf = ws.wind_field_from_landscape(ls, speed=5.0, direction=270.0,
                                      top_margin=120.0)
    assert wf.shape == ls.shape
    assert wf.west == 500000.0 and wf.north == 5000000.0

    midflame = wf.to_midflame(ls, wind_reduction_factor=0.4)
    assert midflame.shape == ls.shape
    out = pyflam.basic_fire_behavior(
        ls, wind_midflame=midflame,
        m_1h=0.06, m_10h=0.07, m_100h=0.08,
        m_live_herb=0.60, m_live_woody=0.90)
    assert np.all(np.isfinite(out["rate_of_spread"]))
    assert np.all(out["rate_of_spread"] > 0)


def test_landscape_without_elevation_raises():
    ls = Landscape(fuel_model=np.ones((4, 4), np.int16),
                   slope=np.zeros((4, 4), np.int16),
                   cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=0.0)
    with pytest.raises(ValueError, match="elevation"):
        ws.wind_field_from_landscape(ls, speed=5.0, direction=270.0)
