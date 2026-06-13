"""End-to-end RANS tests that actually run OpenFOAM.

Skipped automatically where OpenFOAM is not installed, so the suite stays green
in CI. Where present, these run real (small, fast) solves and check the physics.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import cfd
from pyflam.cfd import run

pytestmark = pytest.mark.skipif(
    run.find_openfoam() is None, reason="OpenFOAM not installed")


def _hill(n=12, relief=80.0, sigma=3.0):
    yy, xx = np.mgrid[0:n, 0:n]
    c = n / 2 - 0.5
    return 1000.0 + relief * np.exp(-(((xx - c) ** 2 + (yy - c) ** 2)
                                      / (2 * sigma ** 2)))


def test_neutral_ridge_speedup():
    n = 12
    wf = cfd.solve_rans(_hill(n), 30.0, speed=6.0, direction=270.0,
                        nz=8, iterations=250, domain_height=400.0)
    c = n // 2
    assert wf.speed[c, c] > wf.speed[c, 0]          # crest faster than upwind
    veer = abs(((wf.direction.mean() - 270.0 + 180.0) % 360.0) - 180.0)
    assert veer < 25.0                              # roughly westerly overall
    assert np.all(wf.speed >= 0)


def test_neutral_flat_preserves_wind():
    """Flat terrain: the ABL wind stays ~uniform and keeps its direction
    (a Richards & Hoxey horizontal-homogeneity sanity check)."""
    flat = np.full((8, 8), 1000.0)
    wf = cfd.solve_rans(flat, 30.0, speed=6.0, direction=270.0,
                        nz=8, iterations=250, domain_height=400.0)
    # interior cells (avoid inlet/outlet edges)
    interior = wf.speed[1:-1, 1:-1]
    assert interior.std() / interior.mean() < 0.15
    veer = abs(((wf.direction[1:-1, 1:-1].mean() - 270.0 + 180.0) % 360.0) - 180.0)
    assert veer < 5.0


def test_buoyant_responds_to_surface_heat_flux():
    """Day (heated) vs night (cooled) surface gives different near-surface wind,
    confirming the buoyancy coupling is active."""
    hill = _hill(10, relief=50.0, sigma=3.0)
    common = dict(cellsize=30.0, speed=5.0, direction=270.0, buoyant=True,
                  nz=8, iterations=200, domain_height=400.0)
    day = cfd.solve_rans(hill, surface_heat_flux=60.0, **common)
    night = cfd.solve_rans(hill, surface_heat_flux=-60.0, **common)
    assert np.all(np.isfinite(day.speed)) and np.all(np.isfinite(night.speed))
    assert np.abs(day.speed - night.speed).max() > 0.05


def test_neutral_feeds_fire_behavior():
    n = 10
    elev = _hill(n).astype(np.int16)
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 1, np.int16),
        slope=np.zeros((n, n), np.int16), elevation=elev,
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0)
    wf = cfd.wind_field_from_landscape(ls, speed=6.0, direction=270.0,
                                       nz=8, iterations=200, domain_height=400.0)
    midflame = wf.to_midflame(ls, wind_reduction_factor=0.4)
    out = pyflam.basic_fire_behavior(
        ls, wind_midflame=midflame, m_1h=0.06, m_10h=0.07, m_100h=0.08,
        m_live_herb=0.60, m_live_woody=0.90)
    assert np.all(np.isfinite(out["rate_of_spread"]))
    assert np.all(out["rate_of_spread"] > 0)
