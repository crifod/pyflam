"""Tests for crown fire in the pyroconvection march (coupling-plan step 3).

`fire_atmosphere_march(crown=True)` rebuilds a crown-aware spread field each
increment from the current (plume-modified) wind, so crown intensity feeds the
plume and the crown rate of spread drives growth. The wind provider is injected
(no OpenFOAM) so the coupling and its feedback are testable in isolation.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import pyroconvection

SC = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07, m_live_herb=0.6, m_live_woody=0.9)


def _canopy_ls(n=31):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        canopy_base_height=np.full((n, n), 10.0),     # *0.1 -> 1 m
        canopy_bulk_density=np.full((n, n), 20.0),    # *0.01 -> 0.20 kg/m^3
        canopy_height=np.full((n, n), 150.0),         # *0.1 -> 15 m
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def _flat_ls(n=31):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def _const_wind(ls, intensity, active, spd, dirn):
    return pyroconvection._uniform_wind_field(ls, spd, dirn)


def _intensity_plume(ls, intensity, active, spd, dirn):
    """A plume that strengthens the wind in proportion to the active-front intensity
    (a stand-in for the CFD plume, so the feedback is testable without OpenFOAM)."""
    wf = pyroconvection._uniform_wind_field(ls, spd, dirn)
    wf.speed = wf.speed.copy()
    inten = np.asarray(intensity, dtype=float)
    a = np.asarray(active, dtype=bool)
    if a.any():
        wf.speed[a] = spd + 0.004 * inten[a]          # boost ~ fireline intensity
    return wf


def _march(ls, *, crown=False, wind=_const_wind, **kw):
    return pyflam.fire_atmosphere_march(
        ls, [(15, 15)], total_time=30, dt=10, speed=8.0, direction=270.0,
        wind_provider=wind, crown=crown, **kw, **SC)


# --- basic crown march --------------------------------------------------------

def test_crown_march_runs_and_returns_fire_type():
    res = _march(_canopy_ls(), crown=True, foliar_moisture=100.0)
    at = res["arrival_time"]
    assert at[15, 15] == 0.0 and np.isfinite(at).sum() > 1
    assert at[15, 16] < at[15, 24]                    # grows with distance
    ft = res["fire_type"]
    assert ft.shape == _canopy_ls().shape
    assert (ft >= 1).any()                            # some crowning occurred


def test_crown_fire_outruns_surface_fire():
    surf = _march(_canopy_ls(), crown=False)
    crown = _march(_canopy_ls(), crown=True, foliar_moisture=100.0)
    assert np.isfinite(crown["arrival_time"]).sum() > \
        np.isfinite(surf["arrival_time"]).sum()


# --- the feedback: crown intensity -> stronger plume -> bigger fire -----------

def test_crown_plume_feedback_grows_the_fire():
    const = _march(_canopy_ls(), crown=True, foliar_moisture=100.0, wind=_const_wind)
    plume = _march(_canopy_ls(), crown=True, foliar_moisture=100.0,
                   wind=_intensity_plume)
    # the intensity-driven plume boosts the wind on the (high-intensity) crown front
    assert np.isfinite(plume["arrival_time"]).sum() >= \
        np.isfinite(const["arrival_time"]).sum()


def test_crown_intensity_exceeds_surface_in_the_plume():
    """The plume sees crown intensity under crown=True, surface intensity without."""
    surf = _march(_canopy_ls(), crown=False, wind=_intensity_plume, return_history=True)
    crown = _march(_canopy_ls(), crown=True, foliar_moisture=100.0,
                   wind=_intensity_plume, return_history=True)
    # crown fields carry much higher fireline intensity -> a stronger plume coupling
    si = np.nanmax([np.nanmax(f.fireline_intensity) for f in surf["fields"]])
    ci = np.nanmax([np.nanmax(f.fireline_intensity) for f in crown["fields"]])
    assert ci > si


# --- validation ---------------------------------------------------------------

def test_crown_requires_foliar_moisture():
    with pytest.raises(ValueError):
        _march(_canopy_ls(), crown=True)              # no foliar_moisture


def test_crown_requires_canopy_bands():
    with pytest.raises(ValueError):
        _march(_flat_ls(), crown=True, foliar_moisture=100.0)   # no CBH/CBD bands


def test_crown_false_is_unchanged():
    """Default (surface) path must be untouched."""
    a = _march(_canopy_ls(), crown=False)
    b = _march(_flat_ls(), crown=False)
    assert "fire_type" not in a
    np.testing.assert_array_equal(
        np.isfinite(a["arrival_time"]), np.isfinite(b["arrival_time"]))
