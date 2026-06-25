"""Tests for the plume-feedback stabilizers (coupling-plan step 4):
under-relaxation of the wind and a physical cap, so the
crowning -> plume -> wind -> crown loop cannot run away.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import pyroconvection

SC = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07, m_live_herb=0.6, m_live_woody=0.9)
AMBIENT = 4.0          # ambient wind speed (m/s)


def _ls(n=41):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def _runaway_plume(gain):
    """A plume that amplifies the wind hugely with the active-front intensity."""
    def provider(ls, intensity, active, spd, dirn):
        wf = pyroconvection._uniform_wind_field(ls, spd, dirn)
        wf.speed = wf.speed.copy()
        a = np.asarray(active, dtype=bool)
        if a.any():
            wf.speed[a] = spd + gain * np.asarray(intensity, dtype=float)[a]
        return wf
    return provider


def _march(**kw):
    return pyflam.fire_atmosphere_march(
        _ls(), [(20, 20)], total_time=60, dt=10, speed=AMBIENT, direction=270.0,
        return_history=True, **kw, **SC)


def _max_speed(res):
    return max(float(np.nanmax(w.speed)) for w in res["winds"])


# --- the cap ------------------------------------------------------------------

def test_max_wind_factor_caps_speed():
    res = _march(wind_provider=_runaway_plume(0.5), max_wind_factor=2.5)
    assert _max_speed(res) <= 2.5 * AMBIENT + 1e-6        # never exceeds the cap


def test_without_cap_the_plume_blows_up():
    """Sanity: the runaway plume really does drive the wind far past the cap."""
    res = _march(wind_provider=_runaway_plume(0.5))       # defaults: no cap
    assert _max_speed(res) > 5.0 * AMBIENT


# --- under-relaxation ---------------------------------------------------------

def test_under_relaxation_damps_the_wind():
    fast = _march(wind_provider=_runaway_plume(0.05), wind_relax=1.0)
    damp = _march(wind_provider=_runaway_plume(0.05), wind_relax=0.3)
    assert _max_speed(damp) < _max_speed(fast)            # damping lowers the peak


# --- boundedness --------------------------------------------------------------

def test_guards_bound_a_runaway_fire():
    free = _march(wind_provider=_runaway_plume(0.2))                      # no guards
    held = _march(wind_provider=_runaway_plume(0.2),
                  wind_relax=0.4, max_wind_factor=2.5)
    assert (np.isfinite(held["arrival_time"]).sum()
            <= np.isfinite(free["arrival_time"]).sum())


def test_mean_wind_history_present_and_sized():
    res = _march(wind_provider=_runaway_plume(0.05), max_wind_factor=3.0)
    assert len(res["mean_wind"]) == len(res["winds"]) == len(res["times"])


# --- defaults are a no-op -----------------------------------------------------

def test_defaults_unchanged():
    """wind_relax=1 + max_wind_factor=inf must reproduce the un-stabilized run."""
    a = _march(wind_provider=_runaway_plume(0.03))
    b = _march(wind_provider=_runaway_plume(0.03), wind_relax=1.0,
               max_wind_factor=float("inf"))
    np.testing.assert_array_equal(a["arrival_time"], b["arrival_time"])
