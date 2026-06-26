"""Tests for pyroconvection wiring in fire_atmosphere_march.

With pyroconvection=True the fireline intensity fed to the plume is scaled by the
convective plume factor (profile-aware), so a dry/unstable inverted-V atmosphere
drives a stronger plume -> stronger wind feedback -> larger fire, while a stable
atmosphere damps it. The wind provider is injected so the coupling is testable
without OpenFOAM.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import pyroconvection
from pyflam.atmosphere import AtmosphericProfile, AtmosphericState

SC = dict(m_1h=0.05, m_10h=0.06, m_100h=0.07, m_live_herb=0.6, m_live_woody=0.9)


def _ls(n=41):
    return pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int),
        slope=np.zeros((n, n)), aspect=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")


def _intensity_plume(gain=0.02):
    """Wind boost proportional to the (plume-factored) active-front intensity."""
    def provider(ls, intensity, active, spd, dirn):
        wf = pyroconvection._uniform_wind_field(ls, spd, dirn)
        wf.speed = wf.speed.copy()
        a = np.asarray(active, dtype=bool)
        if a.any():
            wf.speed[a] = spd + gain * np.asarray(intensity, dtype=float)[a]
        return wf
    return provider


def _iv_profile():
    return AtmosphericProfile.from_rh(
        [1000, 850, 700, 600, 500], [36, 22, 12, 4, -6], [15, 12, 20, 75, 60])


def _iv_state():
    return AtmosphericState(wind_speed=8, wind_direction=270, temperature=36,
                            relative_humidity=15, boundary_layer_height=3200)


def _stable_state():
    return AtmosphericState(wind_speed=8, wind_direction=270, temperature=12,
                            relative_humidity=70, sensible_heat_flux=-40.0)


def _march(state=None, profile=None, pyro=False, gain=0.02, **kw):
    return pyflam.fire_atmosphere_march(
        _ls(), [(20, 20)], total_time=40, dt=10, speed=8.0, direction=270.0,
        wind_provider=_intensity_plume(gain), pyroconvection=pyro, state=state,
        profile=profile, return_history=True, **kw, **SC)


# --- plume factor recorded ----------------------------------------------------

def test_plume_factor_is_one_when_off():
    res = _march(state=_iv_state(), profile=_iv_profile(), pyro=False)
    assert all(pf == 1.0 for pf in res["plume_factor"])


def test_plume_factor_boosts_under_inverted_v():
    res = _march(state=_iv_state(), profile=_iv_profile(), pyro=True)
    assert all(pf > 1.3 for pf in res["plume_factor"])      # inverted-V boost


def test_plume_factor_damps_under_stable():
    res = _march(state=_stable_state(), pyro=True)
    assert all(pf < 1.0 for pf in res["plume_factor"])      # stable damping


# --- effect on the fire -------------------------------------------------------

def test_pyroconvective_atmosphere_grows_a_larger_fire():
    off = _march(state=_iv_state(), profile=_iv_profile(), pyro=False, gain=0.04)
    on = _march(state=_iv_state(), profile=_iv_profile(), pyro=True, gain=0.04)
    assert np.isfinite(on["arrival_time"]).sum() > \
        np.isfinite(off["arrival_time"]).sum()


def test_stable_atmosphere_does_not_grow_it():
    off = _march(state=_stable_state(), pyro=False, gain=0.04)
    on = _march(state=_stable_state(), pyro=True, gain=0.04)
    # damped plume (factor 0.7) must not produce a larger fire than no scaling
    assert np.isfinite(on["arrival_time"]).sum() <= \
        np.isfinite(off["arrival_time"]).sum()


# --- diagnostics in the output ------------------------------------------------

def test_output_carries_pyroconvection_and_pft():
    res = _march(state=_iv_state(), profile=_iv_profile(), pyro=True)
    assert "pyroconvection" in res
    assert res["pyroconvection"].plume_dominated_favorable is True
    assert "pyrocb_firepower_threshold" in res
    assert res["pyrocb_firepower_threshold"] > 0.0


def test_no_pft_without_profile():
    res = _march(state=_iv_state(), pyro=True)            # no profile
    assert "pyroconvection" in res                        # still computes (surface-only)
    assert "pyrocb_firepower_threshold" not in res


def test_spatial_per_cell_pyroconvection():
    """A gridded (half dry / half moist) atmosphere drives a per-cell plume factor:
    the dry half (high LCL) is boosted more than the moist half under one column."""
    from pyflam.atmosphere import AtmosphereProvider

    class _HalfDry(AtmosphereProvider):
        def state_at(self, lat, lon, time=None):
            return AtmosphericState(wind_speed=8, wind_direction=270,
                                    temperature=32, relative_humidity=20)

        def field_on(self, ls, time=None, *, latlon=None):
            T = np.full(ls.shape, 32.0)
            RH = np.full(ls.shape, 85.0)
            RH[:, :ls.shape[1] // 2] = 15.0           # dry left half
            return AtmosphericState(
                wind_speed=np.full(ls.shape, 8.0),
                wind_direction=np.full(ls.shape, 270.0),
                temperature=T, relative_humidity=RH)

    ls = _ls()
    prof = AtmosphericProfile.from_rh([1000, 850, 700, 600, 500],
                                      [32, 22, 12, 4, -6], [20, 12, 20, 75, 60])
    res = pyflam.fire_atmosphere_march(
        ls, [(20, 10)], total_time=40, dt=10, atmosphere=_HalfDry(), spatial=True,
        plume=False, wind_provider=_intensity_plume(0.02), pyroconvection=True,
        profile=prof, return_history=True, m_live_herb=0.6, m_live_woody=0.9)
    assert res["arrival_time"][20, 10] == 0.0
    # the per-cell factor itself: dry cells boosted above moist cells
    n = ls.shape[1]
    fld = _HalfDry().field_on(ls)
    pf = pyflam.convective_plume_factor(fld, profile=prof)
    assert pf.shape == ls.shape
    assert pf[:, :n // 2].mean() > pf[:, n // 2:].mean()


def test_off_is_backward_compatible():
    """pyroconvection=False must reproduce the un-scaled run exactly."""
    a = _march(state=_iv_state(), profile=_iv_profile(), pyro=False)
    b = pyflam.fire_atmosphere_march(
        _ls(), [(20, 20)], total_time=40, dt=10, speed=8.0, direction=270.0,
        wind_provider=_intensity_plume(0.02), return_history=True, **SC)
    np.testing.assert_array_equal(a["arrival_time"], b["arrival_time"])
    assert "pyroconvection" not in b
