"""Tests for the vertical-profile pyroconvection diagnostics (LCL, Continuous
Haines, inverted-V detector, pyroconvection_potential).

The physics: high-convective / pyroCb-prone fire weather is a deep dry mixed layer
(high LCL) capped by moisture aloft -- the inverted-V sounding -- NOT high surface
CAPE. The diagnostics must flag that geometry and reject a moist, stable column.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam.atmosphere import (
    AtmosphericProfile, AtmosphericState, briggs_buoyancy_flux, briggs_plume_rise,
    brunt_vaisala_squared, continuous_haines, convective_plume_factor,
    dewpoint_from_rh, inverted_v, lcl_height_m, pyrocb_firepower_threshold,
    pyroconvection_potential,
)


# --- LCL + dewpoint -----------------------------------------------------------

def test_dewpoint_inverts_rh():
    # round-trip: RH -> Td -> RH via the Magnus relation
    from pyflam.atmosphere import relative_humidity_from_dewpoint
    td = dewpoint_from_rh(30.0, 40.0)
    assert relative_humidity_from_dewpoint(30.0, td) == pytest.approx(40.0, abs=0.5)


def test_lcl_higher_when_drier():
    assert lcl_height_m(35, 15) > lcl_height_m(35, 60) > lcl_height_m(20, 90)
    assert lcl_height_m(20, 100) == pytest.approx(0.0, abs=1.0)   # saturated -> ~0


def test_lcl_vectorized():
    h = lcl_height_m(np.array([35.0, 20.0]), np.array([15.0, 90.0]))
    assert h.shape == (2,) and h[0] > h[1]


# --- Continuous Haines --------------------------------------------------------

def _profile(P, T, RH):
    return AtmosphericProfile.from_rh(P, T, RH)


def test_chaines_higher_for_dry_unstable():
    P = [1000, 850, 700, 600, 500]
    dry_unstable = _profile(P, [36, 22, 10, 2, -8], [15, 12, 18, 30, 40])
    moist_stable = _profile(P, [22, 17, 12, 6, -2], [85, 80, 75, 70, 60])
    assert continuous_haines(dry_unstable) > continuous_haines(moist_stable)
    assert continuous_haines(moist_stable) < 4.0


def test_chaines_interpolates_levels():
    # profile given at non-standard levels still resolves 850/700
    p = _profile([1013, 900, 800, 650, 500], [30, 20, 14, 6, -5], [30, 25, 20, 40, 55])
    ch = continuous_haines(p)
    assert np.isfinite(ch)


# --- inverted-V ---------------------------------------------------------------

def test_inverted_v_detects_dry_low_moist_aloft():
    iv = _profile([1000, 850, 700, 600, 500], [36, 22, 12, 4, -6], [15, 12, 20, 75, 60])
    flag, sfc_depr, mid_rh = inverted_v(iv)
    assert flag is True
    assert sfc_depr > 20.0 and mid_rh >= 50.0


def test_inverted_v_rejects_moist_boundary_layer():
    mp = _profile([1000, 850, 700, 600, 500], [22, 16, 10, 4, -4], [80, 75, 70, 30, 25])
    flag, _, _ = inverted_v(mp)
    assert flag is False                        # moist low + dry aloft = not inverted-V


# --- pyroconvection_potential -------------------------------------------------

def _hot_dry_state(pbl=3200):
    return AtmosphericState(wind_speed=4, wind_direction=270, temperature=36,
                            relative_humidity=15, cape=20, boundary_layer_height=pbl)


def test_potential_flags_pyrocb_prone_column():
    iv = _profile([1000, 850, 700, 600, 500], [36, 22, 12, 4, -6], [15, 12, 20, 75, 60])
    pp = pyroconvection_potential(_hot_dry_state(), profile=iv)
    assert pp.plume_dominated_favorable is True
    assert pp.lcl_height_m > 2000.0 and pp.inverted_v is True
    assert pp.continuous_haines is not None


def test_potential_rejects_moist_stable_column():
    mp = _profile([1000, 850, 700, 600, 500], [22, 16, 10, 4, -4], [80, 75, 70, 65, 60])
    st = AtmosphericState(wind_speed=4, wind_direction=270, temperature=22,
                          relative_humidity=80, cape=100, boundary_layer_height=700)
    pp = pyroconvection_potential(st, profile=mp)
    assert pp.plume_dominated_favorable is False


def test_potential_not_driven_by_surface_cape():
    """A high-surface-CAPE but moist/shallow column must NOT flag pyroconvection
    (the literature's key point: surface CAPE is a poor pyroCb predictor)."""
    mp = _profile([1000, 850, 700, 600, 500], [24, 18, 12, 6, -2], [85, 80, 78, 72, 65])
    st = AtmosphericState(wind_speed=3, wind_direction=270, temperature=24,
                          relative_humidity=85, cape=3000, boundary_layer_height=900)
    pp = pyroconvection_potential(st, profile=mp)
    assert pp.plume_dominated_favorable is False     # despite CAPE=3000


def test_potential_surface_only_fallback():
    pp = pyroconvection_potential(_hot_dry_state())   # no profile
    assert pp.continuous_haines is None and pp.inverted_v is None
    assert pp.lcl_height_m > 2000.0
    assert pp.plume_dominated_favorable is True       # high LCL + deep PBL
    assert "surface-only" in pp.notes


def test_potential_rejects_array_state():
    st = AtmosphericState(wind_speed=np.zeros((2, 2)), wind_direction=np.zeros((2, 2)),
                          temperature=np.full((2, 2), 30.0),
                          relative_humidity=np.full((2, 2), 20.0))
    with pytest.raises(ValueError):
        pyroconvection_potential(st)


# --- convective_plume_factor re-weighting -------------------------------------

def _iv_profile():
    return _profile([1000, 850, 700, 600, 500], [36, 22, 12, 4, -6], [15, 12, 20, 75, 60])


def _moist_profile():
    return _profile([1000, 850, 700, 600, 500], [22, 16, 10, 4, -4], [80, 75, 70, 65, 60])


def test_plume_factor_backward_compatible_without_profile():
    """No profile -> unchanged CAPE/stability behaviour."""
    capey = AtmosphericState(wind_speed=3, wind_direction=270, temperature=30,
                             relative_humidity=20, cape=3000.0, sensible_heat_flux=300.0)
    calm = AtmosphericState(wind_speed=3, wind_direction=270, temperature=20,
                            relative_humidity=40)
    assert convective_plume_factor(capey) > convective_plume_factor(calm)
    assert 0.5 <= convective_plume_factor(capey) <= 3.0


def test_plume_factor_profile_boosts_inverted_v_not_moist():
    st = AtmosphericState(wind_speed=3, wind_direction=270, temperature=36,
                          relative_humidity=15, boundary_layer_height=3200)
    base = convective_plume_factor(st)
    boosted = convective_plume_factor(st, profile=_iv_profile())
    moist = convective_plume_factor(st, profile=_moist_profile())
    assert boosted > base                       # inverted-V profile raises the factor
    assert moist == pytest.approx(base, abs=1e-6)   # moist column adds nothing
    assert boosted <= 3.0                        # still bounded


# --- Briggs plume rise + PFT --------------------------------------------------

def test_buoyancy_flux_linear_in_heat():
    assert briggs_buoyancy_flux(2e8) == pytest.approx(2 * briggs_buoyancy_flux(1e8))


def test_briggs_rise_monotonic():
    Q = 5e8
    assert briggs_plume_rise(2 * Q, 5, stability_s2=3e-4) > \
        briggs_plume_rise(Q, 5, stability_s2=3e-4)            # hotter -> higher
    assert briggs_plume_rise(Q, 5, stability_s2=3e-4) > \
        briggs_plume_rise(Q, 12, stability_s2=3e-4)          # more wind -> lower
    assert briggs_plume_rise(Q, 5, stability_s2=1e-4) > \
        briggs_plume_rise(Q, 5, stability_s2=5e-4)           # less stable -> higher


def test_briggs_neutral_needs_distance_and_grows():
    Q = 5e8
    with pytest.raises(ValueError):
        briggs_plume_rise(Q, 5)                              # neutral, no distance
    assert briggs_plume_rise(Q, 5, distance_m=4000) > \
        briggs_plume_rise(Q, 5, distance_m=1000)


def test_brunt_vaisala_positive_for_stable_layer():
    assert brunt_vaisala_squared(_iv_profile()) > 0.0


def test_pft_higher_when_lcl_higher():
    """A drier surface (higher LCL) raises the firepower needed to reach condensation."""
    prof = _iv_profile()
    dry = AtmosphericState(wind_speed=5, wind_direction=270, temperature=36,
                           relative_humidity=12)
    less_dry = AtmosphericState(wind_speed=5, wind_direction=270, temperature=36,
                                relative_humidity=35)
    assert pyrocb_firepower_threshold(dry, prof) > \
        pyrocb_firepower_threshold(less_dry, prof)
    assert pyrocb_firepower_threshold(dry, prof) > 0.0


def test_pft_higher_with_stronger_cap():
    """A more stable capping layer needs a more powerful fire to punch through."""
    st = AtmosphericState(wind_speed=5, wind_direction=270, temperature=36,
                          relative_humidity=15)
    weak = _profile([1000, 850, 700, 600, 500], [36, 22, 10, 0, -12], [15, 12, 20, 70, 60])
    strong = _profile([1000, 850, 700, 600, 500], [36, 22, 14, 12, 8], [15, 12, 20, 70, 60])
    assert brunt_vaisala_squared(strong) > brunt_vaisala_squared(weak)
    assert pyrocb_firepower_threshold(st, strong) > pyrocb_firepower_threshold(st, weak)


def test_pft_degenerate_target_is_inf():
    st = AtmosphericState(wind_speed=5, wind_direction=270, temperature=20,
                          relative_humidity=100)            # LCL ~ 0
    assert pyrocb_firepower_threshold(st, _iv_profile(), target_height_m=0.0) == float("inf")
