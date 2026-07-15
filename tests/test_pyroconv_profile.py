# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Cristiano Foderi <cristiano.foderi@gmail.com>

"""The profile-based pyroconvection path: Rib ABL, Bolton LCL, shear, ladders.

Companion to test_pyroconvection_type.py (which validates the default three-diagnostic
ladder against the paper's labelled cases). Here we check the machinery that the
profile ladders add on top: a bulk-Richardson ABL over a grid, an exact LCL, the
three shear-height implementations and the ladder selection rules -- on synthetic
columns whose answers are known by construction, so no GRIB download is needed.
"""

from __future__ import annotations

import numpy as np
import pytest

from pyflam.atmosphere import (
    DEFAULT_PYROCONV_THRESHOLDS, bulk_richardson_abl_grid, lcl_height_bolton_m,
    lcl_height_m, pyroconvection_score, pyroconvection_type,
    pyroconvection_type_castellnou, pyroconvection_type_noshear,
    pyroconvection_type_shear, saturation_vapour_pressure_pa, shear_distance_ratio,
    shear_height_adaptive, shear_height_none, shear_height_window,
    specific_humidity_from_rh, virtual_potential_temperature,
)


# --- moist thermodynamics ------------------------------------------------------

def test_saturation_vapour_pressure_reference_values():
    """Bolton (1980) es: ~611 Pa at 0 C, ~2339 Pa at 20 C, ~4246 Pa at 30 C."""
    es = saturation_vapour_pressure_pa(np.array([273.15, 293.15, 303.15]))
    assert es[0] == pytest.approx(611.2, rel=0.01)
    assert es[1] == pytest.approx(2339.0, rel=0.02)
    assert es[2] == pytest.approx(4246.0, rel=0.02)


def test_specific_humidity_and_virtual_theta():
    """q rises with RH; theta_v exceeds theta for moist air and equals it when dry."""
    q_dry = specific_humidity_from_rh(0.0, 300.0, 100000.0)
    q_moist = specific_humidity_from_rh(80.0, 300.0, 100000.0)
    assert q_dry == pytest.approx(0.0, abs=1e-9)
    assert 0.005 < q_moist < 0.030
    assert virtual_potential_temperature(300.0, q_dry) == pytest.approx(300.0)
    assert virtual_potential_temperature(300.0, q_moist) > 300.0


def test_bolton_lcl_saturated_and_dry():
    """A saturated parcel condenses at the surface; a dry one lifts to a high LCL."""
    assert float(lcl_height_bolton_m(300.0, 300.0, 100000.0)) == pytest.approx(0.0, abs=1.0)
    z = float(lcl_height_bolton_m(308.15, 283.15, 100000.0))   # 35 C / 10 C dewpoint
    assert 2500.0 < z < 3600.0                                  # ~125 m per C of spread


def test_bolton_agrees_with_espy_within_a_few_hundred_metres():
    """The exact LCL and the Espy rule (125 m per C) should not diverge wildly."""
    t_c, rh = 30.0, 35.0
    espy = lcl_height_m(t_c, rh)
    from pyflam.atmosphere import dewpoint_from_rh
    td_c = dewpoint_from_rh(t_c, rh)
    bolton = float(lcl_height_bolton_m(t_c + 273.15, td_c + 273.15, 100000.0))
    assert abs(bolton - espy) < 400.0


# --- bulk-Richardson ABL over a grid -------------------------------------------

def _column_to_grid(values, ny=2, nx=3):
    """Broadcast a 1-D level profile to a (nlev, ny, nx) stack."""
    a = np.asarray(values, float)
    return np.repeat(np.repeat(a[:, None, None], ny, 1), nx, 2)


def test_rib_grid_finds_the_capping_inversion():
    """Neutral layer to 1200 m under a strong inversion -> ABL at the base of the cap.

    theta_v is exactly constant through the mixed layer (Rib stays 0), so the first
    Rib = 0.33 crossing is forced by the inversion itself and must land just above
    1200 m. Note how little stratification it takes: even +0.5 K at 1200 m under a
    light wind pushes Rib well past the threshold, which is why the mixed layer here
    has to be strictly neutral.
    """
    z = np.array([50, 300, 700, 1200, 1800, 2600, 3400], float)
    thv = np.array([305, 305, 305, 305, 311, 317, 323], float)   # cap above 1200 m
    u = np.array([1, 2, 3, 4, 5, 6, 7], float)
    v = np.zeros_like(u)

    abl = bulk_richardson_abl_grid(
        _column_to_grid(z), _column_to_grid(thv), _column_to_grid(u), _column_to_grid(v),
        theta_v_surface=np.full((2, 3), 305.0),
        wind_u_surface=np.zeros((2, 3)), wind_v_surface=np.zeros((2, 3)))
    assert abl.shape == (2, 3)
    assert np.all((abl > 1200.0) & (abl < 1500.0))
    assert np.allclose(abl, abl[0, 0])          # identical columns -> identical ABL


def test_rib_grid_shallower_abl_for_a_stable_profile():
    """A profile stable right from the surface caps the ABL far lower than a mixed one."""
    z = np.array([50, 300, 700, 1200, 1800, 2600, 3400], float)
    mixed = np.array([305, 305, 305, 305.5, 311, 317, 323], float)
    stable = np.array([305, 309, 314, 320, 326, 332, 338], float)
    u = np.array([1, 2, 3, 4, 5, 6, 7], float)
    v = np.zeros_like(u)
    kw = dict(theta_v_surface=np.full((2, 3), 305.0),
              wind_u_surface=np.zeros((2, 3)), wind_v_surface=np.zeros((2, 3)))

    a_mixed = bulk_richardson_abl_grid(*[_column_to_grid(x) for x in (z, mixed, u, v)], **kw)
    a_stable = bulk_richardson_abl_grid(*[_column_to_grid(x) for x in (z, stable, u, v)], **kw)
    assert np.all(a_stable < a_mixed)
    assert np.all(a_stable >= DEFAULT_PYROCONV_THRESHOLDS.min_abl_m)


def test_rib_grid_ignores_underground_levels():
    """A column whose bottom level is below ground must still get a real ABL.

    Over the Apennines the 1000 hPa level is underground, so it arrives as nan. A nan
    level must not advance the Rib bracket -- if it does, it poisons the interpolation
    for every level above it and the whole hill silently comes back unclassified.
    Regression: exactly that bug.
    """
    z = np.array([np.nan, 300, 700, 1200, 1800, 2600, 3400], float)   # 1st underground
    thv = np.array([np.nan, 305, 305, 305, 311, 317, 323], float)
    u = np.array([np.nan, 2, 3, 4, 5, 6, 7], float)
    v = np.zeros_like(u); v[0] = np.nan

    abl = bulk_richardson_abl_grid(
        *[_column_to_grid(x) for x in (z, thv, u, v)],
        theta_v_surface=np.full((2, 3), 305.0),
        wind_u_surface=np.zeros((2, 3)), wind_v_surface=np.zeros((2, 3)))
    assert np.all(np.isfinite(abl)), "underground level poisoned the ABL"
    assert np.all((abl > 1200.0) & (abl < 1500.0))


def test_rib_grid_clips_and_falls_back():
    """A never-crossing column takes the fallback; results stay inside the clip range."""
    z = np.array([50, 300, 700, 1200, 1800], float)
    thv = np.full(5, 305.0)                       # perfectly neutral: Rib never crosses
    u = np.array([1, 2, 3, 4, 5], float)
    v = np.zeros_like(u)
    abl = bulk_richardson_abl_grid(
        *[_column_to_grid(x) for x in (z, thv, u, v)],
        theta_v_surface=np.full((2, 3), 305.0),
        wind_u_surface=np.zeros((2, 3)), wind_v_surface=np.zeros((2, 3)),
        fallback_m=np.full((2, 3), 900.0))
    assert np.allclose(abl, 900.0)
    th = DEFAULT_PYROCONV_THRESHOLDS
    assert np.all((abl >= th.min_abl_m) & (abl <= th.max_abl_m))


# --- the three shear-height implementations ------------------------------------

def test_shear_window_locates_a_planted_shear_layer():
    """A jet-like kink at ~2000 m is found by the moving-window maximum."""
    z = np.linspace(0.0, 6000.0, 61)                  # 61 levels: well resolved
    u = np.where(z < 2000.0, 2.0, 2.0 + 0.02 * (z - 2000.0))   # shear starts at 2 km
    u = np.where(z > 2600.0, 2.0 + 0.02 * 600.0, u)            # ...and stops at 2.6 km
    v = np.zeros_like(z)
    zs = shear_height_window(z, u, v)
    assert 1900.0 <= zs <= 2700.0


def test_shear_none_always_declines():
    assert np.isnan(shear_height_none([1, 2, 3], [1, 2, 3], [1, 2, 3]))


def test_shear_adaptive_refuses_a_coarse_profile_but_uses_a_rich_one():
    """The whole point of variant 3: 5 levels -> nan; 61 levels -> a real height.

    This is what keeps the ICON-2I product off the five-diagnostic ladder.
    """
    z_coarse = np.array([110.0, 770.0, 1457.0, 3012.0, 5574.0])     # ICON-2I open data
    u_coarse = np.array([2.0, 4.0, 7.0, 15.0, 25.0])
    v_coarse = np.zeros_like(z_coarse)
    assert np.isnan(shear_height_adaptive(z_coarse, u_coarse, v_coarse))

    z = np.linspace(0.0, 6000.0, 61)
    u = np.where(z < 2000.0, 2.0, 2.0 + 0.02 * (z - 2000.0))
    assert np.isfinite(shear_height_adaptive(z, u, np.zeros_like(z)))


def test_shear_distance_ratio_is_zero_on_the_abl_top():
    assert shear_distance_ratio(1500.0, 1500.0, 900.0) == pytest.approx(0.0)
    assert shear_distance_ratio(3000.0, 1000.0, 900.0) == pytest.approx(2.0)
    assert np.isnan(shear_distance_ratio(np.nan, 1000.0, 900.0))


# --- ladder selection ----------------------------------------------------------

DEEP = dict(lcl_abl_ratio=0.9, ml_theta_gradient=5.0e-4, gamma_theta=3.9e-3)


def test_default_ladder_is_castellnou():
    """The published three-diagnostic ladder stays the default: no silent change."""
    assert pyroconvection_type(**DEEP) == pyroconvection_type_castellnou(**DEEP)
    assert pyroconvection_type(**DEEP) == "deep_pyrocu_pyrocb"


def test_profile_ladders_are_stricter_about_deep_pyrocb():
    """The same column that Castellnou calls deep is NOT deep for the profile ladders
    when the ABL top is dry -- the moisture requirement is the difference."""
    assert pyroconvection_type_castellnou(**DEEP) == "deep_pyrocu_pyrocb"
    dry = pyroconvection_type_noshear(**DEEP, rh_top_abl=40.0)
    assert dry != "deep_pyrocu_pyrocb"
    moist = pyroconvection_type_noshear(**DEEP, rh_top_abl=90.0)
    assert moist == "deep_pyrocu_pyrocb"


def test_shear_ladder_needs_the_shear_layer_close():
    """With everything else favourable, a distant shear maximum blocks class 4."""
    kw = dict(**DEEP, rh_top_abl=90.0)
    assert pyroconvection_type_shear(**kw, shear_distance=0.1) == "deep_pyrocu_pyrocb"
    far = pyroconvection_type_shear(**kw, shear_distance=2.0)
    assert far == "resilient_pyrocu"


def test_adaptive_picks_the_ladder_the_data_support():
    """Five diagnostics -> shear ladder; four -> noshear; three -> castellnou."""
    a = pyroconvection_type(**DEEP, ladder="adaptive", rh_top_abl=90.0, shear_distance=2.0)
    assert a == pyroconvection_type_shear(**DEEP, rh_top_abl=90.0, shear_distance=2.0)

    b = pyroconvection_type(**DEEP, ladder="adaptive", rh_top_abl=90.0)
    assert b == pyroconvection_type_noshear(**DEEP, rh_top_abl=90.0)

    c = pyroconvection_type(**DEEP, ladder="adaptive")
    assert c == pyroconvection_type_castellnou(**DEEP)


def test_explicit_ladder_refuses_missing_diagnostics():
    """Asking for a ladder you cannot feed is an error, not a silent downgrade."""
    with pytest.raises(ValueError, match="rh_top_abl"):
        pyroconvection_type(**DEEP, ladder="noshear")
    with pytest.raises(ValueError, match="shear_distance"):
        pyroconvection_type(**DEEP, ladder="shear", rh_top_abl=90.0)
    with pytest.raises(ValueError, match="unknown ladder"):
        pyroconvection_type(**DEEP, ladder="nonsense")


def test_fire_power_gate_applies_to_every_ladder():
    """The gate is what makes the map 'expected' rather than 'potential'."""
    weak = dict(fireline_intensity_kw=3_000)
    assert pyroconvection_type(**DEEP, **weak) == "surface_plume"
    assert pyroconvection_type_noshear(**DEEP, rh_top_abl=90.0, **weak) == "surface_plume"
    assert pyroconvection_type_shear(**DEEP, rh_top_abl=90.0, shear_distance=0.1,
                                     **weak) == "surface_plume"
    strong = dict(fireline_intensity_kw=25_000)
    assert pyroconvection_type_shear(**DEEP, rh_top_abl=90.0, shear_distance=0.1,
                                     **strong) == "deep_pyrocu_pyrocb"


# --- parcel mixing depth vs Rib ABL (the entrainment-jump problem) -------------

def test_parcel_depth_finds_the_base_of_the_inversion():
    """The parcel method returns the top of the *well-mixed* layer, at the inversion base.

    Note it is NOT universally below the Rib height: which of the two is higher depends on
    the shear (in ICON columns Rib runs 300-600 m deeper; in degraded soundings it runs
    shallower). What is invariant is what the parcel top *means* -- the last height a
    surface parcel mixes to, before the entrainment jump.
    """
    from pyflam.atmosphere import parcel_mixing_depth_grid

    # Well-mixed to 1200 m, then a 6 K inversion (the jump), then free troposphere.
    z = np.array([50, 300, 700, 1200, 1800, 2600, 3400], float)
    th = np.array([305, 305, 305, 305, 311, 314, 317], float)
    sfc = np.full((2, 3), 305.0)

    parcel = parcel_mixing_depth_grid(_column_to_grid(z), _column_to_grid(th), sfc)
    assert np.all((parcel >= 1200.0) & (parcel < 1350.0))   # base of the inversion
    assert np.allclose(parcel, parcel[0, 0])


def test_gradient_to_the_rib_top_absorbs_the_entrainment_jump():
    """Why 'surface_to_abl' was superseded: theta at the Rib top is post-jump.

    On a genuinely well-mixed column, measuring the gradient up to a height inside the
    inversion reports the *inversion*, reading as spuriously stable and gating the column
    out of the pyroCu branch. Same column, two upper bounds, opposite verdicts.
    """
    z = np.array([50, 300, 700, 1200, 1800, 2600, 3400], float)
    th = np.array([305, 305, 305, 305, 311, 314, 317], float)
    thr = DEFAULT_PYROCONV_THRESHOLDS.ml_stable

    mixed_top = 1250.0                      # top of the well-mixed layer
    above_jump = 1800.0                     # inside the entrainment inversion
    g_mixed = (np.interp(mixed_top, z, th) - 305.0) / mixed_top
    g_jump = (np.interp(above_jump, z, th) - 305.0) / above_jump

    assert g_mixed <= thr, "the mixed layer must read as pyroCu-capable"
    assert g_jump > thr, "measuring through the jump must read as spuriously stable"


def test_surface_to_parcel_is_a_mixing_depth_proxy_by_construction():
    """The default ML diagnostic is identically ``excess / depth`` -- document the fact.

    The parcel top is *defined* as theta_sfc + 0.5 K, so the bulk gradient to it carries no
    information about theta beyond the depth itself. It is a mixing-depth proxy for ML
    stability, not a measurement of dtheta/dz, and the code must not pretend otherwise:
    at the 1.1e-3 threshold it says exactly "parcel mixing depth >= 455 m".
    """
    from pyflam.atmosphere import _PARCEL_EXCESS_K
    thr = DEFAULT_PYROCONV_THRESHOLDS.ml_stable
    for depth in (300.0, 455.0, 800.0, 2000.0):
        grad = _PARCEL_EXCESS_K / depth
        assert (grad <= thr) == (depth >= _PARCEL_EXCESS_K / thr)
    assert _PARCEL_EXCESS_K / thr == pytest.approx(455.0, abs=1.0)


def test_ml_stable_threshold_matches_the_paper():
    """Castellnou et al. (2022) sec.2.4.1: stable > 1.1e-3 K/m (not 1.0e-3)."""
    assert DEFAULT_PYROCONV_THRESHOLDS.ml_stable == pytest.approx(1.1e-3)


# --- mixed-layer gradient method ----------------------------------------------

def test_surface_to_abl_gradient_sees_the_mixing_the_window_misses():
    """On a coarse profile, the bulk surface->ABL gradient reports a well-mixed layer
    as neutral/unstable; the reference 0.2-0.8*ABL window does not.

    With ~900 m level spacing the window's bounds straddle the entrainment zone, so it
    returns the capping inversion's gradient (strongly stable) for a column that is in
    fact well mixed -- which then gates the cell out of the pyroCu-capable branch. This
    pins the reason the pipeline defaults to ``ml_method="surface_to_abl"``.
    """
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pyflam_gui.core.pyroconv import profile_diagnostics

    levels = (1000, 925, 850, 700, 500)
    ny, nx = 2, 2
    # Hot surface (307 K) over a well-mixed layer: theta is ~constant to ~1500 m, then
    # a firm inversion above. ICON-2I's real level spacing.
    z = np.array([110.0, 770.0, 1457.0, 3012.0, 5574.0])
    T = np.array([304.0, 298.2, 292.6, 283.0, 262.0])
    RH = np.array([45.0, 55.0, 85.0, 45.0, 35.0])
    U = np.array([2.0, 4.0, 6.0, 12.0, 20.0])
    V = np.zeros(5)
    g = lambda a: np.repeat(np.repeat(a[:, None, None], ny, 1), nx, 2)[:, None, ...]
    d = dict(levels=levels, z_asl=g(z), T=g(T), RH=g(RH), U=g(U), V=g(V),
             T2m=np.full((1, ny, nx), 307.0), Td2m=np.full((1, ny, nx), 290.0),
             U10=np.full((1, ny, nx), 3.0), V10=np.zeros((1, ny, nx)),
             PS=np.full((1, ny, nx), 100500.0),
             hsurf=np.zeros((ny, nx)), frland=np.ones((ny, nx)))

    bulk = profile_diagnostics(d, 0, ml_method="surface_to_abl")["ml_grad"]
    window = profile_diagnostics(d, 0, ml_method="mid_layer")["ml_grad"]
    assert np.nanmedian(bulk) < np.nanmedian(window), (
        "the surface->ABL gradient must not read more stable than the mid-layer window "
        "on a well-mixed column")
    # The bulk gradient sees a pyroCu-capable (neutral/unstable) mixed layer.
    assert np.nanmedian(bulk) <= DEFAULT_PYROCONV_THRESHOLDS.ml_stable


def test_unknown_ml_method_raises():
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pyflam_gui.core.pyroconv import profile_diagnostics
    with pytest.raises(ValueError, match="unknown ml_method"):
        profile_diagnostics({}, 0, ml_method="nonsense")


# --- ICON-EU model-level path (the hybrid atmosphere) -------------------------

def _iconeu_core():
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pyflam_gui.core import pyroconv as pc
    return pc


def test_masked_lin_slope_recovers_a_known_gradient():
    """The in-ML least-squares fit must return the true slope, and nan below 3 points."""
    pc = _iconeu_core()
    z = np.linspace(0, 2000, 11)[:, None, None] * np.ones((11, 2, 3))
    slope_true = 3.0e-3
    theta = 300.0 + slope_true * z
    mask = z <= 1500.0
    got = pc._masked_lin_slope(z, theta, mask)
    assert np.allclose(got, slope_true, atol=1e-6)
    # a column with only two masked levels -> nan (cannot fit)
    thin = z <= 200.0
    assert np.all(np.isnan(pc._masked_lin_slope(z, theta, thin)))


def test_regrid_to_is_identity_on_the_same_grid_and_interpolates():
    pc = _iconeu_core()
    lat = np.array([42.0, 43.0, 44.0]); lon = np.array([10.0, 11.0, 12.0])
    f = np.arange(9.0).reshape(3, 3)
    same = pc.regrid_to(f, lat, lon, lat, lon)
    assert np.allclose(same, f)
    # midpoint of a linear ramp interpolates to the mean of its neighbours
    mid = pc.regrid_to(f, lat, lon, np.array([42.5]), np.array([10.5]))
    assert mid[0, 0] == pytest.approx((f[0, 0] + f[0, 1] + f[1, 0] + f[1, 1]) / 4)


def test_iconeu_diagnostics_measures_a_well_mixed_column():
    """On a synthetic model-level column, fit_in_ml recovers a near-zero ML gradient and
    the diagnostics dict is classifiable -- i.e. the hybrid atmosphere path is sound."""
    pc = _iconeu_core()
    ny, nx = 2, 3
    # 24 model levels, fine near the surface. Well-mixed theta to ~1400 m, then a cap.
    z1d = np.array([10, 40, 90, 160, 250, 350, 460, 590, 720, 870, 1030, 1200, 1390,
                    1580, 1790, 2000, 2230, 2460, 2710, 2960, 3230, 3500, 3790, 4080], float)
    nlev = z1d.size
    theta1d = np.where(z1d <= 1400, 305.0, 305.0 + 6.0e-3 * (z1d - 1400))   # ~0 in ML
    # pressure from a rough hydrostatic column; RH moderate; light sheared wind
    p1d = 1000e2 * np.exp(-z1d / 8500.0)
    T1d = theta1d * (p1d / 1e5) ** 0.286
    QV1d = np.full(nlev, 0.008)
    U1d = 2.0 + 0.004 * z1d; V1d = np.zeros(nlev)
    g = lambda a: np.repeat(np.repeat(a[:, None, None], ny, 1), nx, 2)

    d = dict(lat=np.array([43.0, 43.1]), lon=np.array([11.0, 11.1, 11.2]),
             levels=tuple(range(74, 74 - nlev, -1)),
             z=g(z1d), T=g(T1d), QV=g(QV1d), P=g(p1d), U=g(U1d), V=g(V1d),
             T2m=np.full((ny, nx), 306.0), Td2m=np.full((ny, nx), 291.0),
             PS=np.full((ny, nx), 1000e2), U10=np.full((ny, nx), 2.0),
             V10=np.zeros((ny, nx)), frland=np.ones((ny, nx)), orog=np.zeros((ny, nx)))

    diag = pc.iconeu_diagnostics(d, ml_method="fit_in_ml")
    assert diag["n_levels"] >= 8, "should see many levels inside the mixed layer"
    assert np.all(np.isfinite(diag["ml_grad"])), "gradient measurable on model levels"
    assert np.nanmedian(diag["ml_grad"]) < 1.1e-3, "a well-mixed column reads as capable"
    # the dict must be classifiable and regriddable exactly like the ICON-2I one
    cls, used = pc.classify_profile(diag, ladder="noshear")
    assert cls.shape == (ny, nx)
    rg = pc.regrid_diagnostics(diag, d["lat"], d["lon"],
                               np.array([43.05]), np.array([11.05]))
    assert set(rg) >= {"abl", "ml_grad", "lcl_ratio", "gamma", "rh_top", "valid"}


def test_fetch_icon_eu_builds_expected_urls(monkeypatch, tmp_path):
    """No network: capture the URLs and local names fetch_icon_eu would request."""
    import pyflam.atmosphere as atm
    from datetime import datetime

    calls = []

    def fake_urlretrieve(url, *a, **k):
        calls.append(url)
        p = tmp_path / "raw.bz2"
        import bz2
        p.write_bytes(bz2.compress(b"x"))
        return str(p), None

    monkeypatch.setattr(atm.urllib.request if hasattr(atm, "urllib") else __import__(
        "urllib.request", fromlist=["request"]), "urlretrieve", fake_urlretrieve, raising=False)
    import urllib.request as ur
    monkeypatch.setattr(ur, "urlretrieve", fake_urlretrieve)

    out = atm.fetch_icon_eu(datetime(2026, 7, 14), run=0, step=12,
                            cache_dir=str(tmp_path), levels=(74, 73))
    # model levels for 5 vars x 2 levels, HHL for 74,73,75, 5 surface, FR_LAND
    assert f"{atm.ICON_EU_BASE}/00/t/" in "".join(c for c in calls if "_74_T." in c)
    assert any("model-level_2026071400_012_74_T.grib2.bz2" in c for c in calls)
    assert any("time-invariant_2026071400_75_HHL.grib2.bz2" in c for c in calls)
    assert any(c.endswith("FR_LAND.grib2.bz2") for c in calls)
    assert out["T74"].endswith("T74.grib2") and out["HHL75"].endswith("HHL75.grib2")


# --- the continuous score ------------------------------------------------------

def test_score_is_bounded_and_ranks_columns():
    """A favourable column outscores a hostile one; both stay on 0-100."""
    good = pyroconvection_score(lcl_abl_ratio=1.0, ml_theta_gradient=0.0,
                                gamma_theta=2.5e-3, rh_top_abl=95.0, shear_distance=0.0)
    bad = pyroconvection_score(lcl_abl_ratio=3.0, ml_theta_gradient=5.0e-3,
                               gamma_theta=6.0e-3, rh_top_abl=10.0, shear_distance=3.0)
    assert 0.0 <= bad < good <= 100.0
    assert good > 90.0


def test_score_stays_on_scale_without_shear():
    """Dropping the shear term redistributes its weight rather than shrinking the range."""
    kw = dict(lcl_abl_ratio=1.0, ml_theta_gradient=0.0, gamma_theta=2.5e-3,
              rh_top_abl=95.0)
    assert pyroconvection_score(**kw) > 90.0
    assert pyroconvection_score(**kw, shear_distance=np.nan) > 90.0
