"""Validation of the pyroconvection-TYPE classifier against the paper's own cases.

Castellnou et al. (2022, JGR-Atmos 127, e2022JD036920) classify pyroconvection
into an ordered ladder of prototypes (their Fig. 7 / Table 1). This test checks
that :func:`pyflam.pyroconvection_type` reproduces the labels the paper assigns to
its 2021-campaign events from the same conditioning variables the paper reports
(mixed-layer dtheta/dz, LCL/ABL ratio, upper-layer gamma-theta). The values used
below are the ones stated in the paper's text (§3.3) and Table 1.

This is the fidelity solve for the "categories matched to the ladder, not validated"
limitation: the classifier is anchored to the paper's labelled observations.
"""

from __future__ import annotations

import numpy as np
import pytest

from pyflam.atmosphere import (
    AtmosphericProfile, bulk_richardson_abl_height, pyroconvection_type,
    theta_gradient, PYROCONVECTION_TYPE_LEVEL,
)


# --- Table 1 / §3.3 labelled cases --------------------------------------------
# (ml_theta_gradient K/m, lcl_abl_ratio, gamma_theta K/m, expected_type)
CASES = {
    # T21: stable ML (dtheta/dz > 1.1e-3) -> non-pyroCu convection plume only.
    "T21": (3.7e-3, 1.2, 5.0e-3, "convection_plume"),
    # SCQ32: slightly-stable/neutral ML, LCL/ABL > 1 -> brief overshooting pyroCu.
    "SCQ32": (0.8e-3, 1.3, 4.8e-3, "overshooting_pyrocu"),
    # M11: unstable ML, LCL/ABL < 1, gamma-theta 4.2e-3 (above cap) -> resilient.
    "M11": (2.6e-5, 0.83, 4.2e-3, "resilient_pyrocu"),
    # SCQ51: unstable ML, weak cap gamma-theta 3.9e-3 -> deep pyroCu / pyroCb.
    "SCQ51": (5.0e-4, 1.0, 3.9e-3, "deep_pyrocu_pyrocb"),
}


@pytest.mark.parametrize("name", list(CASES))
def test_reproduces_paper_labels(name):
    mlg, ratio, gth, expected = CASES[name]
    got = pyroconvection_type(lcl_abl_ratio=ratio, ml_theta_gradient=mlg,
                              gamma_theta=gth)
    assert got == expected, f"{name}: expected {expected}, got {got}"


def test_scq41_pyrocu_but_not_pyrocb():
    """SCQ41: unstable ML but a STRONG cap (gamma-theta 5.1e-3) inhibits deepening
    -- the paper classifies it pyroCu, not pyroCb."""
    t = pyroconvection_type(lcl_abl_ratio=0.9, ml_theta_gradient=5.0e-4,
                            gamma_theta=5.1e-3)
    assert t in ("resilient_pyrocu", "overshooting_pyrocu")
    assert t != "deep_pyrocu_pyrocb"


def test_ladder_is_monotonic_in_activity():
    """The five prototypes are ordered low -> high pyroconvective activity."""
    order = ["surface_plume", "convection_plume", "overshooting_pyrocu",
             "resilient_pyrocu", "deep_pyrocu_pyrocb"]
    assert [PYROCONVECTION_TYPE_LEVEL[t] for t in order] == [0, 1, 2, 3, 4]


# --- fire-power gate (the "potential vs realized" solve) ----------------------

def test_fire_power_gate_downgrades_to_surface():
    """A pyroCb-prone column with a fire too weak (FLI < 1e4 kW/m) is only a
    surface plume; raising the fire power unlocks the atmospheric type."""
    kw = dict(lcl_abl_ratio=0.9, ml_theta_gradient=5.0e-4, gamma_theta=3.9e-3)
    assert pyroconvection_type(**kw, fireline_intensity_kw=3_000) == "surface_plume"
    assert pyroconvection_type(**kw, fireline_intensity_kw=25_000) == "deep_pyrocu_pyrocb"


def test_potential_mode_ignores_gate():
    """With no FLI given, the classifier reports the potential type."""
    assert pyroconvection_type(lcl_abl_ratio=0.9, ml_theta_gradient=5.0e-4,
                               gamma_theta=3.9e-3) == "deep_pyrocu_pyrocb"


# --- bulk-Richardson ABL height (the paper's method, solving the PBL limitation)

def test_bulk_richardson_finds_inversion_height():
    """A well-mixed layer (theta ~ const) capped by an inversion: Rib crosses 0.25
    at the cap. We inject a potential-temperature profile directly by passing a
    constant pressure so theta == temperature_c + 273.15."""
    z = np.array([0, 200, 400, 800, 1200, 1600, 2000, 2600, 3200], float)
    theta_k = np.array([305, 305, 305, 305, 305, 305, 309, 315, 322], float)
    t = theta_k - 273.15
    p = np.full_like(z, 1000.0)            # makes theta_kelvin return theta_k
    u = np.array([1, 1.5, 2, 3, 4, 5, 6, 7, 8], float)   # sheared wind
    v = np.zeros_like(z)
    h = bulk_richardson_abl_height(z, t, u, v, pressure_hpa=p, surface_start_m=400.0)
    assert 1500.0 <= h <= 2100.0          # ABL top at the base of the cap
    # weaker shear (smaller denominator -> larger Rib -> crosses lower) lowers ABL
    h_calm = bulk_richardson_abl_height(z, t, u * 3.0, v, pressure_hpa=p,
                                        surface_start_m=400.0)
    assert h_calm >= h                     # more shear -> deeper diagnosed ABL


def test_theta_gradient_sign():
    """Stable layer -> positive dtheta/dz; near-adiabatic -> ~0."""
    stable = AtmosphericProfile.from_rh([1000, 850, 700, 600, 500],
                                        [36, 26, 18, 12, 6], [20, 20, 30, 40, 50])
    assert theta_gradient(stable, 700, 500) > 0.0
