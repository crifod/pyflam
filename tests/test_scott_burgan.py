"""Tests for the Scott & Burgan (2005) 40 fuel models.

Covers the two behaviors that distinguish this set from the original 13:
dynamic herbaceous curing and nonburnable models. As with the 13, the
golden-master numbers below are provisional regression locks, not yet verified
against FlamMap/BehavePlus (see tests/REFERENCE.md).
"""

from __future__ import annotations

import math

import pytest

import pyflam
from pyflam.fuel_models import STANDARD_40
from pyflam.units import mph_to_ft_per_min

SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)

BURNABLE = sorted(n for n, fm in STANDARD_40.items() if fm.is_burnable)
NONBURNABLE = sorted(n for n, fm in STANDARD_40.items() if not fm.is_burnable)


def _fb(key, wind_mph=0.0, slope=0.0, **moist):
    fm = pyflam.get_fuel_model(key)
    return pyflam.spread(
        fm, wind_midflame=mph_to_ft_per_min(wind_mph), slope=slope,
        **{**SCENARIO, **moist},
    )


# --- Set membership -----------------------------------------------------------

def test_set_sizes():
    # 40 burnable (GR9 + GS4 + SH9 + TU5 + TL9 + SB4) plus 5 nonburnable (NB).
    assert len(BURNABLE) == 40
    assert len(NONBURNABLE) == 5
    assert len(pyflam.ALL_FUEL_MODELS) == 13 + 45


def test_lookup_by_code_and_number_agree():
    assert pyflam.get_fuel_model("GR1") is pyflam.get_fuel_model(101)
    assert pyflam.get_fuel_model("sb4") is pyflam.get_fuel_model(204)


def test_unknown_model_raises():
    with pytest.raises(KeyError):
        pyflam.get_fuel_model("ZZ9")
    with pytest.raises(KeyError):
        pyflam.get_fuel_model(500)


# --- Physics / property tests -------------------------------------------------

@pytest.mark.parametrize("n", BURNABLE)
def test_burnable_models_spread(n):
    fb = _fb(n, wind_mph=5)
    assert math.isfinite(fb.rate_of_spread)
    assert fb.rate_of_spread > 0.0


@pytest.mark.parametrize("n", NONBURNABLE)
def test_nonburnable_models_do_not_spread(n):
    fb = _fb(n, wind_mph=20, slope=0.5)
    assert fb.rate_of_spread == 0.0
    assert fb.flame_length == 0.0
    assert fb.reaction_intensity == 0.0


@pytest.mark.parametrize("n", BURNABLE)
def test_wind_and_slope_increase_spread(n):
    base = _fb(n).rate_of_spread
    assert _fb(n, wind_mph=10).rate_of_spread > base
    assert _fb(n, slope=0.5).rate_of_spread > base


# --- Dynamic curing -----------------------------------------------------------

DYNAMIC = sorted(n for n in BURNABLE if STANDARD_40[n].dynamic)
STATIC = sorted(n for n in BURNABLE if not STANDARD_40[n].dynamic)


def test_expected_dynamic_models():
    # All grass and grass-shrub models are dynamic; SH1/SH9/TU1/TU3 too.
    assert DYNAMIC == [101, 102, 103, 104, 105, 106, 107, 108, 109,
                       121, 122, 123, 124, 141, 149, 161, 163]


@pytest.mark.parametrize("n", DYNAMIC)
def test_curing_speeds_up_dynamic_models(n):
    """Drier (more cured) live herb transfers more load to dead -> faster."""
    cured = _fb(n, wind_mph=5, m_live_herb=0.30).rate_of_spread
    green = _fb(n, wind_mph=5, m_live_herb=1.20).rate_of_spread
    assert cured > green


@pytest.mark.parametrize("n", STATIC)
def test_static_models_ignore_live_herb_curing(n):
    """Static models have no herb transfer; live herb moisture changes nothing
    via the curing path (these models carry no live herb load anyway)."""
    a = _fb(n, wind_mph=5, m_live_herb=0.30).rate_of_spread
    b = _fb(n, wind_mph=5, m_live_herb=1.20).rate_of_spread
    assert a == pytest.approx(b)


def test_curing_fraction_endpoints():
    from pyflam.rothermel import _cured_fraction
    assert _cured_fraction(0.10) == 1.0   # very dry -> fully cured
    assert _cured_fraction(0.30) == 1.0   # at dry bound
    assert _cured_fraction(1.20) == 0.0   # at wet bound -> fully green
    assert _cured_fraction(2.00) == 0.0   # very wet -> fully green
    assert _cured_fraction(0.75) == pytest.approx(0.5)  # midpoint


# --- Golden-master regression -------------------------------------------------
# Provisional values from this implementation at the SCENARIO above with
# midflame wind = 5 mph, slope = 0. Cross-check against FlamMap/BehavePlus.
GOLDEN_ROS_FT_MIN = {
    "GR1": 18.32, "GR2": 38.25, "GR4": 77.72, "GS2": 25.18,
    "SH5": 63.01, "SH7": 41.20, "TU1": 3.34, "TU5": 10.19,
    "TL3": 1.92, "TL8": 6.73, "SB3": 30.45,
}


@pytest.mark.parametrize("code,expected", sorted(GOLDEN_ROS_FT_MIN.items()))
def test_golden_master_ros(code, expected):
    assert _fb(code, wind_mph=5).rate_of_spread == pytest.approx(expected, abs=0.01)
