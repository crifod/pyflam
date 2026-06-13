"""Tests for the Rothermel surface fire spread model.

Two layers of validation:

1. *Physics / property tests* — assert relationships the model must obey (ROS
   rises with wind and slope, falls with moisture, hits zero at the moisture of
   extinction, etc.). These genuinely catch implementation bugs and don't depend
   on any external tool.

2. *Golden-master regression* — locks the current numeric outputs for the 13
   standard fuel models so future refactors can't silently change behavior.
   These are NOT independently-verified reference values yet: the real
   acceptance test is a cell-by-cell diff against FlamMap "Basic Fire Behavior"
   / BehavePlus on the same inputs. See ``tests/REFERENCE.md`` for that harness.
"""

from __future__ import annotations

import math

import pytest

import pyflam
from pyflam.units import mph_to_ft_per_min

# A common reference scenario used throughout: dry dead fuels, green live fuels.
SCENARIO = dict(m_1h=0.06, m_10h=0.07, m_100h=0.08,
                m_live_herb=0.60, m_live_woody=0.90)


def _ros(n, wind_mph=0.0, slope=0.0, **moist):
    fm = pyflam.get_fuel_model(n)
    kw = {**SCENARIO, **moist}
    return pyflam.spread(
        fm, wind_midflame=mph_to_ft_per_min(wind_mph), slope=slope, **kw
    ).rate_of_spread


# --- Physics / property tests -------------------------------------------------

@pytest.mark.parametrize("n", range(1, 14))
def test_baseline_positive_and_finite(n):
    """Every fuel model spreads at a finite, positive rate with no wind/slope."""
    r = _ros(n)
    assert math.isfinite(r)
    assert r > 0.0


@pytest.mark.parametrize("n", range(1, 14))
def test_wind_increases_spread(n):
    assert _ros(n, wind_mph=10) > _ros(n, wind_mph=0)


@pytest.mark.parametrize("n", range(1, 14))
def test_slope_increases_spread(n):
    assert _ros(n, slope=0.50) > _ros(n, slope=0.0)


@pytest.mark.parametrize("n", range(1, 14))
def test_higher_dead_moisture_slows_spread(n):
    dry = _ros(n, m_1h=0.03, m_10h=0.04, m_100h=0.05)
    wet = _ros(n, m_1h=0.20, m_10h=0.21, m_100h=0.22)
    assert wet < dry


def test_zero_at_dead_moisture_of_extinction():
    """At/above the dead moisture of extinction, a dead-only fuel won't spread."""
    fm = pyflam.get_fuel_model(8)  # closed timber litter, mx_dead = 30%
    fb = pyflam.spread(fm, m_1h=0.30, m_10h=0.30, m_100h=0.30)
    assert fb.rate_of_spread == pytest.approx(0.0, abs=1e-9)


def test_wind_factor_monotonic_in_wind():
    fm = pyflam.get_fuel_model(1)
    factors = [
        pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(u), **SCENARIO).wind_factor
        for u in (0, 2, 5, 10, 20)
    ]
    assert factors == sorted(factors)
    assert factors[0] == 0.0


def test_flame_length_tracks_intensity():
    """Byram flame length is monotonic in fireline intensity."""
    fm = pyflam.get_fuel_model(4)  # chaparral, high intensity
    calm = pyflam.spread(fm, **SCENARIO)
    windy = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(10), **SCENARIO)
    assert windy.fireline_intensity > calm.fireline_intensity
    assert windy.flame_length > calm.flame_length


def test_byram_intensity_consistency():
    """I_B = H_A * R / 60 must hold exactly (units self-consistency)."""
    fm = pyflam.get_fuel_model(10)
    fb = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(6), slope=0.2, **SCENARIO)
    assert fb.fireline_intensity == pytest.approx(
        fb.heat_per_unit_area * fb.rate_of_spread / 60.0, rel=1e-12
    )


def test_no_negative_outputs():
    for n in range(1, 14):
        fb = pyflam.spread(
            pyflam.get_fuel_model(n),
            wind_midflame=mph_to_ft_per_min(8), slope=0.3, **SCENARIO,
        )
        for field in (fb.rate_of_spread, fb.reaction_intensity,
                      fb.fireline_intensity, fb.flame_length):
            assert field >= 0.0


# --- Golden-master regression -------------------------------------------------
# Provisional values produced by this implementation at:
#   m_1h=6% m_10h=7% m_100h=8% live_herb=60% live_woody=90%,
#   midflame wind = 5 mph, slope = 0.
# Replace with FlamMap/BehavePlus-verified numbers once cross-checked.
GOLDEN_ROS_FT_MIN = {
    1: 103.28, 2: 40.22, 3: 129.56, 4: 89.72, 5: 29.05,
    6: 37.23, 7: 32.76, 8: 2.23, 9: 9.65, 10: 10.03,
    11: 6.73, 12: 14.63, 13: 17.66,
}


@pytest.mark.parametrize("n,expected", sorted(GOLDEN_ROS_FT_MIN.items()))
def test_golden_master_ros(n, expected):
    assert _ros(n, wind_mph=5) == pytest.approx(expected, abs=0.01)


# --- Fuel load factor ---------------------------------------------------------

def test_load_factor_unity_is_unchanged():
    fm = pyflam.get_fuel_model(10)
    base = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(5), **SCENARIO)
    same = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(5),
                         load_factor=1.0, **SCENARIO)
    assert same.rate_of_spread == pytest.approx(base.rate_of_spread)
    assert same.fireline_intensity == pytest.approx(base.fireline_intensity)


@pytest.mark.parametrize("n", [1, 2, 3, 4, 10, 104, 145])
def test_load_factor_raises_intensity(n):
    """For well-aerated fuels, more load -> more energy and fireline intensity."""
    fm = pyflam.get_fuel_model(n)
    lo = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(8), **SCENARIO)
    hi = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(8),
                       load_factor=1.3, **SCENARIO)
    assert hi.heat_per_unit_area > lo.heat_per_unit_area
    assert hi.fireline_intensity > lo.fireline_intensity
    assert hi.flame_length > lo.flame_length


def test_load_factor_nonmonotonic_for_compact_litter():
    """Compact litter is already dense: extra load pushes packing past the
    Rothermel optimum, so reaction intensity can *fall* -- correct, if surprising.
    The load factor must still change the result (it isn't ignored)."""
    fm = pyflam.get_fuel_model(183)  # TL3, compact conifer litter
    lo = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(8), **SCENARIO)
    hi = pyflam.spread(fm, wind_midflame=mph_to_ft_per_min(8),
                       load_factor=1.3, **SCENARIO)
    assert hi.reaction_intensity != lo.reaction_intensity
    assert hi.packing_ratio > lo.packing_ratio       # denser bed


def test_load_factor_propagates_to_landscape():
    import numpy as np
    n = 5
    ls = pyflam.Landscape(
        fuel_model=np.full((n, n), 10, dtype=int), slope=np.zeros((n, n)),
        cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=n * 30.0,
        slope_units="degrees")
    base = pyflam.basic_fire_behavior(ls, wind_midflame=mph_to_ft_per_min(5),
                                      **SCENARIO)
    boosted = pyflam.basic_fire_behavior(ls, wind_midflame=mph_to_ft_per_min(5),
                                         load_factor=1.3, **SCENARIO)
    assert (boosted["fireline_intensity"] > base["fireline_intensity"]).all()


def test_load_factor_per_fuel_dict():
    """A {fuel: factor} mapping corrects each fuel type independently."""
    import numpy as np
    fuel = np.array([[1, 1, 104], [104, 1, 104]], dtype=int)
    ls = pyflam.Landscape(fuel_model=fuel, slope=np.zeros((2, 3)),
                          cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=60.0,
                          slope_units="degrees")
    out = pyflam.basic_fire_behavior(
        ls, wind_midflame=mph_to_ft_per_min(5),
        load_factor={1: 1.0, 104: 1.3}, **SCENARIO)
    ref1 = pyflam.basic_fire_behavior(ls, wind_midflame=mph_to_ft_per_min(5),
                                      load_factor=1.0, **SCENARIO)
    ref104 = pyflam.basic_fire_behavior(ls, wind_midflame=mph_to_ft_per_min(5),
                                        load_factor=1.3, **SCENARIO)
    f1 = fuel == 1
    f104 = fuel == 104
    assert np.allclose(out["fireline_intensity"][f1], ref1["fireline_intensity"][f1])
    assert np.allclose(out["fireline_intensity"][f104],
                       ref104["fireline_intensity"][f104])


def test_per_cell_moisture_field():
    """basic_fire_behavior accepts a per-cell dead fuel moisture field."""
    import numpy as np
    n = 6
    ls = pyflam.Landscape(fuel_model=np.full((n, n), 104, dtype=int),
                          slope=np.zeros((n, n)), cellsize_x=30.0, cellsize_y=30.0,
                          west=0.0, north=n * 30.0, slope_units="degrees")
    m1 = np.full((n, n), 0.04)
    m1[:, 3:] = 0.12                       # moist east half
    out = pyflam.basic_fire_behavior(
        ls, wind_midflame=mph_to_ft_per_min(5),
        m_1h=m1, m_10h=m1 + 0.01, m_100h=m1 + 0.02,
        m_live_herb=0.6, m_live_woody=0.9)
    dry = out["rate_of_spread"][:, :3]
    wet = out["rate_of_spread"][:, 3:]
    assert np.allclose(dry, dry[0, 0]) and np.allclose(wet, wet[0, 0])
    assert dry[0, 0] > wet[0, 0]           # drier fuel spreads faster


def test_load_factor_per_cell_raster():
    """A per-cell load-correction raster applies cell by cell."""
    import numpy as np
    n = 6
    fuel = np.full((n, n), 104, dtype=int)
    ls = pyflam.Landscape(fuel_model=fuel, slope=np.zeros((n, n)),
                          cellsize_x=30.0, cellsize_y=30.0, west=0.0,
                          north=n * 30.0, slope_units="degrees")
    lf = np.full((n, n), 1.0)
    lf[:, 3:] = 1.3
    out = pyflam.basic_fire_behavior(ls, wind_midflame=mph_to_ft_per_min(5),
                                     load_factor=lf, **SCENARIO)
    left = out["fireline_intensity"][:, :3]
    right = out["fireline_intensity"][:, 3:]
    assert np.allclose(left, left[0, 0])      # uniform within each region
    assert np.allclose(right, right[0, 0])
    assert right[0, 0] > left[0, 0]           # boosted load -> more intensity
