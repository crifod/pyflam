"""Tests for the Canadian Forest Fire Weather Index System (pyflam.fwi).

The anchor is the canonical single-day worked example from Van Wagner & Pickett
(1985): month April, T=17 C, RH=42 %, wind=25 km/h, no rain, starting from
FFMC=85, DMC=6, DC=15 -> FFMC~87.7, DMC~8.5, DC~19.0, ISI~10.9, BUI~8.5, FWI~10.1.
"""

from __future__ import annotations

from datetime import date

import numpy as np
import pytest

from pyflam import fwi


REF = dict(temperature=17.0, relative_humidity=42.0, wind_kmh=25.0,
           rain_mm=0.0, month=4)


def test_canonical_van_wagner_example():
    F = fwi.ffmc(REF["temperature"], REF["relative_humidity"], REF["wind_kmh"],
                 REF["rain_mm"], 85.0)
    P = fwi.dmc(REF["temperature"], REF["relative_humidity"], REF["rain_mm"],
                6.0, REF["month"])
    D = fwi.dc(REF["temperature"], REF["rain_mm"], 15.0, REF["month"])
    I = fwi.isi(REF["wind_kmh"], F)
    U = fwi.bui(P, D)
    S = fwi.fwi(I, U)
    assert F == pytest.approx(87.7, abs=0.2)
    assert P == pytest.approx(8.5, abs=0.2)
    assert D == pytest.approx(19.0, abs=0.3)
    assert I == pytest.approx(10.9, abs=0.2)
    assert U == pytest.approx(8.5, abs=0.2)
    assert S == pytest.approx(10.1, abs=0.3)


def test_system_step_matches_functions():
    sys = fwi.FWISystem(fwi.FWIState(ffmc=85.0, dmc=6.0, dc=15.0))
    out = sys.step(temperature=17.0, relative_humidity=42.0, wind_kmh=25.0,
                   rain_mm=0.0, month=4)
    assert out.fwi == pytest.approx(10.1, abs=0.3)
    assert out.fine_fuel_moisture == pytest.approx(fwi.ffmc_to_moisture(out.ffmc))
    # the carried state advanced to the new codes
    assert sys.state.ffmc == pytest.approx(out.ffmc)
    assert sys.state.dc == pytest.approx(out.dc)


def test_ffmc_moisture_monotone_and_range():
    # drier fuel <-> higher FFMC: moisture strictly decreases with the code
    ms = [fwi.ffmc_to_moisture(F) for F in (30.0, 60.0, 85.0, 95.0)]
    assert ms[0] > ms[1] > ms[2] > ms[3] >= 0.0
    # the two canonical FFMC scale forms are inverses up to the standard ~0.03
    # constant asymmetry (147.2*(101-F) vs 59.5*(250-m))
    for F in (30.0, 60.0, 85.0, 95.0):
        m = fwi.ffmc_to_moisture(F)
        back = 59.5 * (250.0 - m) / (147.2 + m)
        assert back == pytest.approx(F, abs=0.06)


def test_rain_wets_fine_fuels_and_dry_heat_dries():
    dry = fwi.ffmc(30.0, 15.0, 20.0, 0.0, 85.0)     # hot, dry, windy
    wet = fwi.ffmc(15.0, 90.0, 5.0, 15.0, 85.0)     # cool, humid, rain
    assert dry > 85.0 > wet                          # FFMC up when drying, down when wetting
    # drought code climbs on a hot dry day, holds/rises little with rain
    assert fwi.dc(30.0, 0.0, 100.0, 7) > 100.0
    assert fwi.dc(30.0, 20.0, 100.0, 7) < fwi.dc(30.0, 0.0, 100.0, 7)


def test_array_matches_scalar():
    T = np.array([17.0, 30.0])
    H = np.array([42.0, 15.0])
    W = np.array([25.0, 20.0])
    R = np.array([0.0, 0.0])
    Fv = fwi.ffmc(T, H, W, R, np.array([85.0, 85.0]))
    assert Fv.shape == (2,)
    assert Fv[0] == pytest.approx(fwi.ffmc(17.0, 42.0, 25.0, 0.0, 85.0))
    assert Fv[1] == pytest.approx(fwi.ffmc(30.0, 15.0, 20.0, 0.0, 85.0))
    Iv = fwi.isi(W, Fv)
    assert Iv[1] == pytest.approx(fwi.isi(20.0, Fv[1]))


def test_hemisphere_shifts_day_length():
    le_s, lf_s = fwi.southern_hemisphere_tables()
    # July (index 6) in the south uses January's northern factor
    assert le_s[6] == fwi.DMC_DAY_LENGTH[0]
    assert lf_s[6] == fwi.DC_DAY_LENGTH[0]


def test_fwi_from_records_carries_state():
    recs = [dict(temperature=17.0, relative_humidity=42.0, wind_ms=25.0 / 3.6,
                 rain_mm=0.0, month=4)] * 3
    out = fwi.fwi_from_records(recs)
    # DC (deep drought memory) accumulates day over day under continued drying
    assert out[2].dc > out[1].dc > out[0].dc
    assert out[0].fwi == pytest.approx(10.1, abs=0.3)


def test_fwi_from_provider_point():
    from pyflam.atmosphere import AtmosphereProvider, AtmosphericState

    class _Hot(AtmosphereProvider):
        def state_at(self, lat, lon, time=None):
            return AtmosphericState(wind_speed=25.0 / 3.6, wind_direction=270.0,
                                    temperature=17.0, relative_humidity=42.0)

    out = fwi.fwi_from_provider(_Hot(), 43.0, 11.0, [date(2026, 4, 15)])
    assert out[0].fwi == pytest.approx(10.1, abs=0.3)


def test_spinup_state_builds_drought_memory():
    # a 3-week dry, hot spell should leave a high Drought Code (deep-drought memory)
    hot_dry = dict(temperature=33.0, relative_humidity=22.0, wind_ms=4.0,
                   rain_mm=0.0, month=7)
    st = fwi.spinup_state([hot_dry] * 21)
    assert st.dc > 200.0                       # far above the spring-startup DC=15
    assert st.ffmc > 88.0
    # rain during the spell holds the Drought Code down
    wet_spell = [hot_dry] * 10 + [dict(temperature=18.0, relative_humidity=85.0,
                                        wind_ms=3.0, rain_mm=20.0, month=7)] * 11
    st_wet = fwi.spinup_state(wet_spell)
    assert st_wet.dc < st.dc


def test_gridded_fwi_matches_scalar_and_uses_prior_state():
    prev = fwi.FWIState(ffmc=90.0, dmc=40.0, dc=350.0)      # spun-up drought
    T = np.array([[30.0, 24.0], [30.0, 24.0]])
    H = np.array([[20.0, 55.0], [20.0, 55.0]])
    W = np.full((2, 2), 18.0)
    out = fwi.gridded_fwi(temperature=T, relative_humidity=H, wind_kmh=W,
                          month=7, previous=prev)
    assert out["fwi"].shape == (2, 2)
    # the hot dry column is more dangerous than the cool humid column
    assert np.all(out["fwi"][:, 0] > out["fwi"][:, 1])
    # per-cell equals the scalar single-cell computation with the same prior
    F = fwi.ffmc(30.0, 20.0, 18.0, 0.0, prev.ffmc)
    U = fwi.bui(fwi.dmc(30.0, 20.0, 0.0, prev.dmc, 7),
                fwi.dc(30.0, 0.0, prev.dc, 7))
    assert out["fwi"][0, 0] == pytest.approx(fwi.fwi(fwi.isi(18.0, F), U))
    # the large prior DC propagates into a large BUI everywhere
    assert np.all(out["bui"] > 40.0)


def test_gridded_fwi_on_field_state():
    from pyflam.atmosphere import AtmosphericState
    fld = AtmosphericState(
        wind_speed=np.full((3, 3), 5.0), wind_direction=np.full((3, 3), 270.0),
        temperature=np.full((3, 3), 30.0), relative_humidity=np.full((3, 3), 25.0))
    prev = fwi.FWIState(ffmc=88.0, dmc=30.0, dc=300.0)
    out = fwi.gridded_fwi_on(fld, month=7, previous=prev)
    assert out["ffmc"].shape == (3, 3)
    assert np.all(np.isfinite(out["fwi"]))
