"""Tests for atmospheric forcing: state, fire-input derivation, providers."""

from __future__ import annotations

import numpy as np
import pytest

import pyflam
from pyflam import atmosphere as atm


# --- humidity & equilibrium moisture content ----------------------------------

def test_emc_known_and_monotonic():
    # Hot & dry -> low EMC; cool & humid -> high EMC.
    dry = atm.equilibrium_moisture_content(35.0, 15.0)
    humid = atm.equilibrium_moisture_content(10.0, 80.0)
    assert 2.0 < dry < 6.0
    assert humid > dry
    # Monotonic in RH at fixed T, and falls with T at fixed RH.
    assert (atm.equilibrium_moisture_content(20.0, 70.0)
            > atm.equilibrium_moisture_content(20.0, 30.0))
    assert (atm.equilibrium_moisture_content(40.0, 50.0)
            < atm.equilibrium_moisture_content(5.0, 50.0))


def test_rh_from_dewpoint():
    assert atm.relative_humidity_from_dewpoint(20.0, 20.0) == pytest.approx(100.0)
    assert atm.relative_humidity_from_dewpoint(30.0, 10.0) < 40.0
    assert 0.0 <= atm.relative_humidity_from_dewpoint(35.0, -5.0) <= 100.0


def test_dead_fuel_moisture_structure():
    st = atm.AtmosphericState(wind_speed=5, wind_direction=270,
                              temperature=25, relative_humidity=40)
    m = atm.dead_fuel_moisture(st)
    assert set(m) == {"m_1h", "m_10h", "m_100h"}
    assert m["m_1h"] < m["m_10h"] < m["m_100h"]      # slower classes offset up
    assert all(0.0 < v < 0.5 for v in m.values())


# --- state construction from forecast/reanalysis variables --------------------

def test_from_si_wind_components_and_units():
    st = atm.AtmosphericState.from_si(
        wind_u=-5.0, wind_v=0.0, temperature_K=300.0, relative_humidity=30.0)
    assert st.wind_speed == pytest.approx(5.0)
    assert st.wind_direction == pytest.approx(90.0)   # u<0 -> wind FROM the east
    assert st.temperature == pytest.approx(26.85, abs=0.01)


def test_from_si_dewpoint_to_rh():
    st = atm.AtmosphericState.from_si(
        wind_speed=3.0, wind_direction=180.0, temperature_K=303.15,
        dewpoint_K=283.15)
    assert 0.0 < st.relative_humidity < 60.0


# --- stability / energy flux --------------------------------------------------

def test_obukhov_sign_and_stability():
    base = dict(wind_speed=4.0, wind_direction=270, temperature=25,
                relative_humidity=30)
    unstable = atm.AtmosphericState(**base, sensible_heat_flux=200.0)
    stable = atm.AtmosphericState(**base, sensible_heat_flux=-50.0)
    neutral = atm.AtmosphericState(**base, sensible_heat_flux=0.0)
    assert atm.obukhov_length(unstable) < 0      # daytime convective
    assert atm.obukhov_length(stable) > 0        # nocturnal stable
    assert atm.obukhov_length(neutral) == float("inf")
    assert atm.stability_class(unstable) == "unstable"
    assert atm.stability_class(stable) == "stable"


def test_cape_implies_unstable():
    st = atm.AtmosphericState(wind_speed=4, wind_direction=270, temperature=30,
                              relative_humidity=20, cape=1800.0, cin=10.0)
    assert atm.stability_class(st) == "unstable"


def test_convective_plume_factor():
    calm = atm.AtmosphericState(wind_speed=3, wind_direction=270, temperature=20,
                                relative_humidity=40)
    capey = atm.AtmosphericState(wind_speed=3, wind_direction=270, temperature=30,
                                 relative_humidity=20, cape=3000.0,
                                 sensible_heat_flux=300.0)
    stable = atm.AtmosphericState(wind_speed=3, wind_direction=270, temperature=10,
                                  relative_humidity=70, sensible_heat_flux=-40.0)
    assert atm.convective_plume_factor(capey) > atm.convective_plume_factor(calm)
    assert atm.convective_plume_factor(stable) < 1.0
    assert 0.5 <= atm.convective_plume_factor(capey) <= 3.0   # bounded


# --- integration with the fire model ------------------------------------------

def test_spread_inputs_from_state():
    st = atm.AtmosphericState(wind_speed=10.0, wind_direction=225.0,
                              temperature=30.0, relative_humidity=20.0)
    si = atm.spread_inputs_from_state(st)
    assert set(si) >= {"m_1h", "m_10h", "m_100h", "wind_midflame", "wind_direction"}
    assert si["wind_direction"] == 225.0
    assert si["wind_midflame"] > 0.0
    # plug straight into the surface model
    fb = pyflam.spread(pyflam.get_fuel_model(104),
                       m_live_herb=0.6, m_live_woody=0.9,
                       **{k: si[k] for k in ("m_1h", "m_10h", "m_100h")},
                       wind_midflame=si["wind_midflame"])
    assert fb.rate_of_spread > 0.0


def test_atmospheric_firebrand_physics_scales_with_convection():
    capey = atm.AtmosphericState(wind_speed=5, wind_direction=270, temperature=32,
                                 relative_humidity=15, cape=2500.0,
                                 sensible_heat_flux=300.0)
    base = pyflam.FirebrandPhysics()
    conv = atm.atmospheric_firebrand_physics(capey)
    assert conv.front_length > base.front_length     # unstable -> farther spotting


# --- providers ----------------------------------------------------------------

def test_constant_atmosphere():
    st = atm.AtmosphericState(wind_speed=6, wind_direction=270, temperature=25,
                              relative_humidity=35)
    prov = pyflam.ConstantAtmosphere(st)
    out = prov.state_at(43.5, 11.0)
    assert out.wind_speed == 6 and out.latitude == 43.5


# --- time-lag (Nelson-style) dead fuel moisture -------------------------------

def test_emc_vectorized():
    rh = np.array([10.0, 50.0, 90.0])
    out = atm.equilibrium_moisture_content(np.full(3, 25.0), rh)
    assert out.shape == (3,)
    assert out[0] < out[1] < out[2]                # rises with RH


def test_time_lag_step_approaches_equilibrium():
    # One time-constant toward EMC removes ~63% of the gap.
    m = atm.time_lag_step(0.20, 0.05, dt_hours=1.0, tau_hours=1.0)
    assert m == pytest.approx(0.05 + 0.15 * np.exp(-1.0), abs=1e-6)


def test_dead_fuel_moisture_model_lag_ordering():
    st = atm.AtmosphericState(wind_speed=4, wind_direction=270, temperature=35,
                              relative_humidity=15)
    model = atm.DeadFuelMoistureModel(m_1h=0.12, m_10h=0.12, m_100h=0.12)
    for _ in range(3):
        out = model.update(st, dt_minutes=60)
    emc = atm.equilibrium_moisture_content(35, 15) / 100.0
    # Drying: 1-h nearly at EMC, 100-h barely moved (memory).
    assert out["m_1h"] < out["m_10h"] < out["m_100h"]
    assert abs(out["m_1h"] - emc) < abs(out["m_100h"] - emc)


def test_dead_fuel_moisture_model_equilibrium_init():
    st = atm.AtmosphericState(wind_speed=4, wind_direction=270, temperature=20,
                              relative_humidity=50)
    model = atm.DeadFuelMoistureModel.equilibrium(st)
    assert model.m_1h == model.m_10h == model.m_100h


# --- per-cell atmospheric fields ----------------------------------------------

def test_constant_field_broadcasts():
    import pyflam
    st = atm.AtmosphericState(wind_speed=5, wind_direction=270, temperature=25,
                              relative_humidity=40, sensible_heat_flux=120.0)
    ls = pyflam.Landscape(fuel_model=np.full((6, 8), 104, dtype=int),
                          slope=np.zeros((6, 8)), cellsize_x=30.0, cellsize_y=30.0,
                          west=0.0, north=180.0, slope_units="degrees")
    fld = atm.ConstantAtmosphere(st).field_on(ls)
    assert np.shape(fld.wind_speed) == (6, 8)
    assert np.all(fld.wind_speed == 5)
    si = atm.spread_inputs_from_state(fld)
    assert np.shape(si["m_1h"]) == (6, 8)            # per-cell moisture field


def test_latlon_grid_none_without_crs():
    import pyflam
    ls = pyflam.Landscape(fuel_model=np.zeros((4, 4), int), slope=np.zeros((4, 4)),
                          cellsize_x=30.0, cellsize_y=30.0, west=0.0, north=120.0)
    assert atm.latlon_grid(ls) is None


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("xarray") is None,
    reason="needs xarray")
def test_gridded_field_on_varies_per_cell():
    import pyflam
    import xarray as xr
    lats = np.linspace(43.4, 43.0, 5)
    lons = np.linspace(11.0, 11.5, 5)
    rh = np.tile(np.linspace(60.0, 20.0, 5), (5, 1))      # dry to the east
    ds = xr.Dataset(
        {"u10": (("latitude", "longitude"), np.full((5, 5), -5.0)),
         "v10": (("latitude", "longitude"), np.zeros((5, 5))),
         "t2m": (("latitude", "longitude"), np.full((5, 5), 303.0)),
         "r2": (("latitude", "longitude"), rh)},
        coords={"latitude": lats, "longitude": lons})
    prov = atm.GriddedAtmosphere(ds, {"wind_u": "u10", "wind_v": "v10",
                                      "temperature_K": "t2m", "relative_humidity": "r2"})
    n = 20
    ls = pyflam.Landscape(fuel_model=np.full((n, n), 104, dtype=int),
                          slope=np.zeros((n, n)), cellsize_x=1000.0,
                          cellsize_y=1000.0, west=0.0, north=n * 1000.0,
                          slope_units="degrees")
    lat2d = np.full((n, n), 43.2)
    lon2d = np.tile(np.linspace(11.0, 11.5, n), (n, 1))
    fld = prov.field_on(ls, latlon=(lat2d, lon2d))
    assert fld.relative_humidity[0, 0] > fld.relative_humidity[0, -1]   # drier east
    si = atm.spread_inputs_from_state(fld)
    assert si["m_1h"][0, 0] > si["m_1h"][0, -1]      # moister west


# --- live fetch helpers -------------------------------------------------------

def test_era5_request_structure():
    req = atm.era5_request(date="2024-08-01", time="13:00",
                           area=(44.0, 10.0, 43.0, 12.0))
    assert req["product_type"] == "reanalysis" and req["format"] == "netcdf"
    assert req["date"] == "2024-08-01" and req["time"] == ["13:00"]
    assert req["area"] == [44.0, 10.0, 43.0, 12.0]
    assert "2m_temperature" in req["variable"]


def test_era5_flux_to_watts():
    # A +250 W/m^2 upward sensible flux is stored by ERA5 as accumulated,
    # downward: sshf = -250 * 3600 J/m^2. Convert back to +250 W/m^2 upward.
    assert atm.era5_flux_to_watts(-250.0 * 3600.0) == pytest.approx(250.0)
    assert atm.era5_flux_to_watts(0.0) == pytest.approx(0.0)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("xarray") is None,
    reason="needs xarray")
def test_open_era5_converts_accumulated_flux(tmp_path):
    """open_atmosphere(source='era5') converts accumulated J/m^2 to W/m^2 upward."""
    import xarray as xr
    f = tmp_path / "era5.nc"
    xr.Dataset(
        {"u10": (("latitude", "longitude"), np.full((2, 2), -3.0)),
         "v10": (("latitude", "longitude"), np.zeros((2, 2))),
         "t2m": (("latitude", "longitude"), np.full((2, 2), 305.0)),
         "d2m": (("latitude", "longitude"), np.full((2, 2), 285.0)),
         "sshf": (("latitude", "longitude"), np.full((2, 2), -200.0 * 3600.0))},
        coords={"latitude": [43.0, 43.1], "longitude": [11.0, 11.1]},
    ).to_netcdf(f)
    prov = atm.open_atmosphere(str(f), source="era5")
    st = prov.state_at(43.05, 11.05)
    assert st.sensible_heat_flux == pytest.approx(200.0)   # +W/m^2 upward
    assert atm.stability_class(st) == "unstable"           # daytime heating
    # field_on applies the same conversion per cell.
    ls = pyflam.Landscape(fuel_model=np.full((4, 4), 104, dtype=int),
                          slope=np.zeros((4, 4)), cellsize_x=30.0, cellsize_y=30.0,
                          west=0.0, north=120.0, slope_units="degrees")
    fld = prov.field_on(ls, latlon=(np.full((4, 4), 43.05),
                                    np.full((4, 4), 11.05)))
    assert np.allclose(fld.sensible_heat_flux, 200.0)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("xarray") is None,
    reason="needs xarray")
def test_fetch_era5_uses_cache(tmp_path):
    """A present cache file is opened without any network/cdsapi call."""
    import xarray as xr
    cache = tmp_path / "era5.nc"
    xr.Dataset(
        {"u10": (("latitude", "longitude"), np.full((2, 2), -3.0)),
         "v10": (("latitude", "longitude"), np.zeros((2, 2))),
         "t2m": (("latitude", "longitude"), np.full((2, 2), 300.0)),
         "d2m": (("latitude", "longitude"), np.full((2, 2), 285.0))},
        coords={"latitude": [43.0, 43.1], "longitude": [11.0, 11.1]},
    ).to_netcdf(cache)
    prov = atm.fetch_era5(str(cache), date="2024-08-01", time="13:00",
                          area=(44.0, 10.0, 43.0, 12.0))   # no force -> cache hit
    st = prov.state_at(43.05, 11.05)
    assert st.wind_speed == pytest.approx(3.0)


def test_longitude_wrapping_for_0_360_grids():
    """A 0-360 dataset (GFS) is queried with a negative longitude correctly."""
    pytest.importorskip("xarray")
    import xarray as xr
    lons = np.array([0.0, 120.0, 240.0, 359.0])     # 0-360 convention
    ds = xr.Dataset(
        {"u10": (("latitude", "longitude"), np.tile([1.0, 2.0, 3.0, 4.0], (2, 1))),
         "v10": (("latitude", "longitude"), np.zeros((2, 4))),
         "t2m": (("latitude", "longitude"), np.full((2, 4), 290.0))},
        coords={"latitude": [40.0, 41.0], "longitude": lons})
    prov = atm.GriddedAtmosphere(ds, {"wind_u": "u10", "wind_v": "v10",
                                      "temperature_K": "t2m"})
    # lon=-120 should map to 240 (value 3.0), not to 0 (value 1.0).
    st = prov.state_at(40.0, -120.0)
    assert st.wind_speed == pytest.approx(3.0)


# --- LIVE network fetch (runs only with herbie/cfgrib + network) --------------

@pytest.mark.skipif(
    __import__("importlib").util.find_spec("herbie") is None
    or __import__("importlib").util.find_spec("cfgrib") is None,
    reason="needs herbie + cfgrib")
def test_fetch_gfs_live():
    """Really download a recent GFS run and derive a fire state (no auth needed)."""
    from datetime import datetime, timedelta, timezone
    run = (datetime.now(timezone.utc) - timedelta(days=1)).replace(
        hour=12, minute=0, second=0, microsecond=0).strftime("%Y-%m-%d %H:%M")
    try:
        prov = atm.fetch_gfs(run=run, fxx=0)
    except Exception as exc:                         # network / run not posted
        pytest.skip(f"GFS fetch unavailable: {exc}")
    st = prov.state_at(43.0, 11.0)                    # Tuscany
    assert 0.0 <= st.wind_speed < 80.0
    assert -60.0 < st.temperature < 60.0
    assert 0.0 <= st.wind_direction <= 360.0
    # the derived fire inputs are usable
    si = atm.spread_inputs_from_state(st)
    assert si["wind_midflame"] >= 0.0 and 0.0 < si["m_1h"] < 0.6


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("cdsapi") is None
    or not __import__("os").path.exists(
        __import__("os").path.expanduser("~/.cdsapirc"))
    or __import__("os").environ.get("PYFLAM_LIVE_ERA5") != "1",
    reason="needs cdsapi + ~/.cdsapirc + PYFLAM_LIVE_ERA5=1 (CDS requests queue)")
def test_fetch_era5_live(tmp_path):
    """Really retrieve a small ERA5 slice (opt-in; CDS queues can be slow)."""
    from datetime import datetime, timedelta
    day = (datetime.utcnow() - timedelta(days=7)).strftime("%Y-%m-%d")
    prov = atm.fetch_era5(str(tmp_path / "era5.nc"), date=day, time="12:00",
                          area=(44.0, 10.0, 43.0, 12.0))   # Tuscany box
    st = prov.state_at(43.0, 11.0)
    assert 0.0 <= st.wind_speed < 80.0
    assert -60.0 < st.temperature < 60.0
    # flux converted to W/m^2 upward (daytime -> typically positive)
    assert st.sensible_heat_flux is not None


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("xarray") is None,
    reason="needs xarray")
def test_gridded_atmosphere_nearest_point():
    import xarray as xr
    lats = np.array([40.0, 41.0, 42.0])
    lons = np.array([10.0, 11.0, 12.0])
    ds = xr.Dataset(
        {
            "u10": (("latitude", "longitude"), np.full((3, 3), -4.0)),
            "v10": (("latitude", "longitude"), np.zeros((3, 3))),
            "t2m": (("latitude", "longitude"), np.full((3, 3), 300.0)),
            "r2": (("latitude", "longitude"), np.full((3, 3), 25.0)),
            "cape": (("latitude", "longitude"), np.full((3, 3), 1200.0)),
        },
        coords={"latitude": lats, "longitude": lons},
    )
    prov = atm.GriddedAtmosphere(ds, {
        "wind_u": "u10", "wind_v": "v10", "temperature_K": "t2m",
        "relative_humidity": "r2", "cape": "cape"})
    st = prov.state_at(40.9, 11.2)
    assert st.wind_speed == pytest.approx(4.0)
    assert st.relative_humidity == pytest.approx(25.0)
    assert st.cape == pytest.approx(1200.0)
    assert atm.stability_class(st) == "unstable"
