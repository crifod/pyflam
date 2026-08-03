"""Tests for atmospheric forcing: state, fire-input derivation, providers."""

from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest

import pyflam
from pyflam import atmosphere as atm


# --- timezone normalisation for gridded providers -----------------------------

def test_naive_utc():
    aware = datetime(2021, 7, 25, 12, tzinfo=timezone(__import__("datetime").timedelta(hours=2)))
    naive = atm._naive_utc(aware)
    assert naive.tzinfo is None and naive.hour == 10        # 12:00+02:00 -> 10:00 UTC
    plain = datetime(2021, 7, 25, 12)
    assert atm._naive_utc(plain) is plain                    # already naive: untouched
    assert atm._naive_utc(None) is None


# --- ERA5 zip (new CDS) reading + merge + valid_time + flux transform ----------

def test_open_atmosphere_era5_zip(tmp_path):
    """The new CDS returns ERA5 as a zip of instant+accum NetCDFs on 'valid_time'."""
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    import zipfile

    lat = np.array([42.5, 43.5, 44.5])
    lon = np.array([10.0, 11.0, 12.0])
    vt = np.array(["2021-07-25T12:00:00"], dtype="datetime64[ns]")
    dims = ("valid_time", "latitude", "longitude")
    shape = (1, lat.size, lon.size)

    def da(val):
        return (dims, np.full(shape, val, dtype="float32"))

    instant = xr.Dataset(
        {"u10": da(3.0), "v10": da(-4.0), "t2m": da(300.0), "d2m": da(285.0),
         "sp": da(99000.0), "cape": da(500.0), "blh": da(1200.0)},
        coords={"valid_time": vt, "latitude": lat, "longitude": lon})
    # ERA5 surface fluxes: accumulated J/m^2 over 1 h, positive downward.
    accum = xr.Dataset(
        {"sshf": da(-3600.0 * 150.0), "slhf": da(-3600.0 * 80.0)},
        coords={"valid_time": vt, "latitude": lat, "longitude": lon})

    f_inst = tmp_path / "data_stream-oper_stepType-instant.nc"
    f_acc = tmp_path / "data_stream-oper_stepType-accum.nc"
    instant.to_netcdf(f_inst)
    accum.to_netcdf(f_acc)
    zpath = tmp_path / "era5.zip"
    with zipfile.ZipFile(zpath, "w") as z:
        z.write(f_inst, f_inst.name)
        z.write(f_acc, f_acc.name)

    prov = atm.open_atmosphere(str(zpath), source="era5")
    assert prov.time_name == "valid_time"
    st = prov.state_at(43.5, 11.0, datetime(2021, 7, 25, 12, tzinfo=timezone.utc))
    assert st.temperature == pytest.approx(300.0 - 273.15, abs=1e-3)
    assert st.cape == pytest.approx(500.0)
    assert st.boundary_layer_height == pytest.approx(1200.0)
    # accumulated J/m^2 down -> W/m^2 up: -(-3600*150)/3600 = +150
    assert st.sensible_heat_flux == pytest.approx(150.0, abs=1e-3)
    assert st.latent_heat_flux == pytest.approx(80.0, abs=1e-3)


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


def test_dead_fuel_moisture_model_per_cell_arrays():
    """A gridded state steps the model per cell (arrays), matching scalar cells."""
    T = np.array([[35.0, 20.0], [35.0, 20.0]])
    RH = np.array([[15.0, 60.0], [15.0, 60.0]])
    st = atm.AtmosphericState(wind_speed=np.full((2, 2), 4.0),
                              wind_direction=np.full((2, 2), 270.0),
                              temperature=T, relative_humidity=RH)
    model = atm.DeadFuelMoistureModel(
        m_1h=np.full((2, 2), 0.12), m_10h=np.full((2, 2), 0.12),
        m_100h=np.full((2, 2), 0.12))
    for _ in range(3):
        out = model.update(st, dt_minutes=60)
    assert out["m_1h"].shape == (2, 2)
    # the dry column dries below the humid column across every lag class
    assert np.all(out["m_1h"][:, 0] < out["m_1h"][:, 1])
    # per-cell result equals the equivalent scalar run
    scalar = atm.DeadFuelMoistureModel(m_1h=0.12, m_10h=0.12, m_100h=0.12)
    dry = atm.AtmosphericState(wind_speed=4, wind_direction=270,
                               temperature=35, relative_humidity=15)
    for _ in range(3):
        sc = scalar.update(dry, dt_minutes=60)
    assert out["m_100h"][0, 0] == pytest.approx(sc["m_100h"], abs=1e-9)


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
    # convection fields available in the analysis: HPBL (decoded despite the
    # 'unknown' short name), CIN and surface pressure.
    assert st.boundary_layer_height is not None and st.boundary_layer_height > 0.0
    assert st.cin is not None

    # Heat fluxes are time-mean forecast fields (only fxx >= 3). Fetch a forecast
    # hour and check they come through (positive upward, no transform needed).
    try:
        fc = atm.fetch_gfs(run=run, fxx=6).state_at(43.0, 11.0)
    except Exception as exc:                          # network / run not posted
        pytest.skip(f"GFS forecast hour unavailable: {exc}")
    assert fc.sensible_heat_flux is not None and fc.latent_heat_flux is not None
    assert -200.0 < fc.sensible_heat_flux < 1200.0
    assert fc.boundary_layer_height is not None and fc.boundary_layer_height > 0.0


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


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("xarray") is None,
    reason="needs xarray")
def test_state_at_with_scalar_time_coordinate():
    """A single-forecast field carries time as a scalar coord; passing an explicit
    time must not raise (regression: live GFS fetch returns exactly this layout)."""
    import xarray as xr
    from datetime import datetime
    lats = np.array([40.0, 41.0, 42.0])
    lons = np.array([10.0, 11.0, 12.0])
    ds = xr.Dataset(
        {
            "t2m": (("latitude", "longitude"), np.full((3, 3), 305.0)),
            "r2": (("latitude", "longitude"), np.full((3, 3), 30.0)),
        },
        coords={"latitude": lats, "longitude": lons,
                "time": np.datetime64("2026-06-24T14:00")},   # 0-d scalar coord
    )
    assert "time" in ds.coords and "time" not in ds.dims
    prov = atm.GriddedAtmosphere(ds, {"temperature_K": "t2m",
                                      "relative_humidity": "r2"})
    st = prov.state_at(40.9, 11.2, time=datetime(2026, 6, 24, 14, 0))
    assert st.relative_humidity == pytest.approx(30.0)
    assert st.temperature == pytest.approx(305.0 - 273.15, abs=0.1)


def test_critical_growth_rate_inverts_the_pft():
    """dA/dt_crit must round-trip Tory & Kepert appendix D: FP = alpha * h * w_a * dA/dt."""
    import numpy as np
    from pyflam.atmosphere import critical_growth_rate_grid, capability_margin

    pft_gw, w_a, alpha, h = 139.0, 1.49, 0.7, 15.0e6
    crit = float(critical_growth_rate_grid(np.array([pft_gw]), fuel_load_kg_m2=w_a)[0])

    # feeding the critical rate back through appendix D must return the PFT
    fp_gw = alpha * h * w_a * (crit * 1.0e4 / 3600.0) / 1.0e9
    assert fp_gw == pytest.approx(pft_gw, rel=1e-9)

    # Guissona: 139 GW threshold -> ~3.2 kha/h, and its observed 7869 ha/h clears it
    assert 3100.0 < crit < 3300.0
    assert capability_margin(7869.0, crit) == pytest.approx(0.39, abs=0.02)
    assert capability_margin(358.0, crit) < 0.0        # Santa Coloma does not

    # a richer fuel bed lowers the growth rate the same atmosphere demands
    assert float(critical_growth_rate_grid(np.array([pft_gw]), fuel_load_kg_m2=3.0)[0]) < crit

    # nan PFT (profile too shallow to reach the free-convection height) propagates
    assert np.isnan(critical_growth_rate_grid(np.array([np.nan]))[0])
    assert np.isnan(capability_margin(0.0, crit))      # a fire with no growth has no margin


# --- cost ladder + form field -------------------------------------------------

def _synthetic_column(*, rh=25.0, abl=1800.0, gamma_free=6.0e-3, wind=8.0,
                      sfc_theta=308.0, nlev=40, top_m=12000.0):
    """Well-mixed ABL of depth ``abl`` under a free troposphere with a constant theta lapse.

    Shaped ``(nlev, 1, 1)`` so the gridded diagnostics run on a single column.
    """
    import numpy as np
    z = np.linspace(10.0, top_m, nlev)
    theta = np.where(z <= abl, sfc_theta, sfc_theta + gamma_free * (z - abl))
    ps = 100000.0
    p = ps * np.exp(-z / 8500.0)
    T = theta * (p / ps) ** 0.286
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    e = np.where(z <= abl, rh, max(rh - 15.0, 5.0)) / 100.0 * es
    q = 0.622 * e / np.maximum(p - 0.378 * e, 1.0)
    col = lambda a: a.reshape(-1, 1, 1)
    return dict(height_agl_m=col(z), temperature_k=col(T), pressure_pa=col(p),
                spec_humidity=col(q), wind_u=col(np.full_like(z, wind)),
                wind_v=col(np.zeros_like(z)), surface_pressure_pa=np.full((1, 1), ps),
                abl_m=np.full((1, 1), float(abl)))


def test_density_explicit_firepower_reproduces_eq31():
    """eq 25 with rho written out must collapse onto eq 31 at the density it folds in."""
    from pyflam.atmosphere import (_firepower_required_gw, _PFT_C, _PFT_C_RHO, _PFT_RHO_REF)

    assert _PFT_C_RHO * _PFT_RHO_REF == pytest.approx(_PFT_C, rel=1e-12)

    # Tory & Kepert's Black Saturday 1000 LST case: z_fc 4.8 km, U_ML 20 m/s, dtheta 9 K
    # -> their stated 1240 GW.
    fp = _firepower_required_gw(4800.0, 9.0, 20.0, _PFT_RHO_REF)
    assert float(fp) == pytest.approx(1244.0, abs=1.0)

    # denser air at a lower target costs proportionally more for the same (z, U, dtheta)
    assert float(_firepower_required_gw(4800.0, 9.0, 20.0, 1.10)) == pytest.approx(
        float(fp) * 1.10 / _PFT_RHO_REF, rel=1e-12)


def test_pft_rho_mode_leaves_the_published_form_alone():
    """``rho_mode`` must default to eq 31 verbatim and reject anything it does not implement."""
    import numpy as np
    from pyflam.atmosphere import pyrocb_firepower_threshold_grid, _PFT_RHO_REF

    c = _synthetic_column()
    args = (c["height_agl_m"], c["temperature_k"], c["pressure_pa"], c["spec_humidity"],
            c["wind_u"], c["wind_v"])
    paper = pyrocb_firepower_threshold_grid(*args, surface_pressure_pa=c["surface_pressure_pa"])
    column = pyrocb_firepower_threshold_grid(*args, surface_pressure_pa=c["surface_pressure_pa"],
                                             rho_mode="column")
    assert float(paper["rho_kg_m3"][0, 0]) == pytest.approx(_PFT_RHO_REF)
    # same (z_fc, dtheta, U) either way -- only the density differs
    assert float(paper["z_fc_m"][0, 0]) == pytest.approx(float(column["z_fc_m"][0, 0]))
    ratio = float(column["pft_gw"][0, 0]) / float(paper["pft_gw"][0, 0])
    assert ratio == pytest.approx(float(column["rho_kg_m3"][0, 0]) / _PFT_RHO_REF, rel=1e-9)

    with pytest.raises(ValueError):
        pyrocb_firepower_threshold_grid(*args, surface_pressure_pa=c["surface_pressure_pa"],
                                        rho_mode="whatever")


def test_cost_ladder_is_ordered_and_keeps_nan_as_impossible():
    """escape <= condense <= deep, and an unreachable pyroCb stays nan rather than a big number."""
    import numpy as np
    from pyflam.atmosphere import pyroconvection_cost_grid, PYROCONVECTION_COST_RUNGS

    out = pyroconvection_cost_grid(**_synthetic_column(rh=60.0, abl=900.0, gamma_free=3.0e-3),
                                   fuel_load_kg_m2=1.49)
    esc, con, deep = (float(out[f"cost_{r}_gw"][0, 0]) for r in PYROCONVECTION_COST_RUNGS)
    assert np.isfinite([esc, con, deep]).all()
    assert esc < con <= deep

    # the targets are nested the same way the criteria are
    assert (float(out["z_escape_m"][0, 0]) < float(out["z_condense_m"][0, 0])
            <= float(out["z_deep_m"][0, 0]))

    # a dry, strongly capped column: no beta makes a buoyant cloud to the -20 C level, and the
    # monotone clamp must not launder that nan into "expensive but possible"
    dry = pyroconvection_cost_grid(**_synthetic_column(rh=20.0, abl=2500.0, gamma_free=7.0e-3))
    assert np.isnan(dry["cost_deep_gw"][0, 0])
    assert np.isfinite(dry["cost_escape_gw"][0, 0])


def test_cost_ladder_prices_the_moist_column_lower():
    """A shallow, moist, weakly capped column must cost less on every rung than a dry deep one."""
    import numpy as np
    from pyflam.atmosphere import pyroconvection_cost_grid

    moist = pyroconvection_cost_grid(**_synthetic_column(rh=60.0, abl=900.0, gamma_free=3.0e-3))
    dry = pyroconvection_cost_grid(**_synthetic_column(rh=30.0, abl=2200.0, gamma_free=6.0e-3))
    for rung in ("escape", "condense"):
        assert float(moist[f"cost_{rung}_gw"][0, 0]) < float(dry[f"cost_{rung}_gw"][0, 0])


def test_cost_ha_h_is_the_gw_ladder_through_appendix_d():
    """The ha/h fields must be exactly critical_growth_rate_grid of the GW fields."""
    import numpy as np
    from pyflam.atmosphere import (pyroconvection_cost_grid, critical_growth_rate_grid,
                                   PYROCONVECTION_COST_RUNGS)

    c = _synthetic_column(rh=55.0, abl=1000.0, gamma_free=3.5e-3)
    bare = pyroconvection_cost_grid(**c)
    with_fuel = pyroconvection_cost_grid(**c, fuel_load_kg_m2=2.2)
    for rung in PYROCONVECTION_COST_RUNGS:
        assert f"cost_{rung}_ha_h" not in bare        # no fuel load, no ha/h
        expect = critical_growth_rate_grid(with_fuel[f"cost_{rung}_gw"], fuel_load_kg_m2=2.2)
        assert with_fuel[f"cost_{rung}_ha_h"] == pytest.approx(expect, nan_ok=True)


def test_form_field_is_categorical_and_admits_the_indeterminate_band():
    """Form carries no firepower, no ordinal, and refuses to take a side near ratio 1."""
    import numpy as np
    from pyflam.atmosphere import (pyroconvection_form, pyroconvection_form_grid,
                                   pyroconvection_cost_grid, PYROCONVECTION_FORM_CODE)

    assert pyroconvection_form(0.4) == "resilient"
    assert pyroconvection_form(1.8) == "overshooting"
    assert pyroconvection_form(1.03) == "indeterminate"      # the ABL depth cannot resolve this
    assert pyroconvection_form(0.95) == "indeterminate"
    assert pyroconvection_form(float("nan")) == "undefined"

    codes = pyroconvection_form_grid(np.array([0.4, 0.95, 1.03, 1.8, np.nan]))
    assert codes.tolist() == [PYROCONVECTION_FORM_CODE["resilient"],
                              PYROCONVECTION_FORM_CODE["indeterminate"],
                              PYROCONVECTION_FORM_CODE["indeterminate"],
                              PYROCONVECTION_FORM_CODE["overshooting"],
                              PYROCONVECTION_FORM_CODE["undefined"]]
    # "undefined" (no ratio) is distinct from "indeterminate" (diagnosed, inside the band)
    assert PYROCONVECTION_FORM_CODE["undefined"] != PYROCONVECTION_FORM_CODE["indeterminate"]

    # a wider band swallows a call the default resolves
    assert pyroconvection_form(1.15, band=(0.5, 1.5)) == "indeterminate"

    # the field rides along with the cost ladder without entering any of its rungs
    c = _synthetic_column()
    plain = pyroconvection_cost_grid(**c)
    labelled = pyroconvection_cost_grid(**c, lcl_abl_ratio=np.array([[1.6]]))
    assert int(labelled["form"][0, 0]) == PYROCONVECTION_FORM_CODE["overshooting"]
    assert int(plain["form"][0, 0]) == PYROCONVECTION_FORM_CODE["undefined"]
    for rung in ("escape", "condense", "deep"):
        assert plain[f"cost_{rung}_gw"] == pytest.approx(labelled[f"cost_{rung}_gw"], nan_ok=True)


def test_escape_rung_is_the_briggs_inversion_without_double_counting_alpha():
    """The escape rung must be Briggs' stable rise solved for the source, alpha applied once."""
    import numpy as np
    from pyflam.atmosphere import (pyroconvection_cost_grid, briggs_plume_rise,
                                   critical_growth_rate_grid, _BRIGGS_C)

    # 120 levels to 12 km = 100 m spacing, so the ~240 m entrainment zone is resolved and the
    # rung returns a number rather than declining (see the resolution-guard test below)
    c = _synthetic_column(rh=35.0, abl=1600.0, gamma_free=4.5e-3, wind=8.0, nlev=120)
    out = pyroconvection_cost_grid(**c)
    cost_gw = float(out["cost_escape_gw"][0, 0])
    s = float(out["stability_escape_s2"][0, 0])
    u = float(out["u_escape_ms"][0, 0])
    abl = float(out["z_escape_m"][0, 0])
    assert s > 0.0 and np.isfinite(cost_gw)

    # round trip: feed the cost back through the forward Briggs law and recover the ABL top.
    # briggs_plume_rise takes the *convective* flux, so the total must be scaled by alpha --
    # which is precisely the step the rung must NOT have taken internally.
    alpha = 0.7
    rise = briggs_plume_rise(cost_gw * 1.0e9 * alpha, u, stability_s2=s)
    # the rung uses the column's own rho/T at the ABL top where briggs_plume_rise uses
    # reference values, so allow the resulting offset rather than an exact identity
    assert 0.75 * abl < rise < 1.35 * abl

    # alpha appears exactly once, in the GW -> ha/h conversion. If the rung had already
    # divided by it, this round trip would come back 1/alpha too large.
    ha_h = critical_growth_rate_grid(np.array([cost_gw]), fuel_load_kg_m2=1.49)[0]
    back = alpha * 15.0e6 * 1.49 * (ha_h * 1.0e4 / 3600.0) / 1.0e9
    assert back == pytest.approx(cost_gw, rel=1e-9)

    # the threshold must scale as Briggs says: linear in wind, cubic in the ABL depth
    fast = pyroconvection_cost_grid(**_synthetic_column(rh=35.0, abl=1600.0,
                                                        gamma_free=4.5e-3, wind=16.0,
                                                        nlev=120))
    assert float(fast["cost_escape_gw"][0, 0]) == pytest.approx(2.0 * cost_gw, rel=0.02)


def test_escape_rung_declines_when_the_entrainment_zone_is_unresolved():
    """A gradient across a layer with no level in it is an interpolation, not a measurement."""
    import numpy as np
    from pyflam.atmosphere import (pyroconvection_cost_grid, entrainment_zone_levels_grid,
                                   _EZ_MIN_LEVELS)

    fine = _synthetic_column(rh=35.0, abl=1600.0, gamma_free=4.5e-3, nlev=120, top_m=12000.0)
    coarse = _synthetic_column(rh=35.0, abl=1600.0, gamma_free=4.5e-3, nlev=6, top_m=12000.0)

    # the zone is max(0.15*1600, 100) = 240 m; 120 levels to 12 km resolve it, 6 do not
    n_fine = entrainment_zone_levels_grid(fine["height_agl_m"], fine["abl_m"])
    n_coarse = entrainment_zone_levels_grid(coarse["height_agl_m"], coarse["abl_m"])
    assert int(n_fine[0, 0]) >= _EZ_MIN_LEVELS
    assert int(n_coarse[0, 0]) < _EZ_MIN_LEVELS

    ok = pyroconvection_cost_grid(**fine)
    bad = pyroconvection_cost_grid(**coarse)
    assert np.isfinite(ok["cost_escape_gw"][0, 0])
    assert np.isnan(bad["cost_escape_gw"][0, 0])          # declines rather than fabricating
    assert int(bad["levels_in_ez"][0, 0]) == int(n_coarse[0, 0])

    # the refusal is reported, not silent, and the rungs above are unaffected by it
    assert np.isfinite(bad["cost_condense_gw"][0, 0])
    # and it can be overridden deliberately, for a caller that accepts the interpolation
    forced = pyroconvection_cost_grid(**coarse, ez_min_levels=0)
    assert np.isfinite(forced["cost_escape_gw"][0, 0])


def test_energy_level_collapses_the_geometry_pair_and_leaves_storage_alone():
    """The raster encoding is frozen; the ordinal used for statistics must not count geometry."""
    from pyflam.atmosphere import (PYROCONVECTION_TYPE_LEVEL, PYROCONVECTION_ENERGY_LEVEL,
                                   PYROCONVECTION_TYPES)

    # storage encoding: unchanged, because published GeoTIFFs are read back through it
    assert [PYROCONVECTION_TYPE_LEVEL[t] for t in PYROCONVECTION_TYPES] == [0, 1, 2, 3, 4]

    # energy ordinal: overshooting and resilient are one level, separated only by LCL/ABL
    assert (PYROCONVECTION_ENERGY_LEVEL["overshooting_pyrocu"]
            == PYROCONVECTION_ENERGY_LEVEL["resilient_pyrocu"])
    assert [PYROCONVECTION_ENERGY_LEVEL[t] for t in PYROCONVECTION_TYPES] == [0, 1, 2, 2, 3]
    assert set(PYROCONVECTION_ENERGY_LEVEL) == set(PYROCONVECTION_TYPES)

    # a column that flips overshooting <-> resilient must not move the severity statistic
    a = PYROCONVECTION_ENERGY_LEVEL["overshooting_pyrocu"]
    b = PYROCONVECTION_ENERGY_LEVEL["resilient_pyrocu"]
    assert a - b == 0
    assert PYROCONVECTION_TYPE_LEVEL["resilient_pyrocu"] - \
        PYROCONVECTION_TYPE_LEVEL["overshooting_pyrocu"] == 1     # the step being retired


def test_plume_top_inverts_the_ladder_and_matches_the_validated_scalar_path():
    """The gridded solver must be the one the MISR validation scored, not a lookalike."""
    import os, sys
    import numpy as np
    from pyflam.atmosphere import plume_top_height_grid, PLUME_TOP_FORMS
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
    from escape_calibration import _predicted_top

    rng = np.random.default_rng(0)
    for _ in range(8):
        c = _synthetic_column(rh=float(rng.uniform(20, 65)), abl=float(rng.uniform(700, 2600)),
                              gamma_free=float(rng.uniform(1.5e-3, 9e-3)),
                              wind=float(rng.uniform(2, 18)), nlev=120)
        fp = float(rng.uniform(1, 300))
        col = [c[k][:, 0, 0] for k in ("height_agl_m", "temperature_k", "pressure_pa",
                                       "wind_u", "wind_v")]
        for form in PLUME_TOP_FORMS:
            grid = float(plume_top_height_grid(
                c["height_agl_m"], c["temperature_k"], c["pressure_pa"], c["wind_u"],
                c["wind_v"], firepower_gw=np.array([[fp]]), abl_m=c["abl_m"], form=form)[0, 0])
            assert grid == pytest.approx(
                _predicted_top(*col, float(c["abl_m"][0, 0]), fp, form=form), abs=1.0)

    c = _synthetic_column(rh=35.0, abl=1600.0, gamma_free=4.5e-3, wind=8.0, nlev=120)
    args = (c["height_agl_m"], c["temperature_k"], c["pressure_pa"], c["wind_u"], c["wind_v"])
    top = lambda fp, **kw: float(plume_top_height_grid(
        *args, firepower_gw=np.array([[fp]]), abl_m=c["abl_m"], **kw)[0, 0])

    # more firepower never buys a lower plume, and the floor is the source layer
    heights = [top(fp) for fp in (0.01, 1.0, 5.0, 20.0, 100.0, 500.0)]
    assert all(a <= b for a, b in zip(heights, heights[1:]))
    # A vanishing fire still reaches the top of the *well-mixed* layer: theta there equals
    # theta_ML, so the demand is zero and rising through a neutral layer costs nothing. The
    # first height that costs anything is the base of the cap.
    assert 1400.0 < heights[0] <= 1600.0
    assert heights[-1] > 3000.0                                   # and a big fire punches out

    # a stronger cap costs height for the same fire
    capped = _synthetic_column(rh=35.0, abl=1600.0, gamma_free=1.1e-2, wind=8.0, nlev=120)
    weak = _synthetic_column(rh=35.0, abl=1600.0, gamma_free=1.5e-3, wind=8.0, nlev=120)
    t = lambda col: float(plume_top_height_grid(
        col["height_agl_m"], col["temperature_k"], col["pressure_pa"], col["wind_u"],
        col["wind_v"], firepower_gw=np.array([[50.0]]), abl_m=col["abl_m"])[0, 0])
    assert t(capped) < t(weak)

    # nan firepower propagates rather than defaulting to a height
    assert np.isnan(plume_top_height_grid(*args, firepower_gw=np.array([[np.nan]]),
                                          abl_m=c["abl_m"])[0, 0])
    with pytest.raises(ValueError):
        top(10.0, form="whatever")
