# SPDX-License-Identifier: AGPL-3.0-or-later
"""Build a full-field ICON-2I atmosphere provider (wind + temp + dewpoint) and a
pressure-level temperature profile, from the MISTRAL open GRIB archive.

The GUI's build_provider only opens T2M; here we merge U10/V10/T2M/TD2M into one
xarray Dataset so wind, humidity and fuel moisture all come from ICON-2I.
"""
from __future__ import annotations
import datetime as _dt
import numpy as np
import xarray as xr
import pyflam
from pyflam import atmosphere as atm


def _open_one(path, want):
    """Open a single-field ICON GRIB, return (DataArray renamed to `want`)."""
    ds = xr.open_dataset(path, engine="cfgrib",
                         backend_kwargs={"indexpath": ""})
    # the sole data variable in the file
    name = [v for v in ds.data_vars][0]
    da = ds[name]
    # normalise coord names
    ren = {}
    for c in da.coords:
        cl = c.lower()
        if cl in ("latitude", "lat"):
            ren[c] = "latitude"
        elif cl in ("longitude", "lon"):
            ren[c] = "longitude"
    da = da.rename(ren)
    # forecast axis is `step`/`valid_time`; make valid_time the time dimension.
    if "step" in da.dims and "valid_time" in da.coords:
        da = da.swap_dims({"step": "valid_time"})
    # drop the scalar reference-run 'time' coord so valid_time can take the name
    if "time" in da.coords and "time" not in da.dims:
        da = da.drop_vars("time", errors="ignore")
    if "valid_time" in da.coords:
        da = da.rename({"valid_time": "time"})
    # keep only time/latitude/longitude; drop leftover scalar coords.
    keep = {"time", "latitude", "longitude"}
    drop = [c for c in da.coords if c not in keep and c not in da.dims]
    if drop:
        da = da.drop_vars(drop, errors="ignore")
    return da.rename(want), name


# surface + pressure-level fields, incl. RELHUM aloft for a real sounding.
FIELDS = [
    ("T850", "T", "isobaricInhPa", 850), ("T700", "T", "isobaricInhPa", 700),
    ("T500", "T", "isobaricInhPa", 500),
    ("RH850", "RELHUM", "isobaricInhPa", 850),
    ("RH700", "RELHUM", "isobaricInhPa", 700),
    ("RH500", "RELHUM", "isobaricInhPa", 500),
    ("T2M", "T_2M", "heightAboveGround", 2),
    ("TD2M", "TD_2M", "heightAboveGround", 2),
    ("U10", "U_10M", "heightAboveGround", 10),
    ("V10", "V_10M", "heightAboveGround", 10),
]


def build_icon2i(date, run=0, cache_dir="/tmp/pyflam_icon2i"):
    files = pyflam.fetch_icon2i_mistral(date, run=run, cache_dir=cache_dir,
                                        fields=FIELDS)
    das = {}
    raw = {}
    for local, want in (("U10", "wind_u"), ("V10", "wind_v"),
                        ("T2M", "temperature_K"), ("TD2M", "dewpoint_K")):
        da, rawname = _open_one(files[local], want)
        das[want] = da
        raw[want] = rawname
    ds = xr.Dataset(das)
    var_map = {"wind_u": "wind_u", "wind_v": "wind_v",
               "temperature_K": "temperature_K", "dewpoint_K": "dewpoint_K"}
    prov = atm.GriddedAtmosphere(ds, var_map, lat_name="latitude",
                                 lon_name="longitude", time_name="time")
    # coarse vertical sounding: T + RH at 850/700/500 (real ICON aloft humidity).
    levels = {}
    for tloc, rhloc, hpa in (("T850", "RH850", 850), ("T700", "RH700", 700),
                             ("T500", "RH500", 500)):
        try:
            tda, _ = _open_one(files[tloc], "t")
            rda, _ = _open_one(files[rhloc], "rh")
            levels[hpa] = (tda, rda)
        except Exception:
            pass
    return prov, ds, levels, files


def profile_at(prov, levels, surface_state, lat, lon, time):
    """Build a coarse AtmosphericProfile (surface + 850/700/500) at a point/time.

    Surface from the ICON 2 m state; aloft T + RH from the pressure-level grids.
    """
    from pyflam import AtmosphericProfile
    from pyflam.atmosphere import _naive_utc
    t = _naive_utc(time)
    pres = [980.0]                       # ~surface (elev ~250 m)
    temp = [float(surface_state.temperature)]
    rh = [float(surface_state.relative_humidity)]
    for hpa in (850, 700, 500):
        if hpa not in levels:
            continue
        tda, rda = levels[hpa]
        tv = float(tda.sel(latitude=lat, longitude=lon % 360.0 if float(tda.longitude.max()) > 180 else lon,
                           method="nearest").sel(time=t, method="nearest")) - 273.15
        rv = float(rda.sel(latitude=lat, longitude=lon % 360.0 if float(rda.longitude.max()) > 180 else lon,
                           method="nearest").sel(time=t, method="nearest"))
        pres.append(float(hpa)); temp.append(tv); rh.append(max(1.0, min(100.0, rv)))
    return AtmosphericProfile.from_rh(np.array(pres), np.array(temp), np.array(rh))


if __name__ == "__main__":
    import sys
    lat, lon = 43.93561, 11.09635
    day = _dt.date(2026, 7, 3)
    run = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    prov, ds, levels, files = build_icon2i(day, run=run)
    print("ICON-2I files:", {k: v.split('/')[-1] for k, v in files.items()})
    print("dataset vars:", list(ds.data_vars), "dims:", dict(ds.sizes))
    print("time coord:", np.asarray(ds["time"].values).ravel()[:3])
    # sample the state through the run window (14:00-20:00 local = 12:00-18:00 UTC)
    for h in range(12, 19):
        t = _dt.datetime(2026, 7, 3, h, tzinfo=_dt.timezone.utc)
        try:
            st = prov.state_at(lat, lon, t)
            m = atm.dead_fuel_moisture(st)
            print(f"{h:02d}Z  wind {st.wind_speed:4.1f} m/s @ {st.wind_direction:3.0f}  "
                  f"T {st.temperature:4.1f}C  RH {st.relative_humidity:4.1f}%  "
                  f"m1h {m['m_1h']*100:4.1f}%  m10h {m['m_10h']*100:4.1f}%  "
                  f"m100h {m['m_100h']*100:4.1f}%")
        except Exception as e:
            print(f"{h:02d}Z  state_at failed: {type(e).__name__}: {e}")
    print("pressure levels fetched:", sorted(levels))
    from pyflam import convective_plume_factor, continuous_haines, inverted_v, lcl_height_m
    t0 = _dt.datetime(2026, 7, 3, 12, tzinfo=_dt.timezone.utc)   # 14:00 local
    st0 = prov.state_at(lat, lon, t0)
    prof = profile_at(prov, levels, st0, lat, lon, t0)
    ch = continuous_haines(prof)
    iv, depr, midrh = inverted_v(prof)
    lcl = lcl_height_m(st0.temperature, st0.relative_humidity)
    pf = convective_plume_factor(st0, profile=prof)
    print(f"\nSounding @14:00 local: Continuous-Haines {ch:.1f}  inverted-V {iv}  "
          f"sfc dewpt-depr {depr:.1f}C  mid-RH {midrh:.0f}%  LCL {lcl:.0f} m")
    print(f"convective_plume_factor (with profile) = {pf:.2f}  "
          f"(vs {convective_plume_factor(st0):.2f} surface-only)")
