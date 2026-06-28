"""Atmosphere-source picker: turn GUI selections into a pyflam provider.

``build_provider`` is the pure mapping from a *choice dict* to a
:class:`pyflam.atmosphere.AtmosphereProvider`; ``pick_atmosphere`` renders the
Streamlit widgets that produce that dict. Splitting them keeps the provider
construction unit-testable (the ``constant`` path needs no network) and lets each
page restrict which sources are offered:

* Page 1 (preview): any source.
* Page 2 (real-time):   GFS / ICON-2I only (present-date forcing).
* Page 3 (reanalysis):  ERA5 only (past real weather).
* Page 4 (burn prob):   ERA5 (past) or GFS / ICON-2I (present) for the context.
"""

from __future__ import annotations

from datetime import datetime

# Canonical source ids -> label shown in the UI.
SOURCES = {
    "icon2i": "ICON-2I 2.2 km (MISTRAL, open)",
    "gfs": "GFS 0.25 deg (NOAA, open)",
    "era5_file": "ERA5 reanalysis (downloaded file)",
    "era5_fetch": "ERA5 reanalysis (fetch from Copernicus CDS)",
    "constant": "Constant / idealized (offline)",
}

# Dependency / credential notes surfaced next to each source.
NOTES = {
    "icon2i": "Needs `xarray`+`cfgrib` (the `atmos` extra). Open data, no login.",
    "gfs": "Needs `herbie-data`+`cfgrib`. Open data, no login.",
    "era5_file": "Needs `xarray`. Point at a NetCDF/GRIB you already downloaded.",
    "era5_fetch": "Needs `cdsapi` and a `~/.cdsapirc` with CDS credentials.",
    "constant": "No network; one uniform state everywhere. Good for smoke tests.",
}

# Which sources each page may offer (see module docstring).
PRESENT_SOURCES = ("gfs", "icon2i")
PAST_SOURCES = ("era5_file", "era5_fetch")


def build_provider(choice: dict):
    """Construct an :class:`AtmosphereProvider` from a ``pick_atmosphere`` dict.

    ``choice`` carries a ``"source"`` id plus the parameters that source needs.
    Heavy/optional imports happen lazily inside each branch so the ``constant``
    path works with only the base install.
    """
    from pyflam import atmosphere as atm

    src = choice["source"]
    if src == "constant":
        state = atm.AtmosphericState(
            wind_speed=float(choice.get("wind_speed", 5.0)),
            wind_direction=float(choice.get("wind_direction", 270.0)),
            temperature=float(choice.get("temperature", 30.0)),
            relative_humidity=float(choice.get("relative_humidity", 25.0)),
            cape=choice.get("cape"),
            boundary_layer_height=choice.get("boundary_layer_height"),
            sensible_heat_flux=choice.get("sensible_heat_flux"),
            latent_heat_flux=choice.get("latent_heat_flux"),
            time=choice.get("time"),
        )
        return atm.ConstantAtmosphere(state)

    if src == "icon2i":
        date = choice["date"]
        if isinstance(date, str):
            date = datetime.strptime(date, "%Y-%m-%d")
        files = atm.fetch_icon2i_mistral(
            date, run=int(choice.get("run", 0)),
            cache_dir=choice.get("cache_dir", "."))
        # fetch_icon2i_mistral returns a dict of GRIB paths; the GUI's pyroconv
        # path consumes those directly. For the single-column providers we open
        # the surface fields as a gridded provider.
        return atm.open_atmosphere(files["T2M"], source="gfs") \
            if isinstance(files, dict) else files

    if src == "gfs":
        return atm.fetch_gfs(run=choice["run"], fxx=int(choice.get("fxx", 0)),
                             cache_dir=choice.get("cache_dir"))

    if src == "era5_file":
        return atm.open_atmosphere(choice["path"], source="era5")

    if src == "era5_fetch":
        return atm.fetch_era5(
            choice["cache_path"], date=choice["date"], time=choice["time"],
            area=choice["area"], force=bool(choice.get("force", False)))

    raise ValueError(f"unknown atmosphere source {src!r}")


def pick_atmosphere(allowed=None, *, key: str = "atmo", aoi=None) -> dict:
    """Render the source picker and return a ``build_provider`` choice dict.

    ``allowed`` restricts the offered sources (default: all). ``aoi`` (an
    :class:`AreaOfInterest`) seeds the ERA5 fetch area.
    """
    import streamlit as st

    ids = list(allowed or SOURCES.keys())
    src = st.selectbox(
        "Weather / reanalysis source", ids,
        format_func=lambda s: SOURCES[s], key=f"{key}_src")
    st.caption(NOTES[src])
    choice: dict = {"source": src}

    if src == "constant":
        c1, c2 = st.columns(2)
        choice["wind_speed"] = c1.number_input("Wind speed (m/s)", 0.0, 60.0, 5.0, key=f"{key}_ws")
        choice["wind_direction"] = c2.number_input("Wind dir (deg FROM)", 0.0, 360.0, 270.0, key=f"{key}_wd")
        choice["temperature"] = c1.number_input("Temperature (C)", -20.0, 55.0, 30.0, key=f"{key}_t")
        choice["relative_humidity"] = c2.number_input("Relative humidity (%)", 0.0, 100.0, 25.0, key=f"{key}_rh")
        with st.expander("Convection / energy-flux (optional)"):
            choice["cape"] = st.number_input("CAPE (J/kg)", 0.0, 6000.0, 800.0, key=f"{key}_cape")
            choice["boundary_layer_height"] = st.number_input("Boundary-layer height (m)", 0.0, 6000.0, 1800.0, key=f"{key}_blh")
            choice["sensible_heat_flux"] = st.number_input("Sensible heat flux (W/m2)", -100.0, 800.0, 250.0, key=f"{key}_shf")
    elif src == "icon2i":
        c1, c2 = st.columns(2)
        choice["date"] = str(c1.date_input("Run date (UTC)", key=f"{key}_date"))
        choice["run"] = c2.selectbox("Run hour (Z)", [0, 12], key=f"{key}_run")
        choice["cache_dir"] = st.text_input("GRIB cache dir", "/tmp/pyflam_icon2i", key=f"{key}_cache")
    elif src == "gfs":
        c1, c2 = st.columns(2)
        choice["run"] = str(c1.text_input("Run (YYYY-MM-DD HH:MM, UTC)", key=f"{key}_run"))
        choice["fxx"] = c2.number_input("Forecast hour", 0, 384, 6, key=f"{key}_fxx")
        st.caption("Surface heat fluxes are time-mean forecast fields — use a "
                   "forecast hour ≥ 3 to get sensible/latent flux (the f000 "
                   "analysis has BLH/CAPE/CIN but no fluxes).")
        choice["cache_dir"] = st.text_input("Herbie cache dir (optional)", "", key=f"{key}_cache") or None
    elif src == "era5_file":
        choice["path"] = st.text_input("Path to downloaded ERA5 file (.nc/.grib)", key=f"{key}_path")
    elif src == "era5_fetch":
        c1, c2 = st.columns(2)
        choice["date"] = str(c1.date_input("Date", key=f"{key}_date"))
        choice["time"] = c2.text_input("Hour (HH:MM)", "12:00", key=f"{key}_time")
        choice["cache_path"] = st.text_input("Cache file path", "/tmp/era5_pyflam.nc", key=f"{key}_cache")
        if aoi is not None:
            choice["area"] = list(aoi.bbox)
            st.caption(f"Fetch area (N,W,S,E): {choice['area']}")
    return choice
