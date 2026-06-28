"""Shared run + results panel for the propagation pages (real-time & reanalysis).

Pages 2 and 3 differ only in which weather sources they offer and their default
dates, so the whole body lives here in :func:`render_page`. It covers ignition
selection (map click / manual / centre, multi-point), the up-front run-duration
scenario, crown / plume / pyroconvection toggles, the :func:`pyflam.run_realtime`
call, a metrics table, arrival + perimeter maps, the energy-flux time series and a
full output bundle (GeoTIFF + GeoJSON + meteo CSV + PDF).
"""

from __future__ import annotations

import csv
import os
from datetime import datetime, timedelta, timezone

import numpy as np
import streamlit as st

from . import atmosphere as atmo
from . import maps, outputs, plotting, state


def ignition_rowcol(ls, lat, lon):
    """Map a lat/lon to a landscape ``(row, col)``; centre on failure."""
    import pyflam
    try:
        x, y = lon, lat
        if ls.crs is not None:
            from pyproj import Transformer, CRS
            crs = CRS.from_user_input(ls.crs)
            if not crs.is_geographic:
                x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
        return pyflam.ignition_from_xy(ls, x, y)
    except Exception:
        return (ls.shape[0] // 2, ls.shape[1] // 2)


def propagation_metrics(product, ls) -> dict:
    """Headline metrics from a :class:`pyflam.nrt.RunProduct`."""
    arr = np.asarray(product.arrival_time)
    burned = np.isfinite(arr)
    n = int(burned.sum())
    cell_ha = ls.cellsize_x * ls.cellsize_y / 1.0e4
    out = {
        "burned cells": n,
        "burned area (ha)": round(n * cell_ha, 2),
        "max arrival (min)": round(float(arr[burned].max()), 1) if n else 0.0,
        "run length (min)": product.total_time,
    }
    if product.operative is not None:
        secs = product.operative.sectors
        if secs:
            out["mean head ROS (ft/min)"] = round(
                next((s.mean_ros for s in secs if s.name == "head"), secs[0].mean_ros), 1)
            out["fire heading (deg)"] = round(product.operative.heading, 0)
    return out


def _select_ignitions(ls, aoi, *, key: str) -> list:
    """Ignition picker: map click / manual lat-lon / AOI centre (multi-point)."""
    modes = ["Click on map", "Manual lat/lon", "AOI centre"]
    mode = st.radio("Ignition mode", modes, horizontal=True, key=f"{key}_mode")
    latlons: list = []
    if mode == "Click on map":
        pts = maps.pick_points(aoi, key=key)
        latlons = pts or []
        if not latlons:
            st.caption("No point clicked yet — will use the AOI centre if you run now.")
    elif mode == "Manual lat/lon":
        txt = st.text_area("One `lat, lon` per line",
                           value=f"{aoi.center[0]:.4f}, {aoi.center[1]:.4f}", key=f"{key}_txt")
        for line in txt.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                la, lo = (float(x) for x in line.split(","))
                latlons.append((la, lo))
            except ValueError:
                st.warning(f"Could not parse ignition line: {line!r}")
    if not latlons:
        latlons = [aoi.center]
    return [ignition_rowcol(ls, la, lo) for la, lo in latlons]


def _write_outputs(kind, product, ls, aoi):
    run_dir = outputs.make_run_dir(state.get_output_base(), kind)
    arr = np.asarray(product.arrival_time, dtype=float)
    arr_finite = np.where(np.isfinite(arr), arr, np.nan)       # inf (unburned) -> nodata
    outputs.write_geotiff(ls, arr_finite, os.path.join(run_dir, "arrival_time.tif"))
    burned = np.isfinite(arr).astype("float32")
    outputs.write_geotiff(ls, burned, os.path.join(run_dir, "burned_mask.tif"),
                          dtype="float32", nodata=0)
    geojson = None
    if product.operative is not None:
        geojson = product.operative.to_geojson(ls, to_wgs84=True)
        outputs.write_geojson(product.operative, ls,
                              os.path.join(run_dir, "operative.geojson"), to_wgs84=True)
    # meteo time series CSV
    recs = product.meteo.to_records()
    if recs:
        cols = list(recs[0].keys())
        with open(os.path.join(run_dir, "meteo_timeseries.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(recs)
    # markdown -> PDF report
    md = os.path.join(run_dir, "report.md")
    metrics = propagation_metrics(product, ls)
    with open(md, "w") as fh:
        fh.write(f"# pyflam {kind} report\n\n")
        fh.write(f"AOI **{aoi.label}** · start {product.start_time} · "
                 f"{product.total_time:g} min\n\n## Metrics\n\n")
        for k, v in metrics.items():
            fh.write(f"- **{k}**: {v}\n")
        fh.write("\n## Meteo variation\n\n```\n" + product.meteo.summary() + "\n```\n")
        if product.operative is not None:
            fh.write("\n## Operative analysis\n\n```\n" + product.operative.summary() + "\n```\n")
    pdf = outputs.build_pdf(md, os.path.join(run_dir, "report.pdf"))
    return run_dir, geojson, pdf


def render_page(kind: str, *, allowed_sources, default_day=None, default_hour=10,
                default_total=240):
    """Full propagation page body. ``kind`` is ``"realtime"`` / ``"reanalysis"``."""
    if not state.require(aoi=True, landscape=True):
        st.stop()
    aoi, ls = state.get_aoi(), state.get_landscape()

    st.subheader("1 · Ignition")
    ignitions = _select_ignitions(ls, aoi, key=f"{kind}_ign")

    st.subheader("2 · Run-duration scenario (declared before the run)")
    c1, c2, c3 = st.columns(3)
    day = c1.date_input("Start day (UTC)", value=default_day, key=f"{kind}_day")
    start_hour = c2.number_input("Start hour (Z)", 0, 23, default_hour, key=f"{kind}_hh")
    total_time = c3.number_input("Total run time (min)", 30, 1440, default_total, 30, key=f"{kind}_tt")
    c4, c5 = st.columns(2)
    dt = c4.number_input("Time step dt (min)", 1.0, 120.0, 30.0, 1.0, key=f"{kind}_dt")
    meteo_step = c5.number_input("Meteo sampling step (min)", 15.0, 360.0, 60.0, 15.0, key=f"{kind}_ms")

    st.subheader("3 · Weather & coupling")
    choice = atmo.pick_atmosphere(allowed=allowed_sources, key=f"{kind}_atmo", aoi=aoi)
    cc1, cc2, cc3 = st.columns(3)
    crown = cc1.checkbox("Crown fire", value=True, key=f"{kind}_crown")
    pyro = cc2.checkbox("Pyroconvection feedback", value=False, key=f"{kind}_pyro")
    plume = cc3.checkbox("CFD plume (needs DEM/OpenFOAM)", value=False, key=f"{kind}_plume")
    foliar = st.number_input("Foliar moisture (%) — for crown fire", 30.0, 200.0, 100.0,
                             5.0, key=f"{kind}_fol", disabled=not crown)
    spatial = st.checkbox(
        "Per-cell weather (spatial spread forcing)", value=False, key=f"{kind}_spatial",
        help="Samples the weather per cell for the spread wind. The energy-flux "
             "time series is always sampled at the AOI centre, so it is unaffected. "
             "Not compatible with crown fire in this version.")

    if st.button("Run propagation", type="primary", key=f"{kind}_go"):
        try:
            import pyflam
            provider = atmo.build_provider(choice)
            start = datetime(day.year, day.month, day.day, int(start_hour), tzinfo=timezone.utc)
            run_spatial = spatial
            if spatial and crown:
                st.warning("Crown fire and per-cell weather aren't supported together "
                           "in this version — running single-column so crown fire works.")
                run_spatial = False
            march_kwargs = {}
            if crown:
                march_kwargs["crown"] = True
                march_kwargs["foliar_moisture"] = float(foliar)
            if pyro:
                march_kwargs["pyroconvection"] = True
            with st.spinner("Spreading fire through the evolving weather…"):
                product = pyflam.run_realtime(
                    ls, ignitions, atmosphere=provider, location=aoi.center,
                    start_time=start, total_time=float(total_time), dt=float(dt),
                    meteo_step_minutes=float(meteo_step), spatial=run_spatial,
                    use_plume=plume, **march_kwargs)
            st.session_state[f"{kind}_product"] = product
        except Exception as exc:
            st.error(f"Run failed: {exc}")

    product = st.session_state.get(f"{kind}_product")
    if product is None:
        return

    st.success("Run complete.")
    st.write("### Metrics")
    st.table({k: [v] for k, v in propagation_metrics(product, ls).items()})

    tab_map, tab_energy, tab_op, tab_out = st.tabs(
        ["Maps", "Energy-flux dynamics", "Operative analysis", "Outputs"])

    geojson = None
    if product.operative is not None:
        try:
            geojson = product.operative.to_geojson(ls, to_wgs84=True)
        except Exception:
            geojson = None

    with tab_map:
        c1, c2 = st.columns(2)
        with c1:
            st.caption("Arrival time (min)")
            st.pyplot(plotting.raster(product.arrival_time, title="arrival time (min)",
                                      cmap="inferno", label="min"))
        with c2:
            st.caption("Perimeter & driving forces (interactive)")
            maps.result_map(aoi, ls=ls, overlay=product.arrival_time,
                            overlay_cmap="inferno", geojson=geojson, key=f"{kind}_rmap")

    with tab_energy:
        st.caption("Live now-state and the evolution over the declared window.")
        st.pyplot(plotting.meteo_timeseries(
            product.meteo, variables=["temperature", "relative_humidity", "wind_speed",
                                      "sensible_heat_flux", "latent_heat_flux", "cape",
                                      "boundary_layer_height", "plume_factor"]))
        st.dataframe(product.meteo.variation())

    with tab_op:
        if product.operative is not None:
            st.pyplot(plotting.operative_quiver(product.operative))
            st.text(product.operative.summary())
        else:
            st.info("No perimeter established (fire did not spread).")

    with tab_out:
        if st.button("Write outputs & bundle", key=f"{kind}_write"):
            run_dir, _gj, pdf = _write_outputs(kind, product, ls, aoi)
            st.success(f"Wrote outputs to `{run_dir}`"
                       + ("" if pdf else " (PDF skipped: pandoc not available)"))
            st.download_button("Download all (.zip)", outputs.zip_dir(run_dir),
                               file_name=os.path.basename(run_dir) + ".zip",
                               mime="application/zip", key=f"{kind}_dl")
