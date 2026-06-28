"""Page 1 — Fire-weather / fire-danger / pyroconvection preview.

Two products over a chosen area:

* **Meteo variation** — samples any atmosphere provider over the run window and
  charts how the fire-weather drivers change (works offline with a constant state).
* **Pyroconvection type** — ICON-2I 2.2 km columns classified into the Castellnou
  et al. (2022) prototypes, as a *potential* (atmosphere-only) and, when a
  landscape is loaded, a *fuel-gated* map; written as PNG + per-hour GeoTIFF (+ PDF).
"""

from __future__ import annotations

import os
from datetime import datetime, timezone

import pyflam_gui.core  # noqa: F401 -- path bootstrap for `pyflam`
import numpy as np
import streamlit as st

from pyflam_gui.core import atmosphere as atmo
from pyflam_gui.core import outputs, plotting, state
from pyflam_gui.core import pyroconv as pc

st.set_page_config(page_title="Fire-Weather Preview", page_icon="🌡️", layout="wide")
st.title("Fire-weather & pyroconvection preview")


def _write_pyro_outputs(run_dir, pres, flip):
    """PNG panels + per-hour GeoTIFFs for the potential/gated maps."""
    import rasterio
    from rasterio.transform import from_origin
    lat, lon = pres["lat"], pres["lon"]
    dlon = float(abs(lon[1] - lon[0])); dlat = float(abs(lat[1] - lat[0]))
    tr = from_origin(lon.min() - dlon / 2, lat.max() + dlat / 2, dlon, dlat)
    for tag, cats in (("potential", pres["pot"]), ("gated", pres["gate"])):
        if cats is None:
            continue
        fig = plotting.pyroconv_panels(cats, pres["bbox"], pres["hours"],
                                       title=f"{tag} — VALID {pres['valid']}",
                                       flip_for_display=flip)
        plotting.save_png(fig, os.path.join(run_dir, f"pyroconv_{tag}.png"))
        for hi, hour in enumerate(pres["hours"]):
            arr = cats[hi] if flip else cats[hi][::-1]
            with rasterio.open(
                os.path.join(run_dir, f"pyroconv_{tag}_{hour:02d}Z.tif"), "w",
                driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                dtype="int16", crs="EPSG:4326", transform=tr, nodata=-1) as ds:
                ds.write(arr.astype("int16"), 1)

if not state.require(aoi=True):
    st.stop()
aoi = state.get_aoi()
st.caption(f"Area of interest: **{aoi.label}** · centre {aoi.center[0]:.3f}, {aoi.center[1]:.3f}")

# =============================================================================
# A. Meteo variation (any source, offline-friendly)
# =============================================================================
st.header("A · Meteo variation over the window")
with st.form("meteo_form"):
    c1, c2, c3 = st.columns(3)
    day = c1.date_input("Start day (UTC)")
    start_hour = c2.number_input("Start hour (Z)", 0, 23, 0)
    total_h = c3.number_input("Window length (h)", 1, 72, 24)
    step_min = st.slider("Sampling step (min)", 30, 360, 180, 30)
    choice = atmo.pick_atmosphere(key="meteo", aoi=aoi)
    go_meteo = st.form_submit_button("Build meteo report")

if go_meteo:
    try:
        from pyflam.meteo_report import meteo_report
        provider = atmo.build_provider(choice)
        start = datetime(day.year, day.month, day.day, int(start_hour), tzinfo=timezone.utc)
        report = meteo_report(provider, location=aoi.center, start_time=start,
                              total_minutes=float(total_h) * 60.0, step_minutes=float(step_min))
        st.session_state["preview_meteo"] = report
    except Exception as exc:
        st.error(f"Meteo report failed: {exc}")

report = st.session_state.get("preview_meteo")
if report is not None:
    st.pyplot(plotting.meteo_timeseries(report))
    st.subheader("Variation summary")
    st.dataframe(report.variation())
    st.text(report.summary())

# =============================================================================
# B. Pyroconvection type (ICON-2I)
# =============================================================================
st.header("B · Pyroconvection type (ICON-2I 2.2 km)")
st.caption("Potential = atmosphere-only upper bound. Fuel-gated (needs a landscape) "
           "= only where modelled fire power ≥ 10 MW/m on your fuels.")
with st.form("pyro_form"):
    c1, c2, c3 = st.columns(3)
    run_day = c1.date_input("Run day (UTC)", key="pyro_run_day")
    run_hh = c2.selectbox("Run hour (Z)", [0, 12], key="pyro_run_hh")
    valid_day = c3.date_input("Valid day (UTC)", key="pyro_valid_day")
    hours = st.multiselect("Valid hours (Z)", list(range(0, 24, 3)),
                           default=[0, 6, 12, 18])
    cache_dir = st.text_input("GRIB cache dir", "/tmp/pyflam_icon2i")
    go_pyro = st.form_submit_button("Fetch ICON-2I & classify")

if go_pyro:
    try:
        from pyflam.atmosphere import fetch_icon2i_mistral, relative_humidity_from_dewpoint
        run_dt = datetime(run_day.year, run_day.month, run_day.day, int(run_hh))
        valid_dt = datetime(valid_day.year, valid_day.month, valid_day.day)
        with st.spinner("Downloading ICON-2I GRIB (first run can be slow)…"):
            files = fetch_icon2i_mistral(run_dt, run=int(run_hh), cache_dir=cache_dir)
        d = pc.read_icon2i(files, aoi.bbox, run_dt, valid_dt, sorted(hours))
        lat, lon, idx = d["lat"], d["lon"], d["idx"]
        RH = relative_humidity_from_dewpoint(d["T2m"], d["Td"])

        ls = state.get_landscape()
        lf = burn = None
        if ls is not None and ls.crs is not None:
            try:
                from pyproj import Transformer
                tr = Transformer.from_crs("EPSG:4326", ls.crs, always_xy=True)
                lf, burn = pc.lcp_fields(ls, lat, lon, tr)
            except Exception as exc:
                st.warning(f"Fuel gate skipped ({exc}); showing potential only.")

        pot, gate = [], []
        prog = st.progress(0.0)
        for j, si in enumerate(idx):
            Tl = {p: d["T"][p][si] for p in pc.DEFAULT_LEVELS}
            pot.append(pc.classify(d["T2m"][si], RH[si], Tl))
            if lf is not None:
                wsp = np.hypot(d["U"][si], d["V"][si])
                gate.append(pc.classify(d["T2m"][si], RH[si], Tl,
                                        fli=pc.fli_grid(d["T2m"][si], RH[si], wsp, lf, burn)))
            prog.progress((j + 1) / len(idx))

        st.session_state["preview_pyro"] = dict(
            pot=np.stack(pot), gate=(np.stack(gate) if gate else None),
            lat=lat, lon=lon, hours=sorted(hours), bbox=aoi.lonlat_bbox,
            valid=valid_dt.strftime("%Y-%m-%d"))
    except Exception as exc:
        st.error(f"Pyroconvection classification failed: {exc}")

pres = st.session_state.get("preview_pyro")
if pres is not None:
    flip = pres["lat"][0] > pres["lat"][-1]
    st.subheader("Potential (atmosphere only)")
    st.pyplot(plotting.pyroconv_panels(
        pres["pot"], pres["bbox"], pres["hours"],
        title=f"Pyroconvection — POTENTIAL — VALID {pres['valid']}",
        flip_for_display=flip))
    if pres["gate"] is not None:
        st.subheader("Fuel-gated (expected on your fuels)")
        st.pyplot(plotting.pyroconv_panels(
            pres["gate"], pres["bbox"], pres["hours"],
            title=f"Pyroconvection — FUEL-GATED — VALID {pres['valid']}",
            flip_for_display=flip))

    if st.button("Write outputs (PNG + GeoTIFF) and bundle"):
        base = state.get_output_base()
        run_dir = outputs.make_run_dir(base, "preview")
        _write_pyro_outputs(run_dir, pres, flip)
        st.success(f"Wrote outputs to `{run_dir}`")
        st.download_button("Download all (.zip)", outputs.zip_dir(run_dir),
                           file_name=os.path.basename(run_dir) + ".zip",
                           mime="application/zip")
