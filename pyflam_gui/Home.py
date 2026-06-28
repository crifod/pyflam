"""pyflam GUI suite -- shared setup (Home).

Pick an **area of interest**, a **landscape** and an **output folder** here; every
page reads these back from the session. The four pages in the sidebar then drive
the pyflam pipelines:

1. Fire-weather / pyroconvection preview
2. Real-time fire propagation
3. Reanalysis of a past event
4. Burn probability
"""

from __future__ import annotations

import os

import pyflam_gui.core  # noqa: F401 -- makes `pyflam` importable (adds src/ to path)
import streamlit as st

from pyflam_gui.core import aoi as aoi_mod
from pyflam_gui.core import landscape as lsc
from pyflam_gui.core import state

st.set_page_config(page_title="pyflam GUI", page_icon="🔥", layout="wide")

st.title("pyflam — fire-behavior pipelines")
st.caption("Open, multiplatform wildfire-behavior modelling. Set up the shared "
           "context below, then open a page from the sidebar.")

# --- 1. Area of interest ------------------------------------------------------
st.header("1 · Area of interest")
col_preset, col_info = st.columns([1, 2])
preset = col_preset.selectbox("Preset", ["(draw on map)"] + list(aoi_mod.PRESETS))
if preset != "(draw on map)":
    state.set_aoi(aoi_mod.preset_aoi(preset))

drawn = aoi_mod.draw_aoi(default=state.get_aoi())
if drawn is not None:
    state.set_aoi(drawn)

current = state.get_aoi()
if current is not None:
    col_info.success(
        f"AOI **{current.label}** — N {current.north:.3f}, W {current.west:.3f}, "
        f"S {current.south:.3f}, E {current.east:.3f}  ·  centre "
        f"{current.center[0]:.3f}, {current.center[1]:.3f}")
else:
    col_info.info("Draw a rectangle/polygon on the map, or pick a preset.")

# --- 2. Landscape -------------------------------------------------------------
st.header("2 · Landscape")
src = st.radio("Source", ["FlamMap .lcp file", "GeoTIFF bands", "Synthetic (offline)"],
               horizontal=True)
try:
    if src == "FlamMap .lcp file":
        path = st.text_input("Path to .lcp", key="lcp_path")
        if path and st.button("Load .lcp"):
            state.set_landscape(lsc.load_lcp(path), f".lcp: {os.path.basename(path)}")
    elif src == "GeoTIFF bands":
        st.caption("Provide one GeoTIFF per band; `fuel_model` and `slope` are required.")
        bands = {}
        for name in ("fuel_model", "slope", "elevation", "aspect", "canopy_cover",
                     "canopy_height", "canopy_base_height", "canopy_bulk_density"):
            p = st.text_input(name, key=f"tif_{name}")
            if p:
                bands[name] = p
        if st.button("Load GeoTIFFs") and "fuel_model" in bands and "slope" in bands:
            state.set_landscape(lsc.load_geotiffs(tuple(bands.items())), "GeoTIFF bands")
    else:
        c1, c2 = st.columns(2)
        n = c1.slider("Grid size (cells)", 40, 200, 120, 10)
        cs = c2.slider("Cell size (m)", 10, 60, 30, 5)
        if st.button("Build synthetic landscape"):
            state.set_landscape(lsc.synthetic_landscape(n=n, cellsize=cs),
                                f"synthetic {n}x{n} @ {cs} m")
except Exception as exc:  # surface load errors instead of crashing the page
    st.error(f"Could not load landscape: {exc}")

ls = state.get_landscape()
if ls is not None:
    st.success(f"Landscape loaded — {st.session_state.get(state.LANDSCAPE_SRC, '')}")
    st.json(lsc.describe(ls))
else:
    st.info("No landscape loaded yet. The preview page can run weather-only without one.")

# --- 3. Output folder ---------------------------------------------------------
st.header("3 · Output folder")
base = st.text_input("Base output directory", value=state.get_output_base())
st.session_state[state.OUTPUT_BASE] = base
st.caption(f"Each run writes a timestamped subfolder under `{base}` "
           "(GeoTIFF + GeoJSON + PNG + PDF).")

st.divider()
st.markdown(
    "**Next:** open a page from the sidebar — "
    "`Fire Weather Preview`, `Realtime Propagation`, `Event Reanalysis`, "
    "or `Burn Probability`.")
