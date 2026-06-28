"""Page 4 — Burn probability for a territory.

Grows many fires from random or map-drawn ignitions and reports the per-cell burn
probability plus the connected metrics (FLP, conditional flame length / intensity,
fire-size distribution) via :func:`pyflam.burn_probability`.

The **context scenario** weather is choosable — ERA5 (past), GFS / ICON-2I
(present) or a constant idealized state — and can be a **weather ensemble**: the
same source sampled at several times (with weights) so fire-to-fire weather
variation drives a meaningful flame-length distribution. Optional ember
**spotting** lets fires cross fuel barriers.
"""

from __future__ import annotations

import csv
import os
from datetime import datetime, timezone

import pyflam_gui.core  # noqa: F401
import numpy as np
import pandas as pd
import streamlit as st

from pyflam_gui.core import atmosphere as atmo
from pyflam_gui.core import maps, outputs, plotting, state
from pyflam_gui.core.propagation import ignition_rowcol

st.set_page_config(page_title="Burn Probability", page_icon="🎲", layout="wide")
st.title("Burn probability")

if not state.require(aoi=True, landscape=True):
    st.stop()
aoi, ls = state.get_aoi(), state.get_landscape()


def _burnable_mask(ls):
    import pyflam
    fuel = np.asarray(ls.fuel_model)
    mask = np.zeros(fuel.shape, dtype=bool)
    for num in np.unique(fuel):
        try:
            mask |= (fuel == num) & pyflam.get_fuel_model(int(num)).is_burnable
        except KeyError:
            continue
    return mask


def _random_ignitions(ls, n, seed=0):
    rng = np.random.default_rng(seed)
    rr, cc = np.where(_burnable_mask(ls))
    if rr.size == 0:
        return []
    idx = rng.choice(rr.size, size=min(n, rr.size), replace=False)
    return list(zip(rr[idx].tolist(), cc[idx].tolist()))


# --- ignitions ----------------------------------------------------------------
st.subheader("1 · Ignitions")
ig_mode = st.radio("Source", ["Random", "Draw on map"], horizontal=True)
if ig_mode == "Random":
    c1, c2 = st.columns(2)
    n_fires = c1.number_input("Number of random ignitions", 10, 5000, 200, 10)
    seed = c2.number_input("Random seed", 0, 9999, 0)
    ignitions = _random_ignitions(ls, int(n_fires), int(seed))
else:
    pts = maps.pick_points(aoi, key="bp_ign") or []
    ignitions = [ignition_rowcol(ls, la, lo) for la, lo in pts]
    st.caption(f"{len(ignitions)} ignition(s).")

max_time = st.number_input("Per-fire spread time (min)", 30, 1440, 240, 30)

# --- context weather (single or ensemble) -------------------------------------
st.subheader("2 · Context scenario — weather")
sources = list(atmo.PAST_SOURCES) + list(atmo.PRESENT_SOURCES) + ["constant"]
choice = atmo.pick_atmosphere(allowed=sources, key="bp", aoi=aoi)
st.caption("Sampling times — add rows for a weather **ensemble** (fire-to-fire "
           "variation). Weight sets how many fires draw each scenario.")
default_times = pd.DataFrame({"date": [datetime.now(timezone.utc).date()],
                              "hour": [13], "weight": [1.0]})
times = st.data_editor(default_times, num_rows="dynamic", key="bp_times",
                       use_container_width=True)

c1, c2 = st.columns(2)
m_live_herb = c1.number_input("Live herbaceous moisture (frac)", 0.0, 3.0, 0.7, 0.05)
m_live_woody = c2.number_input("Live woody moisture (frac)", 0.0, 3.0, 0.9, 0.05)

st.subheader("3 · Spotting (optional)")
use_spotting = st.checkbox("Enable ember spotting", value=False)
spot_prob = st.slider("Spot probability", 0.0, 1.0, 0.1, 0.05, disabled=not use_spotting)

go = st.button("Compute burn probability", type="primary")


def _build_field(provider, when, ls):
    """A SpreadField + spotting wind from a provider state at ``when``."""
    import pyflam
    from pyflam import units
    from pyflam.atmosphere import spread_inputs_from_state
    s = provider.state_at(aoi.center[0], aoi.center[1], when)
    si = spread_inputs_from_state(s)
    field = pyflam.spread_field(
        ls, m_1h=si["m_1h"], m_10h=si["m_10h"], m_100h=si["m_100h"],
        m_live_herb=float(m_live_herb), m_live_woody=float(m_live_woody),
        wind_midflame=si["wind_midflame"], wind_direction=si["wind_direction"])
    w20 = units.m_per_s_to_ft_per_min(float(s.wind_speed or 0.0))
    return field, w20, float(s.wind_direction or 0.0)


if go:
    try:
        import pyflam
        if not ignitions:
            raise ValueError("no ignitions (no burnable cells, or none drawn)")
        provider = atmo.build_provider(choice)
        scenarios = []
        for _, row in times.iterrows():
            d = pd.to_datetime(row["date"]).to_pydatetime()
            when = datetime(d.year, d.month, d.day, int(row["hour"]), tzinfo=timezone.utc)
            field, w20, wdir = _build_field(provider, when, ls)
            scenarios.append({"field": field, "weight": float(row["weight"]),
                              "wind_20ft": w20, "wind_direction": wdir})
        field_arg = scenarios[0]["field"] if len(scenarios) == 1 else scenarios

        spotting = wind20 = wdir0 = None
        if use_spotting:
            from pyflam.spotting import SpottingModel
            spotting = SpottingModel(spot_probability=float(spot_prob))
            wind20, wdir0 = scenarios[0]["wind_20ft"], scenarios[0]["wind_direction"]

        kw = dict(max_time=float(max_time), return_metrics=True)
        if spotting is not None:
            kw.update(spotting=spotting, wind_20ft=wind20, wind_direction=wdir0)
        with st.spinner(f"Growing {len(ignitions)} fires"
                        f"{' with spotting' if use_spotting else ''}…"):
            result = pyflam.burn_probability(field_arg, ignitions, **kw)
        st.session_state["bp_result"] = result
        st.session_state["bp_nscen"] = len(scenarios)
    except Exception as exc:
        st.error(f"Burn-probability run failed: {exc}")

result = st.session_state.get("bp_result")
if result is not None:
    st.success(f"Burn probability over {result.n_fires} fires "
               f"({st.session_state.get('bp_nscen', 1)} weather scenario(s)).")
    cell_ha = ls.cellsize_x * ls.cellsize_y / 1.0e4
    burned_frac = float((result.burn_prob > 0).mean())
    st.table({"mean burn prob (burned cells)": [round(float(np.nanmean(
                  result.burn_prob[result.burn_prob > 0])) if burned_frac else 0.0, 3)],
              "burnable area ever reached (%)": [round(100 * burned_frac, 1)],
              "max burn prob": [round(float(np.nanmax(result.burn_prob)), 3)],
              "mean fire size (ha)": [round(float(np.mean(result.fire_sizes)) / 1e4, 2)]})

    tab_map, tab_flame, tab_size, tab_out = st.tabs(
        ["Burn probability", "Flame length", "Fire sizes", "Outputs"])
    bbox = aoi.lonlat_bbox if ls.crs is not None else None

    with tab_map:
        c1, c2 = st.columns(2)
        with c1:
            st.pyplot(plotting.raster(result.burn_prob, title="Burn probability",
                                      cmap="magma", lonlat_bbox=bbox, label="P(burn)"))
        with c2:
            maps.result_map(aoi, ls=ls, overlay=result.burn_prob,
                            overlay_cmap="magma", key="bp_rmap")

    with tab_flame:
        c1, c2 = st.columns(2)
        with c1:
            st.pyplot(plotting.raster(result.conditional_flame_length,
                                      title="Conditional flame length", cmap="inferno",
                                      lonlat_bbox=bbox, label="ft"))
        with c2:
            st.pyplot(plotting.raster(result.conditional_intensity,
                                      title="Conditional fireline intensity",
                                      cmap="inferno", lonlat_bbox=bbox, label="Btu/ft/s"))
        st.pyplot(plotting.flp_bars(result))

    with tab_size:
        st.pyplot(plotting.fire_size_hist(result.fire_sizes, to_ha=1e-4))

    with tab_out:
        if st.button("Write outputs (GeoTIFF + CSV) and bundle", key="bp_write"):
            run_dir = outputs.make_run_dir(state.get_output_base(), "burnprob")
            outputs.write_geotiff(ls, result.burn_prob, os.path.join(run_dir, "burn_prob.tif"))
            outputs.write_geotiff(ls, result.conditional_flame_length,
                                  os.path.join(run_dir, "cond_flame_length.tif"))
            outputs.write_geotiff(ls, result.conditional_intensity,
                                  os.path.join(run_dir, "cond_intensity.tif"))
            with open(os.path.join(run_dir, "fire_sizes_ha.csv"), "w", newline="") as fh:
                w = csv.writer(fh); w.writerow(["fire_size_ha"])
                w.writerows([[s / 1e4] for s in np.asarray(result.fire_sizes).tolist()])
            st.success(f"Wrote outputs to `{run_dir}`")
            st.download_button("Download all (.zip)", outputs.zip_dir(run_dir),
                               file_name=os.path.basename(run_dir) + ".zip",
                               mime="application/zip", key="bp_dl")
