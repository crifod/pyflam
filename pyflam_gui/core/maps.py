"""Interactive Folium maps: click-to-place ignitions and result overlays.

Built on ``streamlit-folium``. Imports of folium are lazy so the rest of the GUI
keeps working when the ``gui`` extra is not installed (the pages fall back to
manual numeric entry).
"""

from __future__ import annotations

import base64
import io

import numpy as np


def _has_folium() -> bool:
    try:
        import folium  # noqa: F401
        import streamlit_folium  # noqa: F401
        return True
    except ImportError:
        return False


def pick_points(aoi, *, key: str, height: int = 420):
    """Click on the map to drop ignition points; returns a list of ``(lat, lon)``.

    Points accumulate in ``st.session_state[key]`` across reruns; a *Clear* button
    resets them. Falls back to ``None`` (caller should offer manual entry) when the
    map component is unavailable.
    """
    import streamlit as st

    store = f"{key}_pts"
    st.session_state.setdefault(store, [])
    if not _has_folium():
        st.info("Interactive picking needs `folium` + `streamlit-folium`.")
        return None

    import folium
    from streamlit_folium import st_folium

    lat0, lon0 = aoi.center
    fmap = folium.Map(location=[lat0, lon0], zoom_start=8, control_scale=True)
    folium.Rectangle(bounds=[[aoi.south, aoi.west], [aoi.north, aoi.east]],
                     color="#d62728", weight=1, fill=False).add_to(fmap)
    for i, (la, lo) in enumerate(st.session_state[store], 1):
        folium.Marker([la, lo], tooltip=f"ignition {i}",
                      icon=folium.Icon(color="red", icon="fire", prefix="fa")).add_to(fmap)
    fmap.fit_bounds([[aoi.south, aoi.west], [aoi.north, aoi.east]])

    out = st_folium(fmap, height=height, use_container_width=True, key=f"{key}_map",
                    returned_objects=["last_clicked"])
    clicked = (out or {}).get("last_clicked")
    if clicked:
        pt = (round(clicked["lat"], 5), round(clicked["lng"], 5))
        if pt not in st.session_state[store]:
            st.session_state[store].append(pt)

    c1, c2 = st.columns([1, 3])
    if c1.button("Clear ignitions", key=f"{key}_clear"):
        st.session_state[store] = []
    c2.caption(f"{len(st.session_state[store])} ignition(s) placed — click the map to add more.")
    return list(st.session_state[store])


def _array_to_png_datauri(array, *, cmap="inferno", vmin=None, vmax=None):
    """Render a 2D array (nan = transparent) to a base64 PNG ``data:`` URI."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.colors import Normalize

    try:
        colormap = matplotlib.colormaps[cmap]          # matplotlib >= 3.5
    except (AttributeError, KeyError):
        import matplotlib.cm as cm
        colormap = cm.get_cmap(cmap)

    a = np.asarray(array, dtype=float)
    finite = np.isfinite(a)
    if vmin is None:
        vmin = float(np.nanmin(a)) if finite.any() else 0.0
    if vmax is None:
        vmax = float(np.nanmax(a)) if finite.any() else 1.0
    norm = Normalize(vmin=vmin, vmax=vmax or 1.0)
    rgba = colormap(norm(a))
    rgba[..., 3] = np.where(finite, 0.75, 0.0)        # transparent where nan
    rgba8 = (rgba * 255).astype(np.uint8)

    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(rgba8).save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _wgs84_bounds(ls):
    """``[[south, west], [north, east]]`` of a landscape in WGS84, or ``None``."""
    from pyflam.atmosphere import latlon_grid
    grid = latlon_grid(ls)
    if grid is None:
        return None
    lat2d, lon2d = grid
    return [[float(lat2d.min()), float(lon2d.min())],
            [float(lat2d.max()), float(lon2d.max())]]


def result_map(aoi, *, ls=None, overlay=None, overlay_cmap="inferno",
               geojson=None, height: int = 480, key: str = "result_map"):
    """Show a result map: optional raster ImageOverlay + GeoJSON perimeter/vectors.

    ``overlay`` is a 2D array on the landscape grid (nan = transparent);
    ``geojson`` an operative FeatureCollection (already in WGS84). No-op message if
    the map component is missing.
    """
    import streamlit as st
    if not _has_folium():
        st.info("Map overlay needs `folium` + `streamlit-folium`.")
        return
    import folium
    from streamlit_folium import st_folium

    lat0, lon0 = aoi.center
    fmap = folium.Map(location=[lat0, lon0], zoom_start=9, control_scale=True)

    bounds = None
    if overlay is not None and ls is not None:
        bounds = _wgs84_bounds(ls)
        if bounds is not None:
            try:
                uri = _array_to_png_datauri(overlay, cmap=overlay_cmap)
                folium.raster_layers.ImageOverlay(
                    image=uri, bounds=bounds, opacity=0.75,
                    mercator_project=True).add_to(fmap)
            except Exception as exc:  # PIL missing etc. — degrade gracefully
                st.caption(f"(raster overlay skipped: {exc})")

    if geojson is not None:
        def _style(feat):
            kind = feat["properties"].get("kind")
            color = {"perimeter": "#ff3300", "arrow": "#3388ff",
                     "sector": "#000000"}.get(kind, "#888888")
            return {"color": color, "weight": 2, "fillOpacity": 0.05}
        folium.GeoJson(geojson, style_function=_style,
                       marker=folium.CircleMarker(radius=3)).add_to(fmap)

    if bounds is not None:
        fmap.fit_bounds(bounds)
    else:
        fmap.fit_bounds([[aoi.south, aoi.west], [aoi.north, aoi.east]])
    st_folium(fmap, height=height, use_container_width=True, key=key,
              returned_objects=[])
