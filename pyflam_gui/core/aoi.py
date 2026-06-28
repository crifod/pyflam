"""Area-of-interest selection: an interactive draw map + pure bbox helpers.

The pure helpers (``bbox_from_geojson``, ``clip_bbox``, ``valid_bbox``) carry no
Streamlit/Folium dependency so they are unit-testable; ``draw_aoi`` renders the
Leaflet map with a rectangle/polygon draw control via ``streamlit-folium`` and
returns the selected :class:`~pyflam_gui.core.state.AreaOfInterest`.
"""

from __future__ import annotations

from .state import AreaOfInterest

# Built-in presets (lon_min, lon_max, lat_min, lat_max) -> handy starting AOIs.
# Tuscany matches the box used by tests/pyroconv_daily.py.
PRESETS: dict[str, tuple[float, float, float, float]] = {
    "Tuscany": (9.6, 12.5, 42.2, 44.6),
    "Italy": (6.5, 18.6, 36.6, 47.1),
    "Catalonia": (0.1, 3.4, 40.5, 42.9),
}


def preset_aoi(name: str) -> AreaOfInterest:
    lon0, lon1, lat0, lat1 = PRESETS[name]
    return AreaOfInterest(north=lat1, west=lon0, south=lat0, east=lon1, label=name)


def valid_bbox(north: float, west: float, south: float, east: float) -> bool:
    """A bbox is valid when north>south, east>west and within WGS84 ranges."""
    return (
        -90.0 <= south < north <= 90.0
        and -180.0 <= west < east <= 180.0
    )


def bbox_from_geojson(geometry: dict) -> tuple[float, float, float, float]:
    """``(north, west, south, east)`` enclosing a GeoJSON Polygon/Rectangle.

    Accepts a GeoJSON geometry dict (as ``streamlit-folium``'s ``Draw`` returns in
    ``last_active_drawing``) with ``coordinates`` of ``[lon, lat]`` rings.
    """
    coords = geometry["coordinates"]
    # Walk to the flat list of [lon, lat] points (Polygon -> [ring], ring -> pts).
    pts = coords
    while pts and isinstance(pts[0], list) and pts[0] and isinstance(pts[0][0], list):
        pts = pts[0]
    lons = [p[0] for p in pts]
    lats = [p[1] for p in pts]
    return (max(lats), min(lons), min(lats), max(lons))


def clip_bbox(
    inner: tuple[float, float, float, float],
    outer: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Clip ``inner`` (n, w, s, e) to ``outer`` (n, w, s, e); both in degrees."""
    n_i, w_i, s_i, e_i = inner
    n_o, w_o, s_o, e_o = outer
    return (min(n_i, n_o), max(w_i, w_o), max(s_i, s_o), min(e_i, e_o))


def aoi_from_geojson(geometry: dict, label: str = "drawn") -> AreaOfInterest:
    n, w, s, e = bbox_from_geojson(geometry)
    return AreaOfInterest(north=n, west=w, south=s, east=e, label=label)


def manual_bbox(default: AreaOfInterest | None = None):
    """Numeric bbox entry -- the fallback when the Folium map is unavailable."""
    import streamlit as st
    d = default or AreaOfInterest(north=44.6, west=9.6, south=42.2, east=12.5)
    c1, c2, c3, c4 = st.columns(4)
    n = c1.number_input("North", -90.0, 90.0, float(d.north), format="%.3f")
    w = c2.number_input("West", -180.0, 180.0, float(d.west), format="%.3f")
    s = c3.number_input("South", -90.0, 90.0, float(d.south), format="%.3f")
    e = c4.number_input("East", -180.0, 180.0, float(d.east), format="%.3f")
    if valid_bbox(n, w, s, e):
        return AreaOfInterest(north=n, west=w, south=s, east=e, label="manual")
    st.warning("Invalid bbox: need north>south and east>west within WGS84 ranges.")
    return None


def draw_aoi(default: AreaOfInterest | None = None, height: int = 480):
    """Render the draw map; return the drawn :class:`AreaOfInterest` or ``None``.

    Streamlit/Folium only -- imported lazily so the pure helpers above stay
    importable without the GUI extra installed. If ``folium`` /
    ``streamlit-folium`` are missing, falls back to :func:`manual_bbox`.
    """
    import streamlit as st
    try:
        import folium
        from folium.plugins import Draw
        from streamlit_folium import st_folium
    except ImportError:
        st.info("Interactive map needs `folium` + `streamlit-folium` "
                "(the `gui` extra). Enter the bounding box manually:")
        return manual_bbox(default)

    if default is not None:
        lat0, lon0 = default.center
        bounds = [[default.south, default.west], [default.north, default.east]]
    else:
        lat0, lon0, bounds = 43.4, 11.0, None

    fmap = folium.Map(location=[lat0, lon0], zoom_start=7, control_scale=True)
    Draw(
        export=False,
        draw_options={"polyline": False, "circle": False, "marker": False,
                      "circlemarker": False, "polygon": True, "rectangle": True},
        edit_options={"edit": False},
    ).add_to(fmap)
    if default is not None:
        folium.Rectangle(bounds=bounds, color="#d62728", weight=2,
                         fill=False, tooltip=f"current AOI ({default.label})").add_to(fmap)
        fmap.fit_bounds(bounds)

    out = st_folium(fmap, height=height, use_container_width=True,
                    returned_objects=["last_active_drawing"])
    drawing = (out or {}).get("last_active_drawing")
    if drawing and drawing.get("geometry"):
        return aoi_from_geojson(drawing["geometry"])
    return None
