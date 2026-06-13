"""Operational analysis of a fire perimeter: what is driving spread, and where.

After a first MTT run, an analyst wants the *operational* picture: split the
perimeter into sectors (head, flanks, tail -- or finer sub-sectors in a complex
scenario) and, for each, the forces pushing the fire out. This module decomposes
the local spread drive into three vectors and their resultant:

1. **Slope** -- the upslope push (Rothermel slope factor phi_s toward upslope,
   from the aspect).
2. **Fuel** -- the push toward more dangerous fuel: the gradient of the intrinsic
   (no wind/slope) spread potential R0, pointing where loads/flammability/burning
   speed increase ahead of the front.
3. **Wind** -- the downwind push (Rothermel wind factor phi_w toward the wind's
   downwind bearing).

Each is reported as a single component vector (its contribution to maximum
spread) and they are summed into a **resultant** -- the net driving force and the
dominant driver per sector. The vectors are the "arrows" a map front-end draws
along the perimeter; this module produces their bearings and magnitudes (and an
optional matplotlib quiver), plus a per-sector text report.

The three magnitudes are dimensionless spread-enhancement terms (phi_w, phi_s,
and a fractional R0 gradient), comparable for situational ranking rather than a
strict physical force balance.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import fuel_models
from .rothermel import kernel_param_groups, surface_kernel

SECTORS = ("head", "right flank", "tail", "left flank")


@dataclass
class ForceVector:
    """A driving-force component: ``magnitude`` and the ``bearing`` it pushes toward."""

    name: str
    magnitude: float
    bearing: float          # degrees clockwise from north (direction of push)

    @property
    def components(self):
        """(east, north) components."""
        r = math.radians(self.bearing)
        return self.magnitude * math.sin(r), self.magnitude * math.cos(r)


@dataclass
class SectorAnalysis:
    name: str
    n_cells: int
    mean_ros: float                  # ft/min over the sector's perimeter cells
    wind: ForceVector
    slope: ForceVector
    fuel: ForceVector
    resultant: ForceVector
    dominant: str                    # 'wind' | 'slope' | 'fuel'
    cells: tuple = ((), ())          # (rows, cols) of this sector's perimeter cells


@dataclass
class OperativeReport:
    heading: float                   # overall fire heading (deg toward)
    centroid: tuple                  # (row, col) of the burned area
    time: float                      # minutes
    sectors: list
    burned: object = None            # boolean burned mask (for polygon tracing)

    def summary(self) -> str:
        lines = [f"Operative analysis at t={self.time:g} min "
                 f"(fire heading {self.heading:.0f} deg):"]
        for s in self.sectors:
            lines.append(
                f"  {s.name:16s} n={s.n_cells:4d} ROS {s.mean_ros:6.1f} ft/min "
                f"| resultant {s.resultant.magnitude:5.2f}@{s.resultant.bearing:3.0f} "
                f"| driver: {s.dominant}")
            lines.append(
                f"      wind {s.wind.magnitude:5.2f}@{s.wind.bearing:3.0f}  "
                f"slope {s.slope.magnitude:5.2f}@{s.slope.bearing:3.0f}  "
                f"fuel {s.fuel.magnitude:5.2f}@{s.fuel.bearing:3.0f}")
        return "\n".join(lines)

    def to_geojson(self, ls, **kwargs) -> dict:
        """Georeferenced GeoJSON (sector points + force arrows); see :func:`to_geojson`."""
        return to_geojson(self, ls, **kwargs)

    def write_geojson(self, ls, path, **kwargs) -> None:
        """Write the GeoJSON to ``path``."""
        write_geojson(self, ls, path, **kwargs)

    def quiver(self, ax=None):  # pragma: no cover - plotting, optional
        """Draw the per-sector force arrows (needs matplotlib). Returns the axes."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots()
        colors = {"wind": "tab:blue", "slope": "tab:green", "fuel": "tab:red",
                  "resultant": "black"}
        for i, s in enumerate(self.sectors):
            for key in ("wind", "slope", "fuel", "resultant"):
                fv = getattr(s, key)
                e, n = fv.components
                ax.quiver(i, 0, e, n, color=colors[key], angles="xy",
                          scale_units="xy", scale=1, label=key if i == 0 else None)
        ax.set_xticks(range(len(self.sectors)))
        ax.set_xticklabels([s.name for s in self.sectors], rotation=30, ha="right")
        ax.legend()
        ax.set_title("Driving-force vectors by sector")
        return ax


def driver_fields(ls, *, m_1h, m_10h, m_100h, m_live_herb=0.0, m_live_woody=0.0,
                  wind_midflame=0.0, wind_direction=0.0, load_factor=1.0):
    """Per-cell driving-force ingredients over a landscape.

    Returns a dict of 2D arrays: ``phi_w`` (wind factor), ``phi_s`` (slope
    factor), ``wind_toward``/``upslope`` (radians), ``r0`` (intrinsic no
    wind/slope ROS = the fuel spread potential), and ``ros_max``/``heading`` (the
    combined max-spread rate and its bearing). Mirrors ``spread_field`` but also
    exposes the wind/slope components and the fuel potential.
    """
    shape = ls.shape
    phi_w = np.zeros(shape)
    phi_s = np.zeros(shape)
    r0 = np.zeros(shape)
    i0 = np.zeros(shape)         # intrinsic (no wind/slope) fireline intensity
    ros_max = np.zeros(shape)
    heading = np.zeros(shape)

    tan_slope = ls.slope_tangent
    fuel = np.asarray(ls.fuel_model)
    wind_mid = np.broadcast_to(np.asarray(wind_midflame, float), shape)
    wind_dir = np.broadcast_to(np.asarray(wind_direction, float), shape)
    wind_toward = np.radians((wind_dir + 180.0) % 360.0)
    if ls.aspect is not None:
        upslope = np.radians((np.asarray(ls.aspect, float) + 180.0) % 360.0)
    else:
        upslope = wind_toward
    moist = dict(m_1h=m_1h, m_10h=m_10h, m_100h=m_100h,
                 m_live_herb=m_live_herb, m_live_woody=m_live_woody)

    for num in np.unique(fuel):
        num = int(num)
        mask = fuel == num
        try:
            fm = fuel_models.get(num)
        except KeyError:
            continue
        if not fm.is_burnable:
            continue
        lf = load_factor if not isinstance(load_factor, dict) \
            else load_factor.get(num, 1.0)
        for sub, p in kernel_param_groups(mask, {"load_factor": lf, **moist}):
            kernel = surface_kernel(fm, **p)
            pw = np.asarray(kernel.wind_factor(wind_mid[sub]), float)
            ps = np.asarray(kernel.slope_factor(tan_slope[sub]), float)
            phi_w[sub] = pw
            phi_s[sub] = ps
            r0[sub] = kernel.r0
            i0[sub] = kernel.heat_per_unit_area * kernel.r0 / 60.0
            vx = pw * np.sin(wind_toward[sub]) + ps * np.sin(upslope[sub])
            vy = pw * np.cos(wind_toward[sub]) + ps * np.cos(upslope[sub])
            ros_max[sub] = kernel.r0 * (1.0 + np.hypot(vx, vy))
            heading[sub] = np.degrees(np.arctan2(vx, vy)) % 360.0

    return {"phi_w": phi_w, "phi_s": phi_s, "wind_toward": wind_toward,
            "upslope": upslope, "r0": r0, "i0": i0,
            "ros_max": ros_max, "heading": heading}


def perimeter_cells(arrival_time, time, burnable=None):
    """Boolean mask of the fire-edge cells at ``time`` (burned, touching unburned)."""
    from scipy.ndimage import binary_erosion
    burned = np.isfinite(arrival_time) & (np.asarray(arrival_time) <= time)
    if not burned.any():
        return burned
    interior = binary_erosion(burned, border_value=0)
    return burned & ~interior


def _circular_mean(bearings_deg, weights):
    r = np.radians(np.asarray(bearings_deg, float))
    w = np.asarray(weights, float)
    s = np.sum(w * np.sin(r))
    c = np.sum(w * np.cos(r))
    return math.degrees(math.atan2(s, c)) % 360.0


def _vector(ex, ny):
    """Build (magnitude, bearing-toward) from mean east/north components."""
    return math.hypot(ex, ny), math.degrees(math.atan2(ex, ny)) % 360.0


def analyze_perimeter(
    ls,
    arrival_time,
    time,
    *,
    drivers=None,
    subsectors: int = 1,
    **driver_kwargs,
) -> OperativeReport:
    """Decompose the spread drive along the perimeter, by operational sector.

    Extracts the perimeter at ``time`` from an MTT ``arrival_time`` raster, splits
    it into head / right flank / tail / left flank (each into ``subsectors``
    angular bins for complex scenarios), and for each sector returns the wind,
    slope and fuel driving-force vectors plus their resultant and the dominant
    driver. ``drivers`` may be a precomputed :func:`driver_fields` dict; otherwise
    pass the run inputs (``wind_midflame``, ``wind_direction``, moistures ...) as
    ``driver_kwargs``.
    """
    if drivers is None:
        drivers = driver_fields(ls, **driver_kwargs)
    edge = perimeter_cells(arrival_time, time)
    rr, cc = np.where(edge)
    if rr.size == 0:
        raise ValueError("no perimeter at this time (nothing burned yet)")

    burned = np.isfinite(arrival_time) & (np.asarray(arrival_time) <= time)
    br, bc = np.where(burned)
    crow, ccol = float(br.mean()), float(bc.mean())

    # Overall heading: ROS-weighted circular mean of the per-cell max-spread az.
    heading = _circular_mean(drivers["heading"][burned], drivers["ros_max"][burned])

    # Each perimeter cell's bearing from the centroid (east=+col, north=-row).
    east = (cc - ccol)
    north = -(rr - crow)
    az_cell = (np.degrees(np.arctan2(east, north))) % 360.0
    delta = (az_cell - heading + 180.0) % 360.0 - 180.0   # [-180, 180]

    # Per-cell force component vectors.
    pw, ps = drivers["phi_w"][rr, cc], drivers["phi_s"][rr, cc]
    wt, up = drivers["wind_toward"][rr, cc], drivers["upslope"][rr, cc]
    wind_e, wind_n = pw * np.sin(wt), pw * np.cos(wt)
    slope_e, slope_n = ps * np.sin(up), ps * np.cos(up)
    # Fuel driver from the intrinsic fireline-intensity potential (rises with
    # load and flammability *and* burning speed); falls back to R0 if absent.
    fuel_potential = drivers.get("i0", drivers["r0"])
    fe, fn, fmag = _fuel_vectors(fuel_potential, ls.cellsize_x, ls.cellsize_y)
    fuel_e, fuel_n = fe[rr, cc], fn[rr, cc]
    ros = drivers["ros_max"][rr, cc]

    sectors = []
    for sec in SECTORS:
        lo, hi = _sector_bounds(sec)
        in_sec = _in_angle(delta, lo, hi)
        for sub in range(subsectors):
            width = (hi - lo) / subsectors
            s_lo, s_hi = lo + sub * width, lo + (sub + 1) * width
            m = in_sec & _in_angle(delta, s_lo, s_hi)
            if not m.any():
                continue
            name = sec if subsectors == 1 else f"{sec} {sub + 1}/{subsectors}"
            sectors.append(_sector_analysis(
                name, m, ros, wind_e, wind_n, slope_e, slope_n,
                fuel_e, fuel_n, rr, cc))

    return OperativeReport(heading=heading, centroid=(crow, ccol), time=float(time),
                           sectors=sectors, burned=burned)


def _fuel_vectors(potential, csx, csy):
    """Fractional gradient of a fuel-hazard ``potential`` (east, north) + magnitude.

    Points toward more dangerous fuel (higher intrinsic fireline intensity = more
    load / flammability / burning speed). The magnitude is the fractional change
    per cell, comparable to phi_w / phi_s.
    """
    gr_row, gr_col = np.gradient(potential)
    denom = np.where(potential > 1e-9, potential, np.inf)
    east = gr_col / denom                    # fractional change per cell, eastward
    north = -gr_row / denom                  # fractional change per cell, northward
    return east, north, np.hypot(east, north)


def _sector_bounds(sector):
    return {"head": (-45.0, 45.0), "right flank": (45.0, 135.0),
            "tail": (135.0, 225.0), "left flank": (-135.0, -45.0)}[sector]


def _in_angle(delta, lo, hi):
    """Whether wrapped angle ``delta`` falls in [lo, hi) (handles the tail wrap)."""
    if hi > 180.0:                            # tail spans the +-180 seam
        return (delta >= lo) | (delta < hi - 360.0)
    return (delta >= lo) & (delta < hi)


def _sector_analysis(name, m, ros, we, wn, se, sn, fe, fn, rr, cc):
    wind = ForceVector("wind", *_vector(we[m].mean(), wn[m].mean()))
    slope = ForceVector("slope", *_vector(se[m].mean(), sn[m].mean()))
    fuel = ForceVector("fuel", *_vector(fe[m].mean(), fn[m].mean()))
    rx = we[m].mean() + se[m].mean() + fe[m].mean()
    ry = wn[m].mean() + sn[m].mean() + fn[m].mean()
    resultant = ForceVector("resultant", *_vector(rx, ry))
    dominant = max((wind, slope, fuel), key=lambda v: v.magnitude).name
    return SectorAnalysis(
        name=name, n_cells=int(m.sum()), mean_ros=float(ros[m].mean()),
        wind=wind, slope=slope, fuel=fuel, resultant=resultant, dominant=dominant,
        cells=(rr[m], cc[m]))


# --- georeferenced export -----------------------------------------------------

def perimeter_rings(burned, ls, *, transform=None, min_vertices=4):
    """Ordered, closed perimeter ring(s) (world coords) traced from a burned mask.

    Uses ``skimage.measure.find_contours`` for subpixel, ordered contours (one
    ring per connected component / hole), converted to the landscape's
    coordinates (or via ``transform``). Returns a list of rings (each a list of
    ``[x, y]``, closed), or ``None`` if scikit-image is unavailable.
    """
    try:
        from skimage import measure
    except ImportError:  # pragma: no cover - exercised only without skimage
        return None
    b = np.asarray(burned, dtype=float)
    if not b.any():
        return []
    padded = np.pad(b, 1)                 # pad so border-touching fires close
    rings = []
    for contour in measure.find_contours(padded, 0.5):
        if len(contour) < min_vertices:
            continue
        ring = []
        for r, c in contour:
            x = ls.west + (c - 1 + 0.5) * ls.cellsize_x   # -1 undoes the pad
            y = ls.north - (r - 1 + 0.5) * ls.cellsize_y
            ring.append(list(transform(x, y)) if transform else [x, y])
        if ring[0] != ring[-1]:
            ring.append(ring[0])
        rings.append(ring)
    return rings


def to_geojson(report, ls, *, arrow_length=1000.0, include_perimeter=True,
               perimeter_geometry="auto", to_wgs84=False) -> dict:
    """Operative analysis as a GeoJSON FeatureCollection for a GIS front-end.

    Emits, in the landscape's coordinates (or WGS84 lon/lat if ``to_wgs84`` and
    ``pyproj`` are available):

    * one **Point** per sector at its perimeter centroid, carrying every force's
      magnitude/bearing, the dominant driver, cell count and mean ROS;
    * one **LineString** "arrow" per force (wind / slope / fuel / resultant) from
      that centroid, its length proportional to magnitude (the longest resultant
      across sectors is scaled to ``arrow_length`` map units);
    * optionally the fire **perimeter** -- a contour-traced ordered
      Polygon/MultiPolygon when ``perimeter_geometry`` is ``"polygon"``/``"auto"``
      and scikit-image is available, else a MultiPoint of edge cells.

    The arrows are exactly the per-sector vectors a map renders along the
    perimeter to show what is driving spread where.
    """
    csx, csy = ls.cellsize_x, ls.cellsize_y
    transform = None
    if to_wgs84 and ls.crs is not None:
        try:
            from pyproj import CRS, Transformer
            transform = Transformer.from_crs(
                CRS.from_user_input(ls.crs), "EPSG:4326", always_xy=True).transform
        except Exception:
            transform = None

    def xy(col, row):
        x = ls.west + (col + 0.5) * csx
        y = ls.north - (row + 0.5) * csy
        return list(transform(x, y)) if transform else [x, y]

    maxmag = max((s.resultant.magnitude for s in report.sectors), default=0.0) or 1.0
    scale = arrow_length / maxmag                    # map units per magnitude unit
    features, all_pts = [], []
    for s in report.sectors:
        rows, cols = s.cells
        rows, cols = np.asarray(rows), np.asarray(cols)
        if rows.size == 0:
            continue
        cr, cc_ = float(rows.mean()), float(cols.mean())
        cxy = xy(cc_, cr)
        all_pts.extend(xy(c, r) for r, c in zip(rows.tolist(), cols.tolist()))
        features.append({
            "type": "Feature", "geometry": {"type": "Point", "coordinates": cxy},
            "properties": {
                "kind": "sector", "sector": s.name, "n_cells": s.n_cells,
                "mean_ros_ft_min": round(s.mean_ros, 2), "dominant": s.dominant,
                **{f"{f}_mag": round(getattr(s, f).magnitude, 3)
                   for f in ("wind", "slope", "fuel", "resultant")},
                **{f"{f}_bearing": round(getattr(s, f).bearing, 1)
                   for f in ("wind", "slope", "fuel", "resultant")}}})
        for f in ("wind", "slope", "fuel", "resultant"):
            fv = getattr(s, f)
            e, n = fv.components
            tip = xy(cc_ + e * scale / csx, cr - n * scale / csy)
            features.append({
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": [cxy, tip]},
                "properties": {"kind": "arrow", "sector": s.name, "force": f,
                               "magnitude": round(fv.magnitude, 3),
                               "bearing": round(fv.bearing, 1)}})

    if include_perimeter:
        rings = None
        if perimeter_geometry in ("polygon", "auto") and report.burned is not None:
            rings = perimeter_rings(report.burned, ls, transform=transform)
        if rings:
            geom = ({"type": "Polygon", "coordinates": [rings[0]]} if len(rings) == 1
                    else {"type": "MultiPolygon", "coordinates": [[r] for r in rings]})
            features.append({"type": "Feature", "geometry": geom,
                             "properties": {"kind": "perimeter",
                                            "time_min": report.time}})
        elif all_pts:                          # fallback: edge cells as points
            features.append({
                "type": "Feature",
                "geometry": {"type": "MultiPoint", "coordinates": all_pts},
                "properties": {"kind": "perimeter", "time_min": report.time}})

    fc = {"type": "FeatureCollection",
          "properties": {"heading": round(report.heading, 1),
                         "time_min": report.time},
          "features": features}
    crs_name = "EPSG:4326" if transform else (str(ls.crs) if ls.crs else None)
    if crs_name:
        fc["crs"] = {"type": "name", "properties": {"name": crs_name}}
    return fc


def write_geojson(report, ls, path, **kwargs) -> None:
    """Write :func:`to_geojson` output to a ``.geojson`` file."""
    import json
    with open(path, "w") as fh:
        json.dump(to_geojson(report, ls, **kwargs), fh)
