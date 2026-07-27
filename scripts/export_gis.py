# SPDX-License-Identifier: AGPL-3.0-or-later
"""Export the run as GIS layers: isochrone contour shapefile + GeoTIFF rasters.

All layers are in the landscape CRS (EPSG:3035). Reads the march checkpoint
(saved by run_fire_11_09_2pm.py) and re-derives the georeferencing from the same
100 m clipped landscape.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.environ.setdefault("PYFLAM_NOCFD", "1")
# Respect the caller's AOI config (must match the run that wrote the checkpoint);
# fall back to the production DS=1 / asymmetric-margin framing when unset.
for _k, _v in dict(PYFLAM_DS="1", PYFLAM_MARGIN_W="16", PYFLAM_MARGIN_E="6",
                   PYFLAM_MARGIN_N="9", PYFLAM_MARGIN_S="11").items():
    os.environ.setdefault(_k, _v)
import run_fire_11_09_2pm as R
from pyflam.operative import perimeter_rings

OUT = R.OUT
CRS = "EPSG:3035"
ls, ign, xy = R.clip_landscape()
d = np.load(os.path.join(OUT, "march_checkpoint.npz"))
arr = d["arrival"]; times = [int(t) for t in d["times"]]
ros = d["ros_max"]; fli = d["fli"]                      # (nsteps, ny, nx), ft/min & Btu/ft/s
fire_type = d["fire_type"]; crown = d["crown_fraction"]

gis = os.path.join(OUT, "gis")
os.makedirs(gis, exist_ok=True)

# ---- rasters -----------------------------------------------------------------
import rasterio
from rasterio.transform import from_origin
transform = from_origin(ls.west, ls.north, ls.cellsize_x, ls.cellsize_y)


def wtif(name, data, dtype="float32", nodata=-9999.0):
    a = np.asarray(data, dtype=dtype)
    path = os.path.join(gis, name)
    with rasterio.open(path, "w", driver="GTiff", height=a.shape[0], width=a.shape[1],
                       count=1, dtype=dtype, crs=CRS, transform=transform,
                       nodata=nodata, compress="deflate") as dst:
        dst.write(a, 1)
    print("wrote", name)


burned = np.isfinite(arr)
FT_M = 0.3048
from pyflam import units
# peak-over-time behaviour experienced at each cell
ros_peak = np.nanmax(ros, axis=0) * FT_M                                  # m/min
fli_peak = units.btu_per_ft_s_to_kw_per_m(np.nanmax(fli, axis=0))         # kW/m
flame_peak = 0.0775 * np.power(np.clip(fli_peak, 0, None), 0.46)          # m (Byram)

wtif("arrival_time_min.tif", np.where(burned, arr, -9999.0))
wtif("ros_peak_mmin.tif", np.where(burned, ros_peak, -9999.0))
wtif("fireline_intensity_peak_kwm.tif", np.where(burned, fli_peak, -9999.0))
wtif("flame_length_peak_m.tif", np.where(burned, flame_peak, -9999.0))
if fire_type.size:
    wtif("fire_type_final.tif", np.where(burned, fire_type, -9999.0),
         dtype="int16", nodata=-9999)     # 0 surface / 1 passive / 2 active crown
if crown.size:
    wtif("crown_fraction_final.tif", np.where(burned, crown, -9999.0))

# ---- isochrone contour shapefile --------------------------------------------
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString

recs = []
for t in times:
    mask = np.isfinite(arr) & (arr <= t)
    rings = perimeter_rings(mask, ls) or []               # world coords (EPSG:3035)
    lines = [LineString(r) for r in rings if len(r) >= 2]
    if not lines:
        continue
    geom = lines[0] if len(lines) == 1 else MultiLineString(lines)
    clock = (R.START + __import__("datetime").timedelta(minutes=t + 120)).strftime("%H:%M")
    recs.append({"t_min": t, "clock_cest": clock, "geometry": geom})

gdf = gpd.GeoDataFrame(recs, crs=CRS)
shp = os.path.join(gis, "isochrones_30min.shp")
gdf.to_file(shp)                          # writes .shp/.shx/.dbf/.prj/.cpg
print("wrote isochrones_30min.shp  (%d isochrones)" % len(gdf))

# ---- ignition point shapefile (regenerated so it tracks the current run) ------
from shapely.geometry import Point
ign_gdf = gpd.GeoDataFrame(
    [{"lat": R.LAT, "lon": R.LON, "row": int(ign[0]), "col": int(ign[1]),
      "start_cest": "2026-07-03 14:00", "geometry": Point(xy[0], xy[1])}], crs=CRS)
ign_gdf.to_file(os.path.join(gis, "ignition.shp"))
print("wrote ignition.shp  (%.5f, %.5f)" % (R.LAT, R.LON))
print("\nAll GIS layers in:", gis)
for f in sorted(os.listdir(gis)):
    print("  ", f)
