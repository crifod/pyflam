"""Tuscany pyroconvection-TYPE forecast maps (Castellnou 2022 ladder) -- 3-hourly,
24 h -- from GFS and from ICON-2I, classified by pyflam.pyroconvection_type.

Reproduces the Catalan "Tipus de piroconvecció" panel: 8 maps (00..21Z) coloured by
pyroconvection prototype (surface plume -> convection plume -> overshooting pyroCu
-> resilient pyroCu -> deep pyroCu/pyroCb), output as an 8-panel PNG + one 2 km
GeoTIFF per timestep.

Atmosphere sources (source-agnostic classifier):
  * gfs    -- GFS 0.25 deg soundings (AWS Open-Data via Herbie), already cached.
  * icon2i -- ICON-2I 2.2 km convection-permitting (Italy), live via open-meteo
              (italia_meteo_arpae_icon_2i), no API key.

Classifier inputs per cell: LCL/ABL ratio, mixed-layer dtheta/dz (1000->850 hPa),
upper-layer gamma-theta (700->500 hPa). ABL from a theta-excess mixing height
(open-meteo ICON-2I exposes no PBL field). Pure-atmospheric potential type (no
fuel gate) -- matching the operational product.
"""
from __future__ import annotations

import sys
import warnings
from datetime import datetime, timedelta

import numpy as np
import requests

import pyflam
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))  # repo root
from pyflam_gui.core.pyroconv import surface_layer_terms, PYROCONV_NODATA
from pyflam.atmosphere import (
    lcl_height_m, theta_kelvin, pyroconvection_type,
    PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL, PYROCONVECTION_TYPE_COLOR,
)

warnings.simplefilter("ignore")

RUN = datetime(2026, 6, 26, 0)
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
CACHE = "/private/tmp/claude-501/-Users-cristianofoderi--softEST-firelab-flammap6-install-0828-2025-pyflam/a9d3d386-f4fb-49aa-aff5-8e25de82a0a4/scratchpad/gfs"
OUTDIR = "/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs"
RASTERDIR = OUTDIR + "/pyroconv_rasters"

# Tuscany + archipelago framing (matches the reference image).
LON0, LON1, LAT0, LAT1 = 9.6, 12.4, 42.2, 44.5
SAMPLE_DEG = 0.05          # ~5 km atmospheric sampling grid
OUT_DEG = 0.02            # ~2 km output raster
LEVELS_HPA = [1000, 925, 850, 700, 600, 500]
PLEV_BARO = np.array([44330.0 * (1.0 - (p / 1013.25) ** 0.1903) for p in LEVELS_HPA])


def grid(spacing):
    lons = np.arange(LON0, LON1 + 1e-9, spacing)
    lats = np.arange(LAT0, LAT1 + 1e-9, spacing)
    return lats, lons


def classify_grid(T2m, RH2m, Tlev, ps_hpa):
    """Per-cell pyroconvection type level (0..4). Arrays shaped [nlat,nlon];
    Tlev shaped [nlev,nlat,nlon] at LEVELS_HPA (C). Pure-atmospheric (no fuel gate)."""
    nlat, nlon = T2m.shape
    out = np.zeros((nlat, nlon), dtype=np.int16)
    z = PLEV_BARO
    th = theta_kelvin(Tlev, np.array(LEVELS_HPA)[:, None, None])
    i850, i700, i500 = LEVELS_HPA.index(850), LEVELS_HPA.index(700), LEVELS_HPA.index(500)
    # The 2 m state sits at the surface pressure, not at 1000 hPa (see
    # pyflam_gui.core.pyroconv.surface_layer_terms); ABL is measured above the ground.
    _, _, ml_grad, abl, usable = surface_layer_terms(
        T2m, {p: th[i] for i, p in enumerate(LEVELS_HPA)}, ps_hpa,
        dict(zip(LEVELS_HPA, z)))
    out[~usable] = PYROCONV_NODATA
    gamma = (th[i500] - th[i700]) / (z[i500] - z[i700])
    lcl = lcl_height_m(T2m, RH2m)
    abl = np.maximum(abl, 100.0)
    ratio = lcl / abl
    for a in range(nlat):
        for b in range(nlon):
            # shallow / nocturnal boundary layer -> no pyroconvective cloud: a
            # surface plume only (this is what renders the quiet white nights).
            if not usable[a, b]:
                continue
            if abl[a, b] < 600.0:
                out[a, b] = 0
                continue
            t = pyroconvection_type(lcl_abl_ratio=float(ratio[a, b]),
                                    ml_theta_gradient=float(ml_grad[a, b]),
                                    gamma_theta=float(gamma[a, b]))
            out[a, b] = PYROCONVECTION_TYPE_LEVEL[t]
    return out


# ---------------- ICON-2I via open-meteo --------------------------------------
def fetch_icon2i(lats, lons):
    """Return dict: 'T2m','RH2m' [h,nlat,nlon]; 'Tlev' [h,nlev,nlat,nlon]."""
    pts = [(la, lo) for la in lats for lo in lons]
    tvars = [f"temperature_{p}hPa" for p in LEVELS_HPA]
    hourly = ["temperature_2m", "relative_humidity_2m", "surface_pressure"] + tvars
    nlat, nlon, nh = len(lats), len(lons), len(HOURS)
    T2m = np.empty((nh, nlat, nlon)); RH2m = np.empty((nh, nlat, nlon))
    PS = np.empty((nh, nlat, nlon))
    Tlev = np.empty((nh, len(LEVELS_HPA), nlat, nlon))
    import time
    CH = 150
    flat = []
    nchunk = (len(pts) + CH - 1) // CH
    for i in range(0, len(pts), CH):
        chunk = pts[i:i+CH]
        la = ",".join(f"{a:.4f}" for a, _ in chunk)
        lo = ",".join(f"{b:.4f}" for _, b in chunk)
        url = ("https://api.open-meteo.com/v1/forecast?"
               f"latitude={la}&longitude={lo}&hourly={','.join(hourly)}"
               "&models=italia_meteo_arpae_icon_2i"
               "&start_date=2026-06-26&end_date=2026-06-26&timezone=UTC")
        sys.stderr.write(f"  ICON-2I chunk {i//CH+1}/{nchunk}\n")
        for attempt in range(5):                       # retry on rate-limit
            r = requests.get(url, timeout=90).json()
            if isinstance(r, dict) and r.get("error"):
                sys.stderr.write(f"    rate-limited ({r.get('reason','')[:50]}); backoff\n")
                time.sleep(8 * (attempt + 1))
                continue
            break
        if isinstance(r, dict):
            r = [r]
        if any("hourly" not in x for x in r):
            raise RuntimeError(f"open-meteo error chunk {i//CH+1}: {r[:1]}")
        flat.extend(r)
        time.sleep(2.0)                                # pace under the per-minute limit
    # map back to grid (pts order = lat-major)
    for idx, rec in enumerate(flat):
        a, b = divmod(idx, nlon)
        h = rec["hourly"]
        for hi, hour in enumerate(HOURS):
            T2m[hi, a, b] = h["temperature_2m"][hour]
            RH2m[hi, a, b] = h["relative_humidity_2m"][hour]
            PS[hi, a, b] = h["surface_pressure"][hour]                 # hPa
            for li, p in enumerate(LEVELS_HPA):
                Tlev[hi, li, a, b] = h[f"temperature_{p}hPa"][hour]
    return dict(T2m=T2m, RH2m=RH2m, Tlev=Tlev, PS=PS)


# ---------------- GFS via cached Herbie ---------------------------------------
def fetch_gfs(lats, lons):
    from herbie import Herbie
    import xarray as xr
    nlat, nlon, nh = len(lats), len(lons), len(HOURS)
    T2m = np.empty((nh, nlat, nlon)); RH2m = np.empty((nh, nlat, nlon))
    PS = np.empty((nh, nlat, nlon))
    Tlev = np.empty((nh, len(LEVELS_HPA), nlat, nlon))
    lon2d, lat2d = np.meshgrid(lons % 360.0, lats)
    lon_da = xr.DataArray(lon2d.ravel(), dims="p"); lat_da = xr.DataArray(lat2d.ravel(), dims="p")
    for hi, hour in enumerate(HOURS):
        sys.stderr.write(f"  GFS fxx={hour:02d}\n")
        prov = pyflam.atmosphere.fetch_gfs(run=RUN, fxx=hour, cache_dir=CACHE)
        sfc = prov.ds.sel(latitude=lat_da, longitude=lon_da, method="nearest")
        T2m[hi] = (np.asarray(sfc["t2m"].values, float).reshape(nlat, nlon) - 273.15)
        RH2m[hi] = np.asarray(sfc["r2"].values, float).reshape(nlat, nlon)
        # real pressure-level T at all of 1000/925/850/700/600/500 hPa
        h = Herbie(RUN, model="gfs", product="pgrb2.0p25", fxx=hour, save_dir=CACHE)
        sp = h.xarray(r":PRES:surface:", remove_grib=False)
        sp = sp[0] if isinstance(sp, list) else sp
        PS[hi] = np.asarray(sp.sel(latitude=lat_da, longitude=lon_da, method="nearest")
                            ["sp"].values, float).reshape(nlat, nlon) / 100.0   # Pa -> hPa
        ds = h.xarray(r":TMP:(1000|925|850|700|600|500) mb:", remove_grib=False)
        if isinstance(ds, list):
            m = ds[0]
            for d in ds[1:]:
                m = xr.merge([m, d], compat="override", combine_attrs="override")
            ds = m
        lev = np.asarray(ds["isobaricInhPa"].values, float)
        col = ds.sel(latitude=lat_da, longitude=lon_da, method="nearest")
        tg = np.asarray(col["t"].values, float).reshape(len(lev), nlat, nlon) - 273.15
        bylev = {int(p): tg[i] for i, p in enumerate(lev)}
        for li, p in enumerate(LEVELS_HPA):
            Tlev[hi, li] = bylev[p]
    return dict(T2m=T2m, RH2m=RH2m, Tlev=Tlev, PS=PS)


# ---------------- render ------------------------------------------------------
def render(cat, lats, lons, source, label):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    import rasterio
    from rasterio.transform import from_origin
    os.makedirs(RASTERDIR, exist_ok=True)
    colors = [PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES]
    cmap = ListedColormap(colors); norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    fig, axes = plt.subplots(1, len(HOURS), figsize=(2.0*len(HOURS), 3.0), constrained_layout=True)
    for hi, hour in enumerate(HOURS):
        ax = axes[hi]
        ax.set_facecolor("#bdbdbd")      # nodata (high ground / no surface p) = grey
        ax.imshow(np.ma.masked_less(cat[hi], 0), origin="lower", extent=[LON0, LON1, LAT0, LAT1],
                  cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_title(f"{RUN:%Y-%m-%d} {hour:02d}Z", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"Tipus de piroconvecció — {label}", fontsize=13)
    # legend
    from matplotlib.patches import Patch
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="gray",
                 label=t.replace("_", " ")) for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    png = f"{OUTDIR}/pyroconv_map_tuscany_{source}_2026-06-26.png"
    fig.savefig(png, dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 2 km GeoTIFFs (EPSG:4326), nearest-upsampled from the sampling grid
    olats, olons = grid(OUT_DEG)
    li = np.clip(np.searchsorted(lats, olats) , 0, len(lats)-1)
    bi = np.clip(np.searchsorted(lons, olons), 0, len(lons)-1)
    transform = from_origin(LON0, olats[-1] + OUT_DEG/2, OUT_DEG, OUT_DEG)
    for hi, hour in enumerate(HOURS):
        up = cat[hi][np.ix_(li, bi)][::-1]   # north-up
        path = f"{RASTERDIR}/pyroconv_{source}_{hour:02d}Z.tif"
        with rasterio.open(path, "w", driver="GTiff", height=up.shape[0],
                           width=up.shape[1], count=1, dtype="int16",
                           crs="EPSG:4326", transform=transform, nodata=-1) as d:
            d.write(up.astype("int16"), 1)
    return png


def main():
    sources = sys.argv[1:] or ["gfs", "icon2i"]
    lats, lons = grid(SAMPLE_DEG)
    sys.stderr.write(f"grid {len(lats)}x{len(lons)} @ {SAMPLE_DEG} deg\n")
    for src in sources:
        sys.stderr.write(f"=== {src} ===\n")
        data = fetch_gfs(lats, lons) if src == "gfs" else fetch_icon2i(lats, lons)
        cat = np.stack([classify_grid(data["T2m"][hi], data["RH2m"][hi], data["Tlev"][hi], data["PS"][hi])
                        for hi in range(len(HOURS))])
        label = ("GFS 0.25° (25 km soundings, output 2 km)" if src == "gfs"
                 else "ICON-2I 2.2 km (output 2 km)")
        png = render(cat, lats, lons, src, label)
        # quick category tally at 12Z
        names = {i: t for t, i in PYROCONVECTION_TYPE_LEVEL.items()}
        u, c = np.unique(cat[HOURS.index(12)], return_counts=True)
        tally = {names.get(int(k), 'nodata'): int(v) for k, v in zip(u, c)}
        print(f"{src}: wrote {png}")
        print(f"   12Z category counts: {tally}")


if __name__ == "__main__":
    main()
