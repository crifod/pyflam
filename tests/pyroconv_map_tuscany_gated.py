"""Tuscany pyroconvection-TYPE maps v2 -- finer ICON-2I sampling + FUEL GATE.

Adds to v1:
 (b) ICON-2I sampled finer (~4 km; the open-meteo free tier minutely=600/hourly~5000
     caps full 2.2 km gridding -- true native needs the MISTRAL ICON-2I GRIB, see
     icon2i_from_grib()).
 (a) FUEL GATE: per cell the realized fireline intensity is computed from the real
     Tuscany .lcp (fuel model + canopy) with the cell's GFS/ICON fuel moisture and
     wind; a cell only keeps a pyroCu category if its FLI >= 10 MW/m (Tedim 2018),
     else it is a surface plume. Sea / non-Tuscany / non-burnable -> surface plume.

Outputs per source: a POTENTIAL panel, a GATED panel, and gated 2 km GeoTIFFs.
"""
from __future__ import annotations

import sys, time, warnings
from datetime import datetime, timedelta
from collections import Counter

import numpy as np
import requests
from pyproj import Transformer

import pyflam
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))  # repo root
from pyflam_gui.core.pyroconv import surface_layer_terms, PYROCONV_NODATA
from pyflam import units, fuel_models
from pyflam.atmosphere import (
    lcl_height_m, theta_kelvin, pyroconvection_type, equilibrium_moisture_content,
    PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL, PYROCONVECTION_TYPE_COLOR,
)

warnings.simplefilter("ignore")

RUN = datetime(2026, 6, 26, 0)
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
CACHE = "/private/tmp/claude-501/-Users-cristianofoderi--softEST-firelab-flammap6-install-0828-2025-pyflam/a9d3d386-f4fb-49aa-aff5-8e25de82a0a4/scratchpad/gfs"
OUTDIR = "/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs"
RASTERDIR = OUTDIR + "/pyroconv_rasters_gated"
LCP = "/Users/cristianofoderi/DATI/FUEL_TOS/pyflam_canopy_tuscany/canopy_tuscany.lcp"

LON0, LON1, LAT0, LAT1 = 9.6, 12.4, 42.2, 44.5
OUT_DEG = 0.02
SAMPLE = {"gfs": 0.05, "icon2i": 0.07}      # icon2i kept budget-safe (open-meteo free hourly cap)
LEVELS_HPA = [1000, 925, 850, 700, 600, 500]
PLEV_BARO = np.array([44330.0 * (1.0 - (p / 1013.25) ** 0.1903) for p in LEVELS_HPA])
FLI_GATE_KW = 1.0e4
_TO3035 = Transformer.from_crs("EPSG:4326", "EPSG:3035", always_xy=True)


def grid(spacing):
    return np.arange(LAT0, LAT1 + 1e-9, spacing), np.arange(LON0, LON1 + 1e-9, spacing)


# ---------------- .lcp fuel sampler (for the gate) ----------------------------
def _is_forest(n):
    return (141 <= n <= 149) or (161 <= n <= 165) or (181 <= n <= 189) or (201 <= n <= 204)


class Lcp:
    def __init__(self, path):
        self.ls = pyflam.Landscape.from_lcp(path)
        n = int(self.ls.fuel_model.max()) + 1
        self.burn = np.array([fuel_models.get(int(k)).is_burnable
                              if _known(k) else False for k in range(n)])

    def fields(self, lats, lons):
        """Per map-cell fuel_num, slope_tan, cbh, cbd, ch (m); -1 fuel outside extent."""
        ls = self.ls
        nlat, nlon = len(lats), len(lons)
        lon2d, lat2d = np.meshgrid(lons, lats)
        x, y = _TO3035.transform(lon2d.ravel(), lat2d.ravel())
        c = ((x - ls.west) / ls.cellsize_x).astype(int)
        r = ((ls.north - y) / ls.cellsize_y).astype(int)
        inside = (r >= 0) & (r < ls.shape[0]) & (c >= 0) & (c < ls.shape[1])
        fnum = np.full(r.shape, -1)
        rr, cc = np.clip(r, 0, ls.shape[0]-1), np.clip(c, 0, ls.shape[1]-1)
        fnum[inside] = ls.fuel_model[rr[inside], cc[inside]]
        slope = ls.slope_tangent[rr, cc]
        cbh = ls.canopy_base_height[rr, cc] / 10.0
        cbd = ls.canopy_bulk_density[rr, cc] / 100.0
        ch = ls.canopy_height[rr, cc] / 10.0
        return dict(fuel=fnum.reshape(nlat, nlon), slope=slope.reshape(nlat, nlon),
                    cbh=cbh.reshape(nlat, nlon), cbd=cbd.reshape(nlat, nlon),
                    ch=ch.reshape(nlat, nlon))


def _known(n):
    try:
        fuel_models.get(int(n)); return True
    except KeyError:
        return False


def fli_grid(T2m, RH2m, wspd, lf, burn_lookup):
    """Per-cell realized fireline intensity (kW/m) from .lcp fuels + cell weather."""
    nlat, nlon = T2m.shape
    out = np.zeros((nlat, nlon))
    for a in range(nlat):
        for b in range(nlon):
            fn = int(lf["fuel"][a, b])
            if fn < 0 or fn >= burn_lookup.size or not burn_lookup[fn]:
                continue
            fuel = fuel_models.get(fn)
            emc = equilibrium_moisture_content(T2m[a, b], RH2m[a, b]) / 100.0
            m1, m10, m100 = emc, emc + 0.01, emc + 0.02
            u = max(float(wspd[a, b]), 0.0)
            wmf = units.m_per_s_to_ft_per_min(u) * 0.4
            surf = pyflam.spread(fuel, m_1h=m1, m_10h=m10, m_100h=m100,
                                 m_live_herb=0.70, m_live_woody=0.90,
                                 wind_midflame=wmf, slope=float(lf["slope"][a, b]))
            fli = units.btu_per_ft_s_to_kw_per_m(surf.fireline_intensity)
            cbh, cbd, ch = lf["cbh"][a, b], lf["cbd"][a, b], lf["ch"][a, b]
            if cbh > 0 and cbd > 0:
                cr = pyflam.crown_fire_behavior(
                    surf, canopy_base_height=cbh, canopy_bulk_density=cbd,
                    foliar_moisture=100.0, wind_20ft_ft_per_min=units.m_per_s_to_ft_per_min(u*0.87),
                    canopy_fuel_load=max(cbd*max(ch-cbh, 0.0), 0.0), m_1h=m1, m_10h=m10,
                    m_100h=m100, m_live_herb=0.70, m_live_woody=0.90, crown_spread="cruz2005")
                fli = max(fli, cr.fireline_intensity)
            out[a, b] = fli
    return out


# ---------------- classifier over the grid ------------------------------------
def classify_grid(T2m, RH2m, Tlev, ps_hpa, fli=None):
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
            if not usable[a, b]:
                continue
            if abl[a, b] < 600.0:
                out[a, b] = 0; continue
            fkw = None if fli is None else float(fli[a, b])
            t = pyroconvection_type(lcl_abl_ratio=float(ratio[a, b]),
                                    ml_theta_gradient=float(ml_grad[a, b]),
                                    gamma_theta=float(gamma[a, b]),
                                    fireline_intensity_kw=fkw, fli_threshold_kw=FLI_GATE_KW)
            out[a, b] = PYROCONVECTION_TYPE_LEVEL[t]
    return out


# ---------------- data sources ------------------------------------------------
def fetch_icon2i(lats, lons):
    pts = [(la, lo) for la in lats for lo in lons]
    tvars = [f"temperature_{p}hPa" for p in LEVELS_HPA]
    hourly = ["temperature_2m", "relative_humidity_2m", "wind_speed_10m", "surface_pressure"] + tvars
    nlat, nlon, nh = len(lats), len(lons), len(HOURS)
    T2m = np.empty((nh, nlat, nlon)); RH2m = np.empty((nh, nlat, nlon)); W = np.empty((nh, nlat, nlon))
    PS = np.empty((nh, nlat, nlon))
    Tlev = np.empty((nh, len(LEVELS_HPA), nlat, nlon))
    CH = 150; flat = []; nchunk = (len(pts)+CH-1)//CH
    for i in range(0, len(pts), CH):
        chunk = pts[i:i+CH]
        la = ",".join(f"{a:.4f}" for a, _ in chunk); lo = ",".join(f"{b:.4f}" for _, b in chunk)
        url = ("https://api.open-meteo.com/v1/forecast?"
               f"latitude={la}&longitude={lo}&hourly={','.join(hourly)}"
               "&models=italia_meteo_arpae_icon_2i&start_date=2026-06-26&end_date=2026-06-26&timezone=UTC")
        sys.stderr.write(f"  ICON-2I chunk {i//CH+1}/{nchunk}\n")
        for attempt in range(6):
            r = requests.get(url, timeout=90).json()
            if isinstance(r, dict) and r.get("error"):
                sys.stderr.write(f"    {r.get('reason','')[:40]}; backoff\n"); time.sleep(20*(attempt+1)); continue
            break
        if isinstance(r, dict):
            r = [r]
        if any("hourly" not in x for x in r):
            raise RuntimeError(f"open-meteo error chunk {i//CH+1}")
        flat.extend(r)
        time.sleep(max(len(chunk)/10.0, 1.0))           # pace <=600 pts/min
    for idx, rec in enumerate(flat):
        a, b = divmod(idx, nlon); h = rec["hourly"]
        for hi, hour in enumerate(HOURS):
            T2m[hi, a, b] = h["temperature_2m"][hour]; RH2m[hi, a, b] = h["relative_humidity_2m"][hour]
            W[hi, a, b] = (h["wind_speed_10m"][hour] or 0.0) / 3.6   # km/h -> m/s
            PS[hi, a, b] = h["surface_pressure"][hour]                 # hPa
            for li, p in enumerate(LEVELS_HPA):
                Tlev[hi, li, a, b] = h[f"temperature_{p}hPa"][hour]
    return dict(T2m=T2m, RH2m=RH2m, W=W, Tlev=Tlev, PS=PS)


def fetch_gfs(lats, lons):
    from herbie import Herbie
    import xarray as xr
    nlat, nlon, nh = len(lats), len(lons), len(HOURS)
    T2m = np.empty((nh, nlat, nlon)); RH2m = np.empty((nh, nlat, nlon)); W = np.empty((nh, nlat, nlon))
    PS = np.empty((nh, nlat, nlon))
    Tlev = np.empty((nh, len(LEVELS_HPA), nlat, nlon))
    lon2d, lat2d = np.meshgrid(lons % 360.0, lats)
    lon_da = xr.DataArray(lon2d.ravel(), dims="p"); lat_da = xr.DataArray(lat2d.ravel(), dims="p")
    for hi, hour in enumerate(HOURS):
        sys.stderr.write(f"  GFS fxx={hour:02d}\n")
        prov = pyflam.atmosphere.fetch_gfs(run=RUN, fxx=hour, cache_dir=CACHE)
        sfc = prov.ds.sel(latitude=lat_da, longitude=lon_da, method="nearest")
        T2m[hi] = np.asarray(sfc["t2m"].values, float).reshape(nlat, nlon) - 273.15
        RH2m[hi] = np.asarray(sfc["r2"].values, float).reshape(nlat, nlon)
        u = np.asarray(sfc["u10"].values, float).reshape(nlat, nlon)
        v = np.asarray(sfc["v10"].values, float).reshape(nlat, nlon)
        W[hi] = np.hypot(u, v)
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
    return dict(T2m=T2m, RH2m=RH2m, W=W, Tlev=Tlev, PS=PS)


# ---------------- render ------------------------------------------------------
def render(cat, lats, lons, source, label, tag):
    import os, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
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
        ax.set_title(f"{RUN:%Y-%m-%d} {hour:02d}Z", fontsize=8); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"Tipus de piroconveccio ({tag}) -- {label}", fontsize=13)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="gray", label=t.replace("_", " "))
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
    png = f"{OUTDIR}/pyroconv_map_tuscany_{source}_{tag}_2026-06-26.png"
    fig.savefig(png, dpi=130, bbox_inches="tight"); plt.close(fig)
    if tag == "gated":
        olats, olons = grid(OUT_DEG)
        li = np.clip(np.searchsorted(lats, olats), 0, len(lats)-1)
        bi = np.clip(np.searchsorted(lons, olons), 0, len(lons)-1)
        transform = from_origin(LON0, olats[-1] + OUT_DEG/2, OUT_DEG, OUT_DEG)
        for hi, hour in enumerate(HOURS):
            up = cat[hi][np.ix_(li, bi)][::-1]
            with rasterio.open(f"{RASTERDIR}/pyroconv_{source}_gated_{hour:02d}Z.tif", "w",
                               driver="GTiff", height=up.shape[0], width=up.shape[1], count=1,
                               dtype="int16", crs="EPSG:4326", transform=transform, nodata=-1) as d:
                d.write(up.astype("int16"), 1)
    return png


def main():
    sources = sys.argv[1:] or ["gfs", "icon2i"]
    lcp = Lcp(LCP)
    for src in sources:
        sys.stderr.write(f"=== {src} (sample {SAMPLE[src]} deg) ===\n")
        lats, lons = grid(SAMPLE[src])
        data = fetch_gfs(lats, lons) if src == "gfs" else fetch_icon2i(lats, lons)
        lf = lcp.fields(lats, lons)
        cat_pot = np.stack([classify_grid(data["T2m"][hi], data["RH2m"][hi], data["Tlev"][hi], data["PS"][hi])
                            for hi in range(len(HOURS))])
        cat_gate = np.stack([classify_grid(data["T2m"][hi], data["RH2m"][hi], data["Tlev"][hi], data["PS"][hi],
                             fli=fli_grid(data["T2m"][hi], data["RH2m"][hi], data["W"][hi], lf, lcp.burn))
                             for hi in range(len(HOURS))])
        res = "2.2 km model (sampled ~8 km)" if src == "icon2i" else "0.25deg (output 2 km)"
        render(cat_pot, lats, lons, src, f"{src.upper()} {res}", "potential")
        png = render(cat_gate, lats, lons, src, f"{src.upper()} {res}", "gated")
        names = {i: t for t, i in PYROCONVECTION_TYPE_LEVEL.items()}
        u, c = np.unique(cat_gate[HOURS.index(12)], return_counts=True)
        print(f"{src}: gated PNG {png}")
        print(f"   12Z gated counts: {dict((names.get(int(k), 'nodata'), int(v)) for k, v in zip(u, c))}")


# native 2.2 km path (free-tier API cannot grid all 16k cells): point pyflam at a
# downloaded MISTRAL ICON-2I GRIB and use it as a GriddedAtmosphere instead.
def icon2i_from_grib(path):
    var_map = {"temperature_K": "t2m", "relative_humidity": "r2",
               "wind_u": "u10", "wind_v": "v10"}
    return pyflam.atmosphere.open_atmosphere(path, source="gfs", var_map=var_map)


if __name__ == "__main__":
    main()
