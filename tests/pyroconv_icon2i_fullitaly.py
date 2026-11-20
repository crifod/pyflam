"""FULL-TUSCANY native 2.2 km pyroconvection map from the full-Italy ICON-2I GRIB.

Source: AgenziaItaliaMeteo / MISTRAL MeteoHub open data (CC-BY), full national 2.2 km
domain -- so it covers ALL of Tuscany (incl. Siena, Grosseto), unlike the Arpae
Emilia-Romagna cutout. Per-variable/per-level GRIB files, subset to the Tuscany bbox,
classified by pyflam.pyroconvection_type. Potential + fuel-gated panels + GeoTIFFs.
"""
from __future__ import annotations

import sys, warnings
import numpy as np
import xarray as xr
from pyproj import Transformer

import pyflam
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))  # repo root
from pyflam_gui.core.pyroconv import surface_layer_terms, PYROCONV_NODATA
from pyflam import units, fuel_models
from pyflam.atmosphere import (
    lcl_height_m, theta_kelvin, pyroconvection_type, equilibrium_moisture_content,
    relative_humidity_from_dewpoint, PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL,
    PYROCONVECTION_TYPE_COLOR,
)

warnings.simplefilter("ignore")

# Inputs/outputs overridable from the environment; defaults are the original June run.
import os as _os
RUN_DATE = _os.environ.get("RUN_DATE", "2026-06-26")
D = _os.environ.get("ICON2I_DIR", "/private/tmp/claude-501/-Users-cristianofoderi--softEST-firelab-flammap6-install-0828-2025-pyflam/a9d3d386-f4fb-49aa-aff5-8e25de82a0a4/scratchpad/icon2i_it")
OUTDIR = _os.environ.get("OUTDIR", "/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs")
RASTERDIR = OUTDIR + "/pyroconv_rasters_native_fulltuscany"
LCP = _os.environ.get("LCP", "/Users/cristianofoderi/DATI/FUEL_TOS/pyflam_canopy_tuscany/canopy_tuscany.lcp")
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
LON0, LON1, LAT0, LAT1 = 9.6, 12.5, 42.2, 44.6        # Tuscany bbox
LEVELS = [850, 700, 500]
# barometric heights for the 5 reference levels used in the classifier (1000 ~ T2m)
Z = {1000: 110.0, 850: 1457.0, 700: 3012.0, 500: 5574.0}
FLI_GATE_KW = 1.0e4
_TO3035 = Transformer.from_crs("EPSG:4326", "EPSG:3035", always_xy=True)


def _ds(fname):
    """Open a single-variable ICON-2I GRIB, subset to Tuscany; return (ds, varname)."""
    ds = xr.open_dataset(f"{D}/{fname}", engine="cfgrib", backend_kwargs={"indexpath": ""})
    la = ds["latitude"].values
    ds = (ds.sel(latitude=slice(LAT1, LAT0), longitude=slice(LON0, LON1)) if la[0] > la[-1]
          else ds.sel(latitude=slice(LAT0, LAT1), longitude=slice(LON0, LON1)))
    return ds, list(ds.data_vars)[0]


def load():
    sys.stderr.write("reading full-Italy ICON-2I GRIB (Tuscany subset)...\n")
    T = {}
    for p, fn in {850: "T850.grib", 700: "T700.grib", 500: "T500.grib"}.items():
        ds, v = _ds(fn); T[p] = ds[v].values - 273.15
    t2, v = _ds("T2M.grib"); T2m = t2[v].values - 273.15
    td, v = _ds("TD2M.grib"); Td2m = td[v].values - 273.15
    u, v = _ds("U10.grib"); U = u[v].values
    vv, vn = _ds("V10.grib"); V = vv[vn].values
    ps, vp = _ds("PS.grib"); PS = ps[vp].values / 100.0     # Pa -> hPa (ICON2I_FIRE_FIELDS)
    lat = t2["latitude"].values; lon = t2["longitude"].values
    steps_h = (t2["step"].values / np.timedelta64(1, "h")).astype(int)
    idx = [int(np.where(steps_h == h)[0][0]) for h in HOURS]
    return dict(lat=lat, lon=lon, idx=idx, T=T, T2m=T2m, Td2m=Td2m, U=U, V=V, PS=PS)


def classify(T2m, RH2m, Tl, ps_hpa, fli=None):
    nlat, nlon = T2m.shape
    out = np.zeros((nlat, nlon), np.int16)
    th = {p: theta_kelvin(Tl[p], float(p)) for p in (850, 700, 500)}
    # 2 m state at the surface pressure, ABL above the ground (surface_layer_terms)
    _, _, ml_grad, abl, usable = surface_layer_terms(T2m, th, ps_hpa, Z, abl_levels=(850, 700))
    out[~usable] = PYROCONV_NODATA
    gamma = (th[500] - th[700]) / (Z[500] - Z[700])
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
            out[a, b] = PYROCONVECTION_TYPE_LEVEL[pyroconvection_type(
                lcl_abl_ratio=float(ratio[a, b]), ml_theta_gradient=float(ml_grad[a, b]),
                gamma_theta=float(gamma[a, b]), fireline_intensity_kw=fkw,
                fli_threshold_kw=FLI_GATE_KW)]
    return out


def _known(n):
    try:
        fuel_models.get(int(n)); return True
    except KeyError:
        return False


def lcp_fields(lat, lon):
    ls = pyflam.Landscape.from_lcp(LCP)
    n = int(ls.fuel_model.max()) + 1
    burn = np.array([fuel_models.get(int(k)).is_burnable if _known(k) else False for k in range(n)])
    lon2d, lat2d = np.meshgrid(lon, lat)
    x, y = _TO3035.transform(lon2d.ravel(), lat2d.ravel())
    c = ((x - ls.west) / ls.cellsize_x).astype(int); r = ((ls.north - y) / ls.cellsize_y).astype(int)
    inside = (r >= 0) & (r < ls.shape[0]) & (c >= 0) & (c < ls.shape[1])
    rr, cc = np.clip(r, 0, ls.shape[0]-1), np.clip(c, 0, ls.shape[1]-1)
    fnum = np.full(r.shape, -1); fnum[inside] = ls.fuel_model[rr[inside], cc[inside]]
    sh = lat2d.shape
    return dict(fuel=fnum.reshape(sh), slope=ls.slope_tangent[rr, cc].reshape(sh),
                cbh=(ls.canopy_base_height[rr, cc]/10.0).reshape(sh),
                cbd=(ls.canopy_bulk_density[rr, cc]/100.0).reshape(sh),
                ch=(ls.canopy_height[rr, cc]/10.0).reshape(sh)), burn


def fli_grid(T2m, RH2m, wspd, lf, burn):
    nlat, nlon = T2m.shape; out = np.zeros((nlat, nlon))
    for a in range(nlat):
        for b in range(nlon):
            fn = int(lf["fuel"][a, b])
            if fn < 0 or fn >= burn.size or not burn[fn]:
                continue
            fuel = fuel_models.get(fn)
            emc = equilibrium_moisture_content(T2m[a, b], RH2m[a, b]) / 100.0
            u = max(float(wspd[a, b]), 0.0)
            surf = pyflam.spread(fuel, m_1h=emc, m_10h=emc+0.01, m_100h=emc+0.02,
                                 m_live_herb=0.70, m_live_woody=0.90,
                                 wind_midflame=units.m_per_s_to_ft_per_min(u)*0.4,
                                 slope=float(lf["slope"][a, b]))
            fli = units.btu_per_ft_s_to_kw_per_m(surf.fireline_intensity)
            if lf["cbh"][a, b] > 0 and lf["cbd"][a, b] > 0:
                cr = pyflam.crown_fire_behavior(surf, canopy_base_height=lf["cbh"][a, b],
                    canopy_bulk_density=lf["cbd"][a, b], foliar_moisture=100.0,
                    wind_20ft_ft_per_min=units.m_per_s_to_ft_per_min(u*0.87),
                    canopy_fuel_load=max(lf["cbd"][a, b]*max(lf["ch"][a, b]-lf["cbh"][a, b], 0), 0),
                    m_1h=emc, m_10h=emc+0.01, m_100h=emc+0.02, m_live_herb=0.70,
                    m_live_woody=0.90, crown_spread="cruz2005")
                fli = max(fli, cr.fireline_intensity)
            out[a, b] = fli
    return out


def render(cats, lat, lon, tag, label):
    import os, matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import rasterio; from rasterio.transform import from_origin
    os.makedirs(RASTERDIR, exist_ok=True)
    flip = lat[0] > lat[-1]; ext = [lon.min(), lon.max(), lat.min(), lat.max()]
    colors = [PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES]
    cmap = ListedColormap(colors); norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    fig, ax = plt.subplots(1, len(HOURS), figsize=(2.1*len(HOURS), 3.0), constrained_layout=True)
    for hi, hour in enumerate(HOURS):
        a = cats[hi][::-1] if flip else cats[hi]
        ax[hi].set_facecolor("#bdbdbd")      # nodata (high ground / no surface p) = grey
        ax[hi].imshow(np.ma.masked_less(a, 0), origin="lower", extent=ext, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax[hi].set_title(f"{RUN_DATE} {hour:02d}Z", fontsize=8); ax[hi].set_xticks([]); ax[hi].set_yticks([])
    fig.suptitle(f"Tipus de piroconveccio ({tag}) -- {label}", fontsize=12)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="gray", label=t.replace("_", " "))
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.03))
    png = f"{OUTDIR}/pyroconv_map_tuscany_icon2i_fulltuscany_{tag}_{RUN_DATE}.png"
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    dlon = float(abs(lon[1]-lon[0])); dlat = float(abs(lat[1]-lat[0]))
    transform = from_origin(lon.min()-dlon/2, lat.max()+dlat/2, dlon, dlat)
    for hi, hour in enumerate(HOURS):
        arr = cats[hi] if flip else cats[hi][::-1]
        with rasterio.open(f"{RASTERDIR}/pyroconv_icon2i_{tag}_{hour:02d}Z.tif", "w", driver="GTiff",
                           height=arr.shape[0], width=arr.shape[1], count=1, dtype="int16",
                           crs="EPSG:4326", transform=transform, nodata=-1) as d:
            d.write(arr.astype("int16"), 1)
    return png


def main():
    d = load(); lat, lon, idx = d["lat"], d["lon"], d["idx"]
    sys.stderr.write(f"Tuscany grid {len(lat)}x{len(lon)} native; lat {lat.min():.2f}-{lat.max():.2f} "
                     f"lon {lon.min():.2f}-{lon.max():.2f}\n")
    lf, burn = lcp_fields(lat, lon)
    RH = relative_humidity_from_dewpoint(d["T2m"], d["Td2m"])
    pot, gate = [], []
    for hi, si in enumerate(idx):
        Tl = {p: d["T"][p][si] for p in LEVELS}
        pot.append(classify(d["T2m"][si], RH[si], Tl, d["PS"][si]))
        wsp = np.hypot(d["U"][si], d["V"][si])
        gate.append(classify(d["T2m"][si], RH[si], Tl, d["PS"][si], fli=fli_grid(d["T2m"][si], RH[si], wsp, lf, burn)))
    pot = np.stack(pot); gate = np.stack(gate)
    lab = "ICON-2I 2.2 km NATIVE full Tuscany (MISTRAL/AgenziaItaliaMeteo GRIB)"
    render(pot, lat, lon, "potential", lab); png = render(gate, lat, lon, "gated", lab)
    names = {i: t for t, i in PYROCONVECTION_TYPE_LEVEL.items()}
    for nm, c in [("potential", pot), ("gated", gate)]:
        u, cnt = np.unique(c[HOURS.index(12)], return_counts=True)
        print(f"fulltuscany {nm} 12Z:", {names.get(int(k), 'nodata'): int(v) for k, v in zip(u, cnt)})
    print("wrote", png)


if __name__ == "__main__":
    main()
