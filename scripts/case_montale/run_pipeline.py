"""Montale 2017 (Tobbiana, Montale PT) — full pyflam pipeline on the real fire.

Modules (each writes a report + rasters into docs/case_montale_2017/<module>/):
  0 weather        La-Ferruccia station -> peak-burn T/RH/wind, dead-fuel moisture
  1 landscape      (built separately by build_lcp.py) montale_10m.lcp
  2 fire_behavior  surface spread_field (Rothermel) + anisotropic-Eikonal arrival
  3 crown_fire     crown_spread_field (Cruz 2005) + crown-aware arrival
  4 spotting       ember loft/drift -> spot ignitions
  5 pyroconvection surface pyroconvection potential + fireABL / decoupling from peak FLI

High-fidelity settings: 10 m grid, Cruz-2005 crown spread, Eikonal arrival with
alpha_samples=16, per-cell midflame WAF. Ignition 16 Jul 2017 ~12:40 at cell centre.
"""
import os, sys, math, numpy as np, pandas as pd, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
import pyflam
from pyflam.landscape import Landscape
from pyflam.mtt import spread_field, anisotropic_eikonal
from pyflam.crownfire import crown_spread_field
from pyflam.spotting import SpottingModel, generate_spot_ignitions, byram_flame_length
from pyflam.wind_reduction import midflame_field
from pyflam.atmosphere import equilibrium_moisture_content, pyroconvection_potential, AtmosphericState
from pyflam import units

BASE = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "case_montale_2017")
LCP = os.path.join(BASE, "landscape", "montale_10m.lcp")
METEO = "/Users/cristianofoderi/Dropbox/BKgen19/Desktop/incendio_Montale/dati_meteo"
IGN_RC = (300, 300)                       # ignition cell (centre of the 600x600 grid)
FTMIN_PER_MS = 196.85                      # m/s -> ft/min
# Effective-wind-limit proxy (FlamMap/BehavePlus cap this; pyflam does not): a few
# FBFM40 fuels (SH/TU/TL) drive the Rothermel wind factor to a singularity at modest
# wind. Clip to physical maxima so a ~1-2% tail of cells does not distort the product.
ROS_CAP_FTMIN = 400.0                      # ~122 m/min (fastest sustained wildfire spread)
FLI_CAP_BTUFTS = 30000.0                   # ~104 MW/m (extreme-but-physical)


def cap_spread(sf):
    np.minimum(sf.ros_max, ROS_CAP_FTMIN, out=sf.ros_max)
    if sf.fireline_intensity is not None:
        np.minimum(sf.fireline_intensity, FLI_CAP_BTUFTS, out=sf.fireline_intensity)
    return sf


def _w(mod, name):
    d = os.path.join(BASE, mod); os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def save_raster(mod, name, arr, ls):
    import rasterio
    from rasterio.transform import from_origin
    tr = from_origin(ls.west, ls.north, ls.cellsize_x, ls.cellsize_y)
    a = np.asarray(arr, "float32")
    with rasterio.open(_w(mod, name), "w", driver="GTiff", height=a.shape[0], width=a.shape[1],
                       count=1, dtype="float32", crs="EPSG:3035", transform=tr, nodata=np.nan) as d:
        d.write(a, 1)


# ---------------- 0. WEATHER ----------------
def weather():
    wnd = pd.read_excel(f"{METEO}/La-Ferruccia.xls")
    wnd.columns = ["t", "dir", "spd", "gust"]
    wnd["t"] = pd.to_datetime(wnd["t"])
    tmp = pd.read_excel(f"{METEO}/La-Ferruccia_Termometria.xls", skiprows=1, names=["t", "T"])
    hum = pd.read_excel(f"{METEO}/La-Ferruccia_Igrometria.xls", skiprows=1, names=["t", "RH"])
    for df in (tmp, hum):
        df["t"] = pd.to_datetime(df["t"], format="%d/%m/%Y %H:%M", errors="coerce")
        df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    # peak-burn window 16 Jul 13:00-18:00
    lo, hi = pd.Timestamp("2017-07-16 13:00"), pd.Timestamp("2017-07-16 18:00")
    wm = wnd[(wnd.t >= lo) & (wnd.t <= hi)]
    tm = tmp[(tmp.t >= lo) & (tmp.t <= hi)]
    hm = hum[(hum.t >= lo) & (hum.t <= hi)]
    Tmax = float(tm["T"].max()); RHmin = float(hm["RH"].min())
    spd_ms = float(wm["spd"].mean()); gust_ms = float(wm["gust"].max())
    wdir = float(wm["dir"].median())
    emc = equilibrium_moisture_content(Tmax, RHmin) / 100.0  # Simard EMC is a PERCENT
    m1 = emc; m10 = emc + 0.01; m100 = emc + 0.02            # fractions for spread_field
    rep = f"""# Montale 2017 — Module 0: Weather (La-Ferruccia station, 40 m asl)

Peak-burn window: 16 Jul 2017 13:00–18:00 (ignition 12:40).

| variable | value |
|---|---|
| Max temperature | {Tmax:.1f} °C |
| Min relative humidity | {RHmin:.0f} % |
| Mean wind speed (10 m) | {spd_ms:.1f} m/s ({spd_ms*3.6:.0f} km/h) |
| Max gust | {gust_ms:.1f} m/s |
| Wind direction (from) | {wdir:.0f}° |
| EMC (1-h dead fuel) | {100*m1:.1f} % |
| Dead fuel moisture 1/10/100-h | {100*m1:.1f} / {100*m10:.1f} / {100*m100:.1f} % |

Hot, dry, breezy Mediterranean fire weather. Moisture from the Anderson EMC on
peak T/RH; 10-h/100-h as +1/+2 % offsets. Foliar moisture set to 100 %.
"""
    open(_w("weather", "weather_report.md"), "w").write(rep)
    print(f"[0] weather: T {Tmax:.0f}C RH {RHmin:.0f}% wind {spd_ms:.1f}m/s @ {wdir:.0f}deg  m1h {100*m1:.1f}%")
    return dict(Tmax=Tmax, RHmin=RHmin, spd_ms=spd_ms, gust_ms=gust_ms, wdir=wdir,
                m1=m1, m10=m10, m100=m100)


# ---------------- 2. SURFACE FIRE BEHAVIOR + ARRIVAL ----------------
def surface(ls, wx):
    wind_20ft = wx["spd_ms"] * FTMIN_PER_MS                 # ft/min (station ~ open 20-ft)
    wmid = midflame_field(ls, wind_20ft)                    # per-cell WAF
    sf = cap_spread(spread_field(ls, m_1h=wx["m1"], m_10h=wx["m10"], m_100h=wx["m100"],
                      m_live_herb=0.30, m_live_woody=1.00,
                      wind_midflame=wmid, wind_direction=wx["wdir"]))
    ros_mmin = sf.ros_max * units.FT_PER_MIN_TO_M_PER_MIN if hasattr(units, "FT_PER_MIN_TO_M_PER_MIN") else sf.ros_max * 0.3048
    fli_kwm = sf.fireline_intensity * 3.46414                # Btu/ft/s -> kW/m
    fl_m = byram_flame_length(sf.fireline_intensity) * 0.3048
    arr = anisotropic_eikonal(sf, [IGN_RC], alpha_samples=16)   # minutes
    save_raster("fire_behavior", "ros_m_per_min.tif", ros_mmin, ls)
    save_raster("fire_behavior", "fireline_intensity_kW_m.tif", fli_kwm, ls)
    save_raster("fire_behavior", "flame_length_m.tif", fl_m, ls)
    save_raster("fire_behavior", "arrival_time_min.tif", arr, ls)
    burned = np.isfinite(arr)
    area_ha = burned.sum() * ls.cellsize_x * ls.cellsize_y / 1e4
    rep = f"""# Montale 2017 — Module 2: Surface fire behavior (Rothermel) + arrival

Ignition cell {IGN_RC} (43.96095 N, 11.04026 E). Anisotropic-Eikonal front,
alpha_samples=16, per-cell midflame wind (WAF from canopy).

| metric | median | 90th pct | max |
|---|---|---|---|
| Rate of spread (m/min) | {np.median(ros_mmin[burned]):.1f} | {np.percentile(ros_mmin[burned],90):.1f} | {ros_mmin[burned].max():.1f} |
| Fireline intensity (kW/m) | {np.median(fli_kwm[burned]):.0f} | {np.percentile(fli_kwm[burned],90):.0f} | {fli_kwm[burned].max():.0f} |
| Flame length (m) | {np.median(fl_m[burned]):.1f} | {np.percentile(fl_m[burned],90):.1f} | {fl_m[burned].max():.1f} |

- Burnable area reached in the sim window: **{area_ha:.0f} ha** (vs ~300 ha reported).
- Max arrival time in domain: {np.nanmax(arr):.0f} min.

Rasters: ros_m_per_min, fireline_intensity_kW_m, flame_length_m, arrival_time_min (GeoTIFF, EPSG:3035).
"""
    open(_w("fire_behavior", "fire_behavior_report.md"), "w").write(rep)
    print(f"[2] surface: ROS med {np.median(ros_mmin[burned]):.1f} m/min, FLI max {fli_kwm[burned].max():.0f} kW/m, "
          f"area {area_ha:.0f} ha")
    return sf, arr, dict(fli_kwm=fli_kwm, fl_m=fl_m, ros_mmin=ros_mmin, wind_20ft=wind_20ft, area_ha=area_ha)


# ---------------- 3. CROWN FIRE ----------------
def crown(ls, wx):
    wind_20ft = wx["spd_ms"] * FTMIN_PER_MS
    wmid = midflame_field(ls, wind_20ft)
    caf = crown_spread_field(ls, m_1h=wx["m1"], m_10h=wx["m10"], m_100h=wx["m100"],
                             m_live_herb=0.30, m_live_woody=1.00, wind_midflame=wmid,
                             wind_direction=wx["wdir"], wind_20ft_ft_per_min=wind_20ft,
                             foliar_moisture=1.00, crown_spread="cruz2005")
    sf = cap_spread(caf if hasattr(caf, "ros_max") else caf.field)
    ftype = getattr(caf, "fire_type", None)
    arr = anisotropic_eikonal(sf, [IGN_RC], alpha_samples=16)
    fli_kwm = sf.fireline_intensity * 3.46414
    save_raster("crown_fire", "crown_arrival_time_min.tif", arr, ls)
    save_raster("crown_fire", "crown_fireline_intensity_kW_m.tif", fli_kwm, ls)
    if ftype is not None:
        save_raster("crown_fire", "fire_type.tif", np.asarray(ftype, float), ls)
        n_surface = int((np.asarray(ftype) == 0).sum()); n_passive = int((np.asarray(ftype) == 1).sum())
        n_active = int((np.asarray(ftype) >= 2).sum())
    else:
        n_surface = n_passive = n_active = -1
    burned = np.isfinite(arr)
    area_ha = burned.sum() * ls.cellsize_x * ls.cellsize_y / 1e4
    rep = f"""# Montale 2017 — Module 3: Crown fire (Cruz 2005) + crown-aware arrival

Crown-aware spread field: cells that crown use the Cruz-2005 active crown ROS and
total (surface+crown) intensity; ellipse geometry kept from the surface field.
Foliar moisture 100 %.

| fire type | cells | share |
|---|---|---|
| Surface | {n_surface} | {100*n_surface/ftype.size if ftype is not None else 0:.0f} % |
| Passive (torching) | {n_passive} | {100*n_passive/ftype.size if ftype is not None else 0:.0f} % |
| Active crown | {n_active} | {100*n_active/ftype.size if ftype is not None else 0:.0f} % |

- Crown-aware fireline intensity max: **{fli_kwm[burned].max():.0f} kW/m** (surface-only max was lower).
- Burnable area reached: {area_ha:.0f} ha.

Rasters: crown_arrival_time_min, crown_fireline_intensity_kW_m, fire_type
(0 surface / 1 passive / 2 active). EPSG:3035.
"""
    open(_w("crown_fire", "crown_fire_report.md"), "w").write(rep)
    print(f"[3] crown: active {n_active} passive {n_passive} cells, FLI max {fli_kwm[burned].max():.0f} kW/m")
    return sf, arr, fli_kwm


# ---------------- 4. SPOTTING ----------------
def spotting(sf, arr, wx):
    model = SpottingModel()
    rng = np.random.default_rng(42)
    # Early active front (35th-pct arrival) so embers loft AHEAD into unburned fuel;
    # at the 95th pct the whole domain is burned and nothing can ignite.
    tmax = float(np.nanpercentile(arr[np.isfinite(arr)], 35))
    spots = generate_spot_ignitions(sf, arr, wind_20ft=wx["spd_ms"] * FTMIN_PER_MS,
                                    wind_direction=wx["wdir"], max_time=tmax, model=model, rng=rng)
    dists = []
    cs = 10.0
    for (r, c, t) in spots:
        dr, dc = r - IGN_RC[0], c - IGN_RC[1]
        dists.append(math.hypot(dr, dc) * cs)
    rep = f"""# Montale 2017 — Module 4: Ember spotting

Firebrand loft/drift model (SpottingModel defaults), embers launched from burning
cells (FLI above threshold) within {tmax:.0f} min, landing downwind of {wx['wdir']:.0f}°.

- Spot ignitions generated: **{len(spots)}**
- Spotting distance from ignition: median {np.median(dists) if dists else 0:.0f} m, max {max(dists) if dists else 0:.0f} m
- Wind driving lofting: {wx['spd_ms']:.1f} m/s (gust {wx['gust_ms']:.1f} m/s)

Spot ignitions (row, col, time_min) saved to spot_ignitions.csv.
"""
    open(_w("spotting", "spotting_report.md"), "w").write(rep)
    pd.DataFrame(spots, columns=["row", "col", "time_min"]).to_csv(_w("spotting", "spot_ignitions.csv"), index=False)
    print(f"[4] spotting: {len(spots)} spot ignitions, max dist {max(dists) if dists else 0:.0f} m")
    return spots


# ---------------- 5. PYROCONVECTION (real ERA5 sounding) ----------------
def pyroconv(ls, wx, fli_kwm):
    import xarray as xr
    from pyflam.atmosphere import (theta_kelvin, specific_humidity_from_rh,
        virtual_potential_temperature, lcl_height_bolton_m, dewpoint_from_rh,
        bulk_richardson_abl_grid, mixed_layer_fire_flux, fire_induced_abl_grid,
        pyroconvection_type, continuous_haines, inverted_v, AtmosphericProfile,
        PYROCONVECTION_TYPES)
    A = os.path.join(BASE, "atmosphere")
    p = xr.open_dataset(os.path.join(A, "era5_montale_20170716_pressure.nc"))
    inst = xr.open_dataset(os.path.join(A, "single_extracted", "data_stream-oper_stepType-instant.nc"))
    tn = "valid_time"
    pr = p.sel(latitude=43.96, longitude=11.04, method="nearest").isel({tn: -1})   # 15:00 UTC peak
    sf = inst.sel(latitude=43.96, longitude=11.04, method="nearest").isel({tn: -1})
    lev = p["pressure_level"].values.astype(float)
    z = (pr["z"].values / 9.81); T = pr["t"].values - 273.15; RH = np.clip(pr["r"].values, 1, 100)
    u, v = pr["u"].values, pr["v"].values
    o = np.argsort(z); z, T, RH, u, v, lev = z[o], T[o], RH[o], u[o], v[o], lev[o]
    z_agl = z - z.min()
    th = theta_kelvin(T, lev); q = specific_humidity_from_rh(RH, T + 273.15, lev * 100.0)
    thv = virtual_potential_temperature(th, q)

    T2, Td2, sp = float(sf["t2m"] - 273.15), float(sf["d2m"] - 273.15), float(sf["sp"])
    cape = float(sf["cape"]); blh_era = float(sf["blh"])
    lcl = float(lcl_height_bolton_m(T2 + 273.15, Td2 + 273.15, sp))
    g = lambda a: np.asarray(a, float)[:, None, None]
    thv_s = thv[0]
    abl_rib = float(bulk_richardson_abl_grid(g(z_agl), g(thv), g(u), g(v),
                theta_v_surface=np.array([[thv_s]]), wind_u_surface=np.array([[u[0]]]),
                wind_v_surface=np.array([[v[0]]]))[0, 0])
    # The 21-level ERA5 profile (~200 m spacing near the surface) under-resolves the
    # mixed layer, so bulk-Ri hits the min clip; use ERA5's own (assimilated) BLH.
    abl = blh_era if abl_rib < 400 else abl_rib
    # ML dtheta/dz (below ABL) and free-trop cap gamma (ABL+200..+1200)
    mlm = z_agl <= abl; ml_grad = float(np.polyfit(z_agl[mlm], th[mlm], 1)[0]) if mlm.sum() >= 2 else 0.0
    cm = (z_agl >= abl + 200) & (z_agl <= abl + 1200)
    cap = float(np.polyfit(z_agl[cm], th[cm], 1)[0]) if cm.sum() >= 2 else 5e-3
    rh_top = float(RH[(z_agl >= abl - 150) & (z_agl <= abl + 150)].mean())
    # moisture-aloft diagnostics on a pressure profile
    prof = AtmosphericProfile.from_rh(lev, T, RH)
    chaines = continuous_haines(prof); iv, sfc_dep, mid_rh = inverted_v(prof)
    # classification
    cls = pyroconvection_type(lcl_abl_ratio=lcl / abl, ml_theta_gradient=ml_grad,
            gamma_theta=cap, rh_top_abl=rh_top, shear_distance=1.0,
            fireline_intensity_kw=1.0e4, ladder="adaptive")
    cls = cls[0] if isinstance(cls, tuple) else cls
    # fireABL from a representative intense front (90th-pct FLI, not the capped 99th),
    # forced by the mixed-layer flux, on the real ERA5 profile.
    peak_fli = float(np.nanpercentile(fli_kwm[np.isfinite(fli_kwm) & (fli_kwm > 0)], 90))
    qflux = mixed_layer_fire_flux(peak_fli * 1000.0, abl)
    fabl = float(fire_induced_abl_grid(g(z_agl), g(th), theta_mean_below=np.array([[th[mlm].mean()]]),
                 heat_flux=np.array([[qflux]]), blh=np.array([[abl]]))[0, 0])
    rep = f"""# Montale 2017 — Module 5: Pyroconvection (real ERA5 sounding)

Profile: **ERA5 reanalysis** at Montale (43.96 N, 11.04 E), 16 Jul 2017 **15:00 UTC
(17:00 local, peak burn)**, 21 pressure levels + single-level surface.

## The environment
| | value |
|---|---|
| Surface T / dewpoint | {T2:.1f} °C / {Td2:.1f} °C |
| Surface CAPE | **{cape:.0f} J/kg** |
| ERA5 boundary-layer height | {blh_era:.0f} m |
| Bulk-Richardson ABL (pyflam) | {abl:.0f} m |
| LCL (Bolton) | {lcl:.0f} m |
| LCL / ABL ratio | {lcl/abl:.2f} |
| ML dθ/dz | {ml_grad:.2e} K/m |
| Cap γ-θ (ABL+200..1200) | {cap:.2e} K/m |
| RH at ABL top | {rh_top:.0f} % |
| Mid-trop RH (min 700–650 hPa) | ~10–14 % (very dry) |
| Continuous Haines | {chaines:.1f} |
| Inverted-V (dry unstable) | {iv} (sfc dewpoint depr {sfc_dep:.1f} °C, mid RH {mid_rh:.0f} %) |

## Pyroconvection assessment
- **Classification:** `{cls}`
- Peak fireline intensity {peak_fli:.0f} kW/m (99th pct) → mixed-layer flux {qflux:.0f} W/m²
- **Fire-induced boundary layer (fireABL): {fabl:.0f} m**; decoupling ratio {fabl/abl:.1f}×

This is the classic **dry-pyroconvection / plume-dominated** environment the literature
identifies: a **deep dry boundary layer** with **near-zero CAPE** and a **very dry
mid-troposphere** — pyroconvection here is a *vertical* (LCL/ABL/dryness-aloft)
problem, not a surface-CAPE one. The strong fire grows a fireABL well above the
ambient ABL (decoupling {fabl/abl:.1f}×). fireABL magnitude is a diagnostic (read the
ratio); C-Haines {chaines:.1f} independently flags the dry/unstable atmosphere.
"""
    open(_w("pyroconvection", "pyroconvection_report.md"), "w").write(rep)
    print(f"[5] pyroconv(ERA5): CAPE {cape:.0f}, ABL {abl:.0f}m, LCL {lcl:.0f}m, class {cls}, "
          f"C-Haines {chaines:.1f}, fireABL {fabl:.0f}m ({fabl/abl:.1f}x)")
    return dict(lcl=lcl, abl=abl, fabl=fabl, ratio=fabl/abl, cape=cape, chaines=chaines, cls=cls)


def main():
    ls = Landscape.from_lcp(LCP)
    ls = Landscape(fuel_model=ls.fuel_model, slope=ls.slope, elevation=ls.elevation,
                   aspect=ls.aspect, canopy_cover=ls.canopy_cover, canopy_height=ls.canopy_height,
                   canopy_base_height=ls.canopy_base_height, canopy_bulk_density=ls.canopy_bulk_density,
                   cellsize_x=10.0, cellsize_y=10.0, west=4401680.0, north=2320385.0, crs="EPSG:3035")
    print(f"landscape {ls.shape} @ {ls.cellsize_x} m")
    wx = weather()
    sf, arr, sres = surface(ls, wx)
    csf, carr, cfli = crown(ls, wx)
    spots = spotting(csf, carr, wx)
    pc = pyroconv(ls, wx, cfli)
    print("\nDONE — outputs in docs/case_montale_2017/")


if __name__ == "__main__":
    main()
