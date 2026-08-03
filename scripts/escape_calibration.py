"""Calibrate the escape rung's prefactor against MISR plume tops -- Briggs' form, MISR's level.

The escape rung prices penetration of the ABL cap by inverting Briggs (1969):

    FP_req = C * (z_abl/2.6)^3 * U * s * pi rho cp T / g          [atmosphere._escape_cost_gw]

`C` is 1 for Briggs' own geometry. Expanding Tory & Kepert eq 25 on the same variables gives
the same law with an effective C ~ 0.34, and neither paper adjudicates. Sofiev et al. (2012)
fitted an injection-height formula to MISR that *could* settle the level, but it carries no
wind term and its stability term is calibrated only for |delta N^2/N0^2| < 0.1 -- it cannot
replace the rung. It can, however, supply the data to calibrate it.

The fit factorises, which is the point of doing it this way:

    predicted escape  <=>  FP_obs >= C * X        X = everything but C, per plume
                      <=>  r >= C                 r = FP_obs / X

So `r` is a *score* and `C` is just an operating point on its ROC curve. That splits the two
questions that are usually confounded:

  * Does the functional form have skill at all?   -- AUC of `r`, independent of C.
  * What is the right level?                      -- the operating point, chosen on the fitted C.

A form with no skill cannot be rescued by any C, and we would rather learn that than tune a
constant into a model that does not order the data.

Stages (each cached, run them in order):

    python scripts/escape_calibration.py --select      # pick the slice, list the ERA5 hours
    python scripts/escape_calibration.py --era5        # fetch reanalysis (slow: CDS queue)
    python scripts/escape_calibration.py --fit         # score, ROC, operating point

`--fit` accepts `--f-rad` (the FRP -> total firepower fraction). This is the largest single
uncertainty in the exercise: Sofiev sec. 4 asserts P_f ~ FRP outright, while the fire-radiative
literature puts the radiative fraction near 0.1-0.17. It is a declared parameter, never a
buried constant, and the fitted C absorbs it linearly -- so C and f_rad are not separately
identifiable from this data alone. What *is* identifiable is their ratio, and the AUC, which is
invariant to both.
"""
import argparse, collections, csv, datetime as dt, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
DATA = os.path.abspath(os.path.join(HERE, "..", "data", "misr"))

from pyflam.atmosphere import (theta_kelvin, _G, _PI, _CP_DRY, _RD, _BRIGGS_C,
                               _BRIGGS_MIN_WIND_MS, _ESCAPE_EZ_FRAC,
                               _ESCAPE_EZ_MIN_M)

LEVELS = [1000, 950, 900, 850, 800, 750, 700, 650]
GRID = 0.5

# Model levels. The pressure-level run could not measure the entrainment zone: its 50 hPa
# spacing is ~500 m below 4 km, against a zone 270 m deep, so 46 % of columns had *no* level
# inside the layer the stability is read over and `s` came out of an interpolation. L137 levels
# 100-137 span 10 m to ~4.3 km at 20-230 m spacing, which resolves it.
ML_LEVELS = list(range(100, 138))
ML_PARAMS = {"t": 130, "q": 133, "u": 131, "v": 132}
AB_TABLE = os.path.abspath(os.path.join(HERE, "..", "data", "era5_l137_ab.csv"))
_RV_OVER_RD_M1 = 0.609133        # Rv/Rd - 1, for virtual temperature

# MODIS IGBP land-cover classes, as MERLIN reports them per plume.
BIOME_GROUPS = {
    "forest": {1, 2, 3, 4, 5},          # evergreen/deciduous needleleaf+broadleaf, mixed
    "grassland": {10},
    "agriculture": {12, 14},            # croplands, cropland/natural-vegetation mosaic
    "savanna": {8, 9},                  # what dominated the first (Africa) slice
    "shrubland": {6, 7},
}

# Named slices. The first Africa run is kept so its result stays reproducible; `boreal` is the
# regime-corrected rerun.
#
# NOTE ON "AFTERNOON": it is not available. MISR flies on Terra, whose sun-synchronous orbit
# fixes the descending node near 10:30 local solar time, and the catalogue confirms it -- LST
# runs 10.1-11.8 with a median of 10.56, and 0.3 % of plumes are past noon. The latest LST the
# instrument can offer comes from high latitudes, where the swath geometry pushes it to
# ~11.3-12.0; that is what `boreal` selects, and it coincides with the forested, high-intensity
# regime. A genuine afternoon sample needs a different instrument, not a different query.
# Per-region high-firepower slices, all biomes. Going global does not widen the marginal
# plume-top spread (global high-FRP IQR 1184 m against the boreal slice's 1125 m, x1.05) --
# what it adds is *explainable between-regime* variance: regional medians run 1305 m (Africa)
# to 2160 m (North America), which a physical model should order and a constant cannot. It
# also supplies the structural variety in the theta profile that the integral and differential
# formulations need in order to disagree at all.
SLICES = {f"r{i}": dict(regions=[i], years=None, biomes=None, min_frp=400.0, n_days=6)
          for i in range(7)}
SLICES.update({
    "africa": dict(regions=[0], years=[2008], biomes=None, min_frp=0.0, n_days=20),
    "boreal": dict(regions=[4], years=None,
                   biomes=sorted(BIOME_GROUPS["forest"] | BIOME_GROUPS["grassland"]
                                 | BIOME_GROUPS["agriculture"]),
                   min_frp=400.0, n_days=40),
})


def _rows():
    with open(os.path.join(DATA, "misr_plumes_usable.csv")) as f:
        return list(csv.DictReader(f))


def _lst(r):
    """Local solar time, h. MISR's overpass is fixed by the orbit, so this is a diagnostic."""
    t = dt.datetime.fromisoformat(r["date"])
    return (t.hour + t.minute / 60.0 + float(r["lon"]) / 15.0) % 24


def select(tag):
    spec = SLICES[tag]
    rows = [r for r in _rows() if int(r["region"]) in spec["regions"]]
    if spec["years"]:
        rows = [r for r in rows if int(r["year"]) in spec["years"]]
    if spec["biomes"]:
        rows = [r for r in rows if int(r["biome"]) in spec["biomes"]]
    rows = [r for r in rows if float(r["frp_mw"]) >= spec["min_frp"]]

    by_day = collections.Counter(r["date"][:10] for r in rows)
    days = [d for d, _ in by_day.most_common(spec["n_days"])]
    keep = [r for r in rows if r["date"][:10] in days]

    lat = [float(r["lat"]) for r in keep]; lon = [float(r["lon"]) for r in keep]
    area = [min(90., max(lat) + 1), max(-180., min(lon) - 1),
            max(-90., min(lat) - 1), min(180., max(lon) + 1)]      # N, W, S, E
    hours = sorted({int(r["hour"]) for r in keep})
    lsts = [_lst(r) for r in keep]

    out = os.path.join(DATA, f"escape_slice_{tag}.json")
    json.dump(dict(tag=tag, spec=spec, levels=LEVELS, grid=GRID, area=area, hours=hours,
                   days=sorted(days), plumes=keep), open(out, "w"), indent=1)
    print(f"[select:{tag}] filters -> {len(rows)} plumes over {len(by_day)} days")
    print(f"[select:{tag}] kept the {len(days)} densest days -> {len(keep)} plumes")
    print(f"[select:{tag}] local solar time p05/50/95 = "
          f"{np.percentile(lsts, [5, 50, 95]).round(2)}  (MISR cannot reach the afternoon)")
    print(f"[select:{tag}] FRP MW p50/p90/max = "
          f"{np.percentile([float(r['frp_mw']) for r in keep], [50, 90]).round(0)} "
          f"{max(float(r['frp_mw']) for r in keep):.0f}")
    print(f"[select:{tag}] area N/W/S/E {np.round(area, 1)}  UTC hours {hours}")
    print(f"[select:{tag}] ERA5: {len({d[:7] for d in days})} year-months x {len(hours)} hours")
    print(f"[select:{tag}] wrote {out}")
    return keep, days


def era5(tag):
    """One CDS request per year-month for the selected days. Needs ~/.cdsapirc."""
    import cdsapi
    sel = json.load(open(os.path.join(DATA, f"escape_slice_{tag}.json")))
    by_month = collections.defaultdict(list)
    for d in sel["days"]:
        by_month[d[:7]].append(d[8:10])          # keyed by YYYY-MM: a slice may span years
    c = cdsapi.Client()
    jobs = [
        ("pl", "reanalysis-era5-pressure-levels",
         dict(variable=["temperature", "u_component_of_wind", "v_component_of_wind",
                        "geopotential"], pressure_level=[str(l) for l in LEVELS])),
        ("sl", "reanalysis-era5-single-levels",
         dict(variable=["boundary_layer_height", "surface_pressure"])),
    ]
    # One request per month, not one for the whole slice: CDS takes month and day as
    # independent lists and returns their cross product, so a slice spanning 6 months and 20
    # distinct day-numbers would download 480 timesteps to use 100.
    for kind, dataset, extra in jobs:
        parts = []
        for ym, mdays in sorted(by_month.items()):
            path = os.path.join(DATA, f"era5_{tag}_{kind}_{ym}.nc")
            parts.append(path)
            if os.path.exists(path):
                print(f"[era5:{tag}] cached {os.path.basename(path)}")
                continue
            req = dict(product_type="reanalysis", format="netcdf",
                       year=ym[:4], month=ym[5:7], day=sorted(set(mdays)),
                       time=[f"{h:02d}:00" for h in sel["hours"]],
                       area=sel["area"], grid=[sel["grid"], sel["grid"]], **extra)
            print(f"[era5:{tag}] requesting {kind} {ym} "
                  f"({len(set(mdays))} days x {len(sel['hours'])} hours)")
            c.retrieve(dataset, req, path)
            print(f"[era5:{tag}] wrote {os.path.basename(path)} "
                  f"({os.path.getsize(path)/1e6:.0f} MB)")
        print(f"[era5:{tag}] {kind}: {len(parts)} monthly files ready")


def era5_ml(tag):
    """Fetch L137 model levels for the slice. `reanalysis-era5-complete` takes MARS keywords."""
    import cdsapi
    sel = json.load(open(os.path.join(DATA, f"escape_slice_{tag}.json")))
    by_month = collections.defaultdict(list)
    for d in sel["days"]:
        by_month[d[:7]].append(d)
    n, w, so, e = sel["area"]
    c = cdsapi.Client()
    for ym, dates in sorted(by_month.items()):
        path = os.path.join(DATA, f"era5_{tag}_ml_{ym}.nc")
        if os.path.exists(path):
            print(f"[era5ml:{tag}] cached {os.path.basename(path)}")
            continue
        req = {"class": "ea", "expver": "1", "stream": "oper", "type": "an",
               "levtype": "ml", "levelist": "/".join(str(l) for l in ML_LEVELS),
               "param": "/".join(str(v) for v in ML_PARAMS.values()),
               "date": "/".join(sorted(set(dates))),
               "time": "/".join(f"{h:02d}:00:00" for h in sel["hours"]),
               "grid": f"{sel['grid']}/{sel['grid']}", "area": f"{n}/{w}/{so}/{e}",
               "format": "netcdf"}
        print(f"[era5ml:{tag}] requesting {ym} ({len(set(dates))} days x "
              f"{len(sel['hours'])} hours x {len(ML_LEVELS)} levels)")
        c.retrieve("reanalysis-era5-complete", req, path)
        print(f"[era5ml:{tag}] wrote {os.path.basename(path)} "
              f"({os.path.getsize(path)/1e6:.0f} MB)")


def _ab():
    """L137 half-level coefficients, indexed by half level 0..137 (0 is the implied model top)."""
    a = np.zeros(138); b = np.zeros(138)
    for row in csv.DictReader(l for l in open(AB_TABLE) if not l.startswith("#")):
        k = int(row["k"]); a[k] = float(row["a_pa"]); b[k] = float(row["b"])
    return a, b


def _column_ml(ds_ml, ds_sl, lat, lon, when, ab=None):
    """Model-level profile at one plume: heights AGL, T, p, wind -- same tuple as `_column`.

    Pressure comes from the half-level coefficients, ``p_half(k) = a(k) + b(k) ps``, with the
    surface pressure taken from the single-level file already fetched for the pressure-level
    run. Heights are the hypsometric integral upward from the surface on virtual temperature,
    which is how the levels are defined -- there is no geopotential field on model levels to
    read instead.
    """
    a, b = ab if ab is not None else _ab()
    kw = dict(latitude=lat, longitude=lon, method="nearest")
    m = ds_ml.sel(**kw).sel(valid_time=np.datetime64(when), method="nearest")
    s = ds_sl.sel(**kw).sel(valid_time=np.datetime64(when), method="nearest")
    ps = float(s["sp"])

    lev = np.asarray(m["model_level"], int)
    t = np.asarray(m["t"], float); q = np.asarray(m["q"], float)
    u = np.asarray(m["u"], float); v = np.asarray(m["v"], float)
    o = np.argsort(-lev)                                    # lowest level (137) first
    lev, t, q, u, v = lev[o], t[o], q[o], u[o], v[o]

    ph = a + b * ps                                         # half-level pressures, k=0..137
    p_full = 0.5 * (ph[lev - 1] + ph[lev])
    tv = t * (1.0 + _RV_OVER_RD_M1 * q)

    # integrate up from the surface: the lowest requested level need not be 137, so start the
    # walk at the lowest level present and report heights relative to it.
    z_half = 0.0
    z = np.empty_like(t)
    for i, k in enumerate(lev):
        dz_to_full = _RD * tv[i] / _G * np.log(ph[k] / p_full[i])
        z[i] = z_half + dz_to_full
        z_half += _RD * tv[i] / _G * np.log(ph[k] / ph[k - 1])
    return z, t, p_full, u, v, float(s["blh"]), q, ps


def _column(ds_pl, ds_sl, lat, lon, when):
    """ERA5 profile at one plume: heights AGL, theta, wind, plus the ABL depth.

    The ground is located by interpolating the geopotential-height profile to the *surface
    pressure*, not by taking the lowest pressure level. Over elevated terrain the 1000 hPa
    level sits below ground and its geopotential is an extrapolation, so referencing heights
    to it shifts the whole profile down by roughly the terrain elevation -- ~1 km over the
    Siberian slice, which is 7 K of temperature error at a given height and puts the profile
    on a different footing from ERA5's `blh`, which is genuinely above ground. Model levels do
    not have this problem: they are integrated up from the surface by construction.
    """
    kw = dict(latitude=lat, longitude=lon, method="nearest")
    p = ds_pl.sel(**kw).sel(valid_time=np.datetime64(when), method="nearest")
    s = ds_sl.sel(**kw).sel(valid_time=np.datetime64(when), method="nearest")
    lev = "pressure_level" if "pressure_level" in p.coords else "level"
    plev = np.asarray(p[lev], float) * 100.0                      # hPa -> Pa
    t = np.asarray(p["t"], float)
    z_msl = np.asarray(p["z"], float) / _G

    o = np.argsort(-plev)                                         # high pressure -> low
    z_ground = float(np.interp(np.log(float(s["sp"])), np.log(plev[o])[::-1],
                               z_msl[o][::-1]))
    z_agl = z_msl - z_ground
    order = np.argsort(z_agl)
    keep = z_agl[order] >= 0.0                                    # drop levels below ground
    return (z_agl[order][keep], t[order][keep], plev[order][keep],
            np.asarray(p["u"], float)[order][keep],
            np.asarray(p["v"], float)[order][keep], float(s["blh"]))


def _interp(z, x, zt):
    return float(np.interp(zt, z, x, left=x[0], right=x[-1]))


def score_plume(z, t, p, u, v, abl):
    """`X` -- the escape threshold in GW at C = 1, and the drivers behind it.

    Mirrors atmosphere._escape_cost_gw exactly, on a single column rather than a grid, so a
    fitted C transfers to the module without a units step in between.
    """
    theta = theta_kelvin(t - 273.15, p / 100.0)
    z_top = abl + max(_ESCAPE_EZ_FRAC * abl, _ESCAPE_EZ_MIN_M)
    th_lo, th_hi = _interp(z, theta, abl), _interp(z, theta, z_top)
    s = max(_G / (0.5 * (th_lo + th_hi)) * (th_hi - th_lo) / max(z_top - abl, 1.0), 0.0)
    inlay = z <= abl
    if not inlay.any():
        inlay = np.zeros_like(z, bool); inlay[0] = True
    u_ml = max(float(np.hypot(u[inlay].mean(), v[inlay].mean())), 0.5)
    rho = _interp(z, p, abl) / (_RD * _interp(z, t, abl))
    x_gw = ((abl / _BRIGGS_C) ** 3 * u_ml * s
            * _PI * rho * _CP_DRY * _interp(z, t, abl) / _G / 1.0e9)
    return x_gw, dict(s=s, u_ml=u_ml, rho=rho, abl=abl)


def roc(score, label):
    """ROC over a scalar score, plus Youden's J. Returns (auc, thresholds, tpr, fpr, j_star)."""
    o = np.argsort(-score)
    y = np.asarray(label, bool)[o]
    tp = np.cumsum(y); fp = np.cumsum(~y)
    tpr = tp / max(y.sum(), 1); fpr = fp / max((~y).sum(), 1)
    auc = float(np.trapezoid(np.r_[0, tpr], np.r_[0, fpr]))
    j = tpr - fpr
    k = int(np.argmax(j))
    return auc, score[o], tpr, fpr, (float(score[o][k]), float(j[k]), float(tpr[k]), float(fpr[k]))


def _scored(tag, vertical="pl"):
    """Score every plume in the slice against its reanalysis column. Cached in memory only.

    `vertical` selects which reanalysis the column is read from: "pl" (pressure levels, the
    first pass) or "ml" (L137 model levels, which resolve the entrainment zone).
    """
    import glob
    import xarray as xr
    sel = json.load(open(os.path.join(DATA, f"escape_slice_{tag}.json")))

    def _open(kind):
        parts = sorted(glob.glob(os.path.join(DATA, f"era5_{tag}_{kind}_*.nc")))
        if not parts:
            sys.exit(f"[fit] no era5_{tag}_{kind}_*.nc -- run --era5 first")
        ds = xr.open_mfdataset(parts, combine="by_coords") if len(parts) > 1 \
            else xr.open_dataset(parts[0])
        return ds.squeeze(drop=True)

    ds_sl = _open("sl")
    if vertical == "ml":
        ds_v, ab = _open("ml"), _ab()
        read = lambda la, lo, wh: _column_ml(ds_v, ds_sl, la, lo, wh, ab=ab)[:6]
    else:
        ds_v = _open("pl")
        read = lambda la, lo, wh: _column(ds_v, ds_sl, la, lo, wh)
    rec = []
    for r in sel["plumes"]:
        try:
            col = read(float(r["lat"]), float(r["lon"]),
                       dt.datetime.fromisoformat(r["date"]))
        except Exception:
            continue
        x_gw, drv = score_plume(*col)
        if not np.isfinite(x_gw) or x_gw <= 0:
            continue
        z, abl = col[0], drv["abl"]
        ez_top = abl + max(_ESCAPE_EZ_FRAC * abl, _ESCAPE_EZ_MIN_M)
        drv["levels_in_ez"] = int(((z >= abl) & (z <= ez_top)).sum())
        rec.append(dict(x_gw=x_gw, fp_raw_mw=float(r["frp_mw"]), biome=int(r["biome"]),
                        escaped=float(r["plume_top_agl_m"]) > drv["abl"], **drv))
    return rec


def fit(tag, f_rad, vertical="pl"):
    rec = _scored(tag, vertical)
    for x in rec:
        x["fp_gw"] = x["fp_raw_mw"] * 1.0e-3 / f_rad     # MW -> GW, radiative -> total
        x["r"] = x["fp_gw"] / x["x_gw"]

    if not rec:
        sys.exit("[fit] no plume could be scored -- check the ERA5 files cover the slice")
    r_ = np.array([x["r"] for x in rec]); esc = np.array([x["escaped"] for x in rec])
    print(f"[fit] scored {len(rec)} plumes  |  escaped the ABL: {esc.sum()} "
          f"({100*esc.mean():.0f} %)")
    if esc.all() or not esc.any():
        sys.exit("[fit] the outcome has one class only -- no ROC to compute")

    auc, srt, tpr, fpr, (c_star, j, tp, fp) = roc(r_, esc)
    print(f"[fit] AUC of the Briggs score = {auc:.3f}   (0.5 = the form has no skill)")
    print(f"[fit] Youden operating point: C = {c_star:.3f}  (J={j:.3f}, TPR={tp:.2f}, "
          f"FPR={fp:.2f})")
    print(f"[fit]   Briggs' own geometry is C = 1;  eq 25 expanded is C ~ 0.34")
    print(f"[fit]   f_rad = {f_rad} -- C scales linearly with it, so only C/f_rad is identified")
    out = os.path.join(DATA, f"escape_fit_{tag}_{vertical}_frad{f_rad}.json")
    json.dump(dict(f_rad=f_rad, n=len(rec), auc=auc, c_star=c_star, youden_j=j,
                   tpr=tp, fpr=fp, escaped_fraction=float(esc.mean())), open(out, "w"), indent=1)
    print(f"[fit] wrote {out}")


def _resolution_guard(rec):
    """Refuse to read a stability AUC off levels that cannot resolve the entrainment zone.

    `s` is a potential-temperature gradient across a layer of depth max(0.15*abl, 100 m) --
    ~270 m for a boreal ABL. ERA5 pressure levels at 50 hPa spacing sit ~500 m apart below
    4 km, so that layer usually contains *no* level and `s` comes out of a linear interpolation
    between levels straddling it. A null AUC for the stability term then means "not measured",
    not "does not matter", and the two must never be reported as the same thing.

    This is the same failure the module already guards against for shear
    (`atmosphere._SHEAR_MIN_LEVELS`). The fix is model levels, not more plumes.
    """
    ez = np.array([max(_ESCAPE_EZ_FRAC * r["abl"], _ESCAPE_EZ_MIN_M) for r in rec])
    n_in = np.array([r.get("levels_in_ez", np.nan) for r in rec], float)
    if not np.isfinite(n_in).any():
        return
    empty = float(np.mean(n_in[np.isfinite(n_in)] == 0))
    print(f"[ablate] entrainment zone {np.median(ez):.0f} m deep; ERA5 levels inside it: "
          f"mean {np.nanmean(n_in):.2f}, none in {100*empty:.0f} % of columns")
    if empty > 0.25:
        print("[ablate] *** the stability term is UNRESOLVED at this vertical resolution. ***")
        print("[ablate] *** Read its AUC as 'not measured', never as 'does not matter'.  ***")


def ablate(tag, vertical="pl"):
    """Decompose the score's AUC -- the step that decides whether the fit means anything.

    `r = FRP / (abl^3 U s)` and the outcome is "MISR plume top > ERA5 ABL". The ABL depth is
    therefore on *both* sides, and a shallow ABL is both easy to exceed and cheap to score. Any
    AUC has to be shown not to be that circularity before a fitted C can be believed, so this
    prints each term's AUC alone and inside strata where the confound cannot operate.
    """
    rec = _scored(tag, vertical)
    frp = np.array([r["fp_raw_mw"] for r in rec]); X = np.array([r["x_gw"] for r in rec])
    abl = np.array([r["abl"] for r in rec]); s = np.array([r["s"] for r in rec])
    u = np.array([r["u_ml"] for r in rec]); esc = np.array([r["escaped"] for r in rec])

    def au(score, name):
        v = np.isfinite(score)
        auc = roc(score[v], esc[v])[0]
        n1, n0 = esc[v].sum(), (~esc[v]).sum()
        q1, q2 = auc / (2 - auc), 2 * auc ** 2 / (1 + auc)
        se = np.sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ** 2)
                      + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0))   # Hanley & McNeil
        print(f"  {name:40s} AUC={auc:.3f} +/- {1.96*se:.3f}")

    print(f"[ablate] n={len(rec)}  escaped={esc.sum()} ({100*esc.mean():.0f} %)")
    _resolution_guard(rec)
    au(frp / X, "r = FRP/(abl^3 U s)  [the model]")
    au(frp, "FRP alone")
    au(1.0 / abl, "1/abl alone   <-- the circularity")
    au(1.0 / (u * s), "1/(U s) alone -- Briggs' atmosphere term")
    au(frp / abl ** 3, "FRP/abl^3 (drop wind and stability)")

    print("\n[ablate] AUC by firepower band (inside a band FRP barely varies, so any"
          "\n         remaining skill is the atmosphere's):")
    print(f"  {'FRP MW':>16s} {'n':>5s} {'esc%':>5s} {'AUC(r)':>7s} {'AUC(FRP)':>9s} "
          f"{'AUC 1/(Us)':>11s}")
    for lo, hi in [(0, 50), (50, 150), (150, 400), (400, 1000), (1000, np.inf)]:
        m = (frp >= lo) & (frp < hi)
        if m.sum() < 60 or esc[m].sum() < 12 or (~esc[m]).sum() < 12:
            print(f"  {lo:6.0f}-{hi:8.0f} {m.sum():5d}   too few / one class")
            continue
        print(f"  {lo:6.0f}-{hi:8.0f} {m.sum():5d} {100*esc[m].mean():5.0f} "
              f"{roc((frp/X)[m], esc[m])[0]:7.3f} {roc(frp[m], esc[m])[0]:9.3f} "
              f"{roc((1.0/(u*s))[m], esc[m])[0]:11.3f}")

    rng = np.random.default_rng(0)
    cs = [roc((frp / X)[i], esc[i])[4][0]
          for i in (rng.integers(0, len(rec), len(rec)) for _ in range(400))]
    lo, hi = np.percentile(cs, [2.5, 97.5])
    print(f"\n[ablate] bootstrap C* (f_rad=1 scale): {roc(frp/X, esc)[4][0]:.1f} "
          f"95% CI [{lo:.1f}, {hi:.1f}]")


def rungs(tag, f_rad):
    """Test all three cost rungs against MISR, each on its own target height.

    The escape test generalises: every rung predicts a firepower to reach a *specific* height,
    and MISR measures the height reached. So the outcome is per rung --

        escape    plume top >= ABL top          (ERA5 blh)
        condense  plume top >= saturation point (the rung's own z_condense)
        deep      plume top >= free convection  (the rung's own z_fc)

    -- and the score is the same ratio in each case. This calls the shipped
    :func:`pyflam.atmosphere.pyroconvection_cost_grid`, not a re-implementation, with every
    plume stacked into one (nlev, 1, N) grid so the ladder is exercised exactly as it runs in
    the pipeline.

    Class balance is reported *before* any AUC: a rung whose target no plume reaches has
    nothing to test, and saying so is the result.
    """
    import glob
    import xarray as xr
    from pyflam.atmosphere import pyroconvection_cost_grid, PYROCONVECTION_COST_RUNGS

    sel = json.load(open(os.path.join(DATA, f"escape_slice_{tag}.json")))
    op = lambda k: xr.open_mfdataset(
        sorted(glob.glob(os.path.join(DATA, f"era5_{tag}_{k}_*.nc"))),
        combine="by_coords").squeeze(drop=True)
    ds_ml, ds_sl, ab = op("ml"), op("sl"), _ab()

    cols, tops, frp = [], [], []
    for r in sel["plumes"]:
        try:
            z, t, p, u, v, blh, q, ps = _column_ml(
                ds_ml, ds_sl, float(r["lat"]), float(r["lon"]),
                dt.datetime.fromisoformat(r["date"]), ab=ab)
        except Exception:
            continue
        cols.append((z, t, p, q, u, v, ps, blh))
        tops.append(float(r["plume_top_agl_m"])); frp.append(float(r["frp_mw"]))
    n = len(cols)
    stack = lambda i: np.stack([c[i] for c in cols], axis=-1)[:, None, :]
    out = pyroconvection_cost_grid(
        stack(0), stack(1), stack(2), stack(3), stack(4), stack(5),
        surface_pressure_pa=np.array([[c[6] for c in cols]]),
        abl_m=np.array([[c[7] for c in cols]]))

    tops = np.array(tops)
    fp_gw = np.array(frp) * 1.0e-3 / f_rad
    print(f"[rungs:{tag}] {n} plumes, observed top AGL p50={np.median(tops):.0f} m, "
          f"p90={np.percentile(tops, 90):.0f} m   (f_rad={f_rad})")
    print(f"  {'rung':9s} {'target m (p50)':>14s} {'reached':>9s} {'cost GW (p50)':>13s} "
          f"{'AUC':>16s}")
    for rung in PYROCONVECTION_COST_RUNGS:
        zt = out[f"z_{rung}_m"][0]; cost = out[f"cost_{rung}_gw"][0]
        good = np.isfinite(zt) & np.isfinite(cost) & (cost > 0)
        if good.sum() < 30:
            print(f"  {rung:9s} {'--':>14s} {'--':>9s} {'--':>13s} "
                  f"{f'unresolved in {n - good.sum()}/{n}':>16s}")
            continue
        hit = tops[good] >= zt[good]
        frac = hit.mean()
        cell = f"{100*frac:.0f} % ({hit.sum()}/{good.sum()})"
        if hit.sum() < 15 or (~hit).sum() < 15:
            verdict = "NO TEST: one class"
        else:
            verdict = f"{roc(fp_gw[good] / cost[good], hit)[0]:.3f}"
        print(f"  {rung:9s} {np.median(zt[good]):14.0f} {cell:>9s} "
              f"{np.median(cost[good]):13.1f} {verdict:>16s}")


def _predicted_top(z, t, p, u, v, blh, fp_gw, *, form="differential"):
    """Solve the inequality system ``FP >= cost(z)`` for the largest admissible height.

    Treating the rungs as a *system* with one unknown, rather than as three separate equations
    each compared against its own target, is what removes the circularity: prediction and
    observation are then both heights and no target is chosen. The rungs become labels on where
    the solution lands.

    Two ways to write the demand, and they are not equivalent:

    ``differential``
        the bulk stability over the rise, ``s(z) = (g/theta_ML)(theta(z)-theta_ML)/z``, i.e. a
        *gradient* built from the endpoint excess. This is Briggs, and (shown by expanding it)
        also Tory & Kepert eq 25 -- with bulk stability the two are the same law.
    ``integral``
        the work actually done against the stratification climbed through,
        ``W(z) = int_0^z g (theta_env - theta_ML)/theta_ML dzeta``, a CIN-like quantity.

    For a linear theta profile ``W(z) = g gamma z^2 / (2 theta)`` and the two coincide exactly;
    they separate only where the column is *structured* -- inversions, residual layers. The
    integral form also averages gradient noise instead of differencing it, which matters
    because a gradient read over a few hundred metres is the noisiest thing in a column.
    """
    th = theta_kelvin(t - 273.15, p / 100.0)
    src = max(blh * 0.5, 200.0)
    th_ml = float(th[z <= src].mean()) if (z <= src).any() else float(th[0])
    rho = p / (_RD * t)
    cu = np.cumsum(u) / np.arange(1, len(u) + 1)
    cv = np.cumsum(v) / np.arange(1, len(v) + 1)
    U = np.maximum(np.hypot(cu, cv), _BRIGGS_MIN_WIND_MS)
    zz = np.maximum(z, 1.0)

    dth = th - th_ml
    k = _PI * _CP_DRY / _G / 1.0e9
    if form == "integral":
        integ = np.maximum(_G * dth / th_ml, 0.0)
        W = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(z))])
        demand = 2.0 * W / zz ** 2                      # equals s_bulk for a linear profile
    else:
        demand = np.maximum(_G / th_ml * dth / zz, 0.0)
    cost = (zz / _BRIGGS_C) ** 3 * U * demand * rho * t * k

    ok = np.isfinite(cost) & (z > src)
    if ok.sum() < 3:
        return np.nan
    zc, cc = z[ok], cost[ok]
    over = np.nonzero(cc > fp_gw)[0]
    if len(over) == 0:
        return float(zc[-1])                            # clears the sampled column
    if over[0] == 0:
        return float(src)
    i = over[0]
    f = (fp_gw - cc[i - 1]) / max(cc[i] - cc[i - 1], 1e-12)
    return float(zc[i - 1] + f * (zc[i] - zc[i - 1]))


def global_fit(f_rad, tags=("r0", "r1", "r2", "r3", "r4", "r5", "r6")):
    """Differential vs integral, pooled across every MISR region.

    The boreal slice could not separate the two forms because a mid-morning boreal column is
    close to linear in theta, where they provably coincide. Pooling the regions supplies both
    the structural variety they need to disagree and the between-regime variance a physical
    model should be able to order -- regional medians span 1305 m (Africa) to 2160 m (North
    America) against a within-region spread of ~1200 m, so a constant that wins inside one
    region has something real to lose here.
    """
    import glob
    import xarray as xr
    from scipy.stats import spearmanr

    ab = _ab()
    rows = []
    for tag in tags:
        sfile = os.path.join(DATA, f"escape_slice_{tag}.json")
        parts = sorted(glob.glob(os.path.join(DATA, f"era5_{tag}_ml_*.nc")))
        if not (os.path.exists(sfile) and parts):
            print(f"[global] {tag}: no data, skipped")
            continue
        sel = json.load(open(sfile))
        op = lambda k: (xr.open_mfdataset(sorted(glob.glob(
            os.path.join(DATA, f"era5_{tag}_{k}_*.nc"))), combine="by_coords").squeeze(drop=True))
        ml, sl = op("ml"), op("sl")
        n0 = len(rows)
        for r in sel["plumes"]:
            try:
                z, t, p, u, v, blh, q, ps = _column_ml(
                    ml, sl, float(r["lat"]), float(r["lon"]),
                    dt.datetime.fromisoformat(r["date"]), ab=ab)
            except Exception:
                continue
            fp = float(r["frp_mw"]) * 1.0e-3 / f_rad
            zd = _predicted_top(z, t, p, u, v, blh, fp, form="differential")
            zi = _predicted_top(z, t, p, u, v, blh, fp, form="integral")
            if not (np.isfinite(zd) and np.isfinite(zi)):
                continue
            rows.append((zd, zi, float(r["plume_top_agl_m"]), float(r["frp_mw"]), blh,
                         int(r["region"]), int(r["biome"])))
        print(f"[global] {tag}: +{len(rows)-n0} plumes")

    zd, zi, tops, frp, abl, reg, bio = map(np.array, zip(*rows))
    print(f"\n[global] n={len(tops)}  regions={len(set(reg))}  biomes={len(set(bio))}")
    print(f"[global] observed top IQR = {np.percentile(tops,25):.0f}-"
          f"{np.percentile(tops,75):.0f} m  (boreal-only was 1656-2651)")
    print(f"[global] the two forms differ by >100 m on "
          f"{100*np.mean(np.abs(zd-zi) > 100):.0f} % of columns")

    const = np.full_like(tops, np.median(tops))
    base = np.median(np.abs(np.log10(np.median(tops) / tops)))
    print(f"\n  {'model':26s} {'rho':>7s} {'bias m':>8s} {'MAE m':>7s} {'<=500m':>7s} "
          f"{'logerr':>7s} {'skill':>7s}")
    for name, v in [("differential (Briggs=eq25)", zd), ("integral (work done)", zi),
                    ("constant = sample median", const),
                    ("ABL depth alone", abl.astype(float)),
                    ("FRP^1/3 alone", frp ** (1 / 3))]:
        rho = spearmanr(v, tops)[0]
        if name.startswith(("ABL", "FRP")):
            print(f"  {name:26s} {rho:+7.3f} {'--':>8s} {'--':>7s} {'--':>7s} "
                  f"{'--':>7s} {'--':>7s}")
            continue
        d = v - tops
        le = np.median(np.abs(np.log10(np.maximum(v, 1) / tops)))
        print(f"  {name:26s} {rho:+7.3f} {np.median(d):8.0f} {np.abs(d).mean():7.0f} "
              f"{100*(np.abs(d)<=500).mean():6.1f}% {le:7.4f} {100*(1-le/base):+6.1f}%")

    print("\n  per-region rank correlation (can the model order the regimes?):")
    for r_ in sorted(set(reg)):
        m = reg == r_
        if m.sum() < 40:
            continue
        print(f"    region {r_}: n={m.sum():4d}  rho_diff={spearmanr(zd[m], tops[m])[0]:+.3f}  "
              f"rho_int={spearmanr(zi[m], tops[m])[0]:+.3f}  "
              f"obs p50={np.median(tops[m]):.0f} m  pred p50={np.median(zd[m]):.0f} m")
    rm = np.array([np.median(tops[reg == r_]) for r_ in sorted(set(reg))])
    pm = np.array([np.median(zd[reg == r_]) for r_ in sorted(set(reg))])
    pi = np.array([np.median(zi[reg == r_]) for r_ in sorted(set(reg))])
    print(f"  ordering the regional medians: rho_diff={spearmanr(pm, rm)[0]:+.3f}  "
          f"rho_int={spearmanr(pi, rm)[0]:+.3f}  (n={len(rm)} regions)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--select", action="store_true")
    ap.add_argument("--era5", action="store_true")
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--ablate", action="store_true")
    ap.add_argument("--f-rad", type=float, default=0.12)
    ap.add_argument("--tag", default="boreal", choices=sorted(SLICES))
    ap.add_argument("--era5-ml", action="store_true", help="fetch L137 model levels")
    ap.add_argument("--vertical", default="pl", choices=["pl", "ml"])
    ap.add_argument("--rungs", action="store_true", help="test all three rungs")
    ap.add_argument("--global-fit", action="store_true", help="pool all regions, diff vs integral")
    a = ap.parse_args()
    if a.select:
        select(a.tag)
    if a.era5:
        era5(a.tag)
    if a.era5_ml:
        era5_ml(a.tag)
    if a.fit:
        fit(a.tag, a.f_rad, a.vertical)
    if a.ablate:
        ablate(a.tag, a.vertical)
    if a.rungs:
        rungs(a.tag, a.f_rad)
    if a.global_fit:
        global_fit(a.f_rad)
    if not (a.select or a.era5 or a.era5_ml or a.fit or a.ablate or a.rungs or a.global_fit):
        ap.print_help()
