"""Specificity test: does the critical growth rate falsely flag the 8 non-pyroCb portal fires?

Reads the ERA5 profiles fetched by scripts/era5_oos.py, computes the PFT from reanalysis alone,
inverts it to a critical growth rate and compares against each fire's observed peak growth.
All 8 fires are observed class 1, so this tests specificity only -- there is no label variance
and no rank correlation to compute. See docs/graf_vs_pyflam_2026-07-26.md sec. 21.9.
"""
import warnings; warnings.simplefilter("ignore")
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import valdata, xarray as xr
from pyflam.atmosphere import (pyrocb_firepower_threshold_grid, critical_growth_rate_grid,
                               capability_margin, _G)
G = lambda a: np.asarray(a, float).reshape(-1, 1, 1)
W_A = 1.49
S = os.path.dirname(os.path.abspath(__file__))
iso = json.load(open("/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/"
                     "pyflam/docs/isochrone_firepower.json"))
oos = {o["slug"]: o for o in json.load(open(os.path.join(valdata.DATA,"oos_fires.json")))}
idx = json.load(open(os.path.join(valdata.DATA,"era5_oos_index.json")))

def era5_pft(nc, lat, lon):
    d = xr.open_dataset(nc)
    d = d.isel(**{k: 0 for k in d.dims if k in ("valid_time", "time", "number", "expver")})
    d = d.sel(latitude=lat, longitude=lon, method="nearest")
    lev = "pressure_level" if "pressure_level" in d.coords else "level"
    z = d["z"].values / _G
    o = np.argsort(z)
    z, T, q = z[o], d["t"].values[o], d["q"].values[o]
    u, v, p = d["u"].values[o], d["v"].values[o], d[lev].values[o] * 100.0
    zg = z - z.min()                                   # ERA5 is ASL; make it AGL-relative
    r = pyrocb_firepower_threshold_grid(G(zg), G(T), G(p), G(q), G(u), G(v),
                                        surface_pressure_pa=np.array([[p.max()]]))
    return float(np.asarray(r["pft_gw"]).ravel()[0]), float(T[-1] - 273.15)

print(f"{'fire':16s} {'obs':>3s} {'PFT GW':>8s} {'crit ha/h':>10s} "
      f"{'obs ha/h':>9s} {'margin':>7s}  verdict")
rows = []
for slug, meta in sorted(idx.items()):
    if not os.path.exists(meta["path"]):
        print(f"{slug[:16]:16s}  -- no ERA5"); continue
    pft, ttop = era5_pft(meta["path"], oos[slug]["lat"], oos[slug]["lon"])
    crit = float(critical_growth_rate_grid(np.array([pft]), fuel_load_kg_m2=W_A)[0])
    pk = iso.get(slug, {}).get("peak_ha_per_h") or oos[slug].get("br_max")
    m = float(capability_margin(pk, crit)) if pk else np.nan
    ok = "ok (below)" if np.isfinite(m) and m <= 0 else (
         "FALSE POSITIVE" if np.isfinite(m) else "n/a")
    print(f"{slug[:16]:16s} {oos[slug]['obs']:3d} {pft:8.0f} {crit:10.0f} "
          f"{(pk if pk else float('nan')):9.0f} {m:+7.2f}  {ok}")
    rows.append(dict(slug=slug, obs=oos[slug]["obs"], pft=pft, crit=crit, peak=pk, margin=m))
fin = [r for r in rows if np.isfinite(r["margin"])]
fp = [r for r in fin if r["margin"] > 0]
print(f"\n{len(fin)} fires tested, all observed class 1 (non-pyroCb)")
print(f"false positives (margin > 0): {len(fp)}"
      + (": " + ", ".join(r['slug'] for r in fp) if fp else ""))
if fin:
    print(f"margin range: {min(r['margin'] for r in fin):+.2f} to "
          f"{max(r['margin'] for r in fin):+.2f}")
json.dump(rows, open(os.path.join(valdata.DATA,"oos_margin.json"), "w"), indent=1)
