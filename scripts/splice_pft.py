import warnings; warnings.simplefilter("ignore")
import numpy as np

import valdata, csv, sys, xarray as xr
from pyflam.atmosphere import pyrocb_firepower_threshold_grid, specific_humidity_from_rh, _G
G = lambda a: np.asarray(a, float).reshape(-1, 1, 1)

def load_sonde(path):
    cols = {"z": "Altitude (m AGL)", "T": "Temperature (C)", "RH": "Relative humidity (%)",
            "P": "Pressure (Pascal)", "s": "Speed (m/s)", "h": "Heading (degrees)"}
    acc = {k: [] for k in cols}
    with open(path) as f:
        rd = csv.reader(f); hdr = [c.strip() for c in next(rd)]
        ix = {k: hdr.index(v) for k, v in cols.items()}
        for r in rd:
            for k in cols:
                try: acc[k].append(float(r[ix[k]]))
                except Exception: acc[k].append(np.nan)
    a = {k: np.array(v, float) for k, v in acc.items()}
    ok = (np.isfinite(a["z"]) & np.isfinite(a["T"]) & np.isfinite(a["P"]) & (a["P"] > 1e4)
          & np.isfinite(a["RH"]) & (a["z"] >= 0))
    a = {k: v[ok] for k, v in a.items()}
    o = np.argsort(a["z"]); a = {k: v[o] for k, v in a.items()}
    a["RH"] = np.clip(a["RH"], 1, 100)
    a["u"] = np.where(np.isfinite(a["s"]), a["s"] * np.sin(np.deg2rad(a["h"])), 0.0)
    a["v"] = np.where(np.isfinite(a["s"]), a["s"] * np.cos(np.deg2rad(a["h"])), 0.0)
    a["q"] = specific_humidity_from_rh(a["RH"], a["T"] + 273.15, a["P"])
    return a

def era5_profile(nc, lat, lon):
    d = xr.open_dataset(nc)
    d = d.isel(**{k: 0 for k in d.dims if k in ("valid_time", "time", "number", "expver")})
    d = d.sel(latitude=lat, longitude=lon, method="nearest")
    lev = "pressure_level" if "pressure_level" in d.coords else "level"
    z = d["z"].values / _G
    o = np.argsort(z)
    return dict(z=z[o], T=d["t"].values[o], q=d["q"].values[o],
                u=d["u"].values[o], v=d["v"].values[o], p=d[lev].values[o] * 100.0)

def pft_spliced(sonde_path, nc, lat, lon):
    s = load_sonde(sonde_path); e = era5_profile(nc, lat, lon)
    zg = e["z"] - e["z"].min()
    ab = zg > s["z"].max()
    zz = np.concatenate([s["z"], zg[ab]]); TT = np.concatenate([s["T"] + 273.15, e["T"][ab]])
    qq = np.concatenate([s["q"], e["q"][ab]]); PP = np.concatenate([s["P"], e["p"][ab]])
    uu = np.concatenate([s["u"], e["u"][ab]]); vv = np.concatenate([s["v"], e["v"][ab]])
    r = pyrocb_firepower_threshold_grid(G(zz), G(TT), G(PP), G(qq), G(uu), G(vv),
                                        surface_pressure_pa=np.array([[s["P"][0]]]))
    return r, s["z"].max(), int(ab.sum()), TT[-1] - 273.15
