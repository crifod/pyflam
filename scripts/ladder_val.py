"""Validate the pyflam ladder against GRAF's own observed fire classes."""
import warnings; warnings.simplefilter("ignore")
import csv, glob, json, os, re
import numpy as np

import valdata
from pyflam.atmosphere import (theta_kelvin, specific_humidity_from_rh, virtual_potential_temperature,
    bulk_richardson_abl_grid, max_rh_abl_grid, parcel_mixing_depth_grid, lcl_height_bolton_m,
    dewpoint_from_rh, pyroconvection_type, PYROCONVECTION_TYPE_LEVEL)

G = lambda a: np.asarray(a, float).reshape(-1, 1, 1)
OBS = {"Surface plume": 0, "Convective plume": 1, "PyroCu": 2, "Overshooting PyroCu": 2,
       "Resilient PyroCu": 3, "Deep PyroCu / PyroCb": 4, "PyroCb": 4}

def load(path):
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
    return a if a["z"].size >= 40 else None

def diagnose(a, abl_kind):
    z, T, RH, P = a["z"], a["T"], a["RH"], a["P"]
    u = np.where(np.isfinite(a["s"]), a["s"] * np.sin(np.deg2rad(a["h"])), 0.0)
    v = np.where(np.isfinite(a["s"]), a["s"] * np.cos(np.deg2rad(a["h"])), 0.0)
    th = theta_kelvin(T, P / 100.0)
    q = specific_humidity_from_rh(RH, T + 273.15, P)
    thv = virtual_potential_temperature(th, q)
    sfc = z <= z.min() + 120
    if abl_kind == "rib":
        abl = float(bulk_richardson_abl_grid(G(z), G(thv), G(u), G(v),
              theta_v_surface=np.array([[thv[sfc].mean()]]),
              wind_u_surface=np.array([[u[sfc].mean()]]),
              wind_v_surface=np.array([[v[sfc].mean()]]))[0, 0])
    else:
        abl = float(max_rh_abl_grid(G(z), G(RH), smooth=9)[0, 0])
    if not np.isfinite(abl) or abl <= 0:
        return None
    lcl = float(lcl_height_bolton_m(T[sfc].mean() + 273.15,
                dewpoint_from_rh(T[sfc].mean(), RH[sfc].mean()) + 273.15, P[sfc].mean()))
    parcel = float(parcel_mixing_depth_grid(G(z), G(th), np.array([[th[sfc].mean()]]))[0, 0])
    inml = (z >= 80) & (z <= max(parcel, 120))
    ml = ((th[inml].mean() - th[sfc].mean()) / max(parcel - 2, 1)) if inml.sum() > 2 else 0.0
    cap = (z >= abl + 200) & (z <= abl + 1200)
    gam = float(np.polyfit(z[cap], th[cap], 1)[0]) if cap.sum() >= 2 else np.nan
    top = (z >= abl - 150) & (z <= abl + 150)
    rht = float(RH[top].mean()) if top.sum() else np.nan
    return dict(abl=abl, lcl=lcl, ratio=lcl / abl, ml=ml,
                gam=gam if np.isfinite(gam) else 3e-3, rht=rht)

def classify(d, ladder):
    kw = dict(lcl_abl_ratio=d["ratio"], ml_theta_gradient=d["ml"], gamma_theta=d["gam"], ladder=ladder)
    if ladder != "castellnou":
        kw["rh_top_abl"] = d["rht"] if np.isfinite(d["rht"]) else None
    return PYROCONVECTION_TYPE_LEVEL[pyroconvection_type(**kw)]
