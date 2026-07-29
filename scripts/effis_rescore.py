"""Capability margin with per-fire EFFIS fuel load, resolution-corrected.

The final state of the sec. 21 chain: every fire-side term measured rather than assumed --
w_a per fire from the EFFIS map (scripts/fuel_load_effis.py, agricultural cells as NFFL 3, the
confirmed assignment for Catalan cereal), dA/dt at native mapping cadence with the measured
resolution correction (scripts/coarsen.py), PFT from the sonde/ERA5 column.

Fires outside the EFFIS map (the 7 Chilean ones) keep the Tuscany default. See
docs/graf_vs_pyflam_2026-07-26.md sec. 21.13.
"""
import warnings; warnings.simplefilter("ignore")
import json, os, sys, numpy as np
from math import comb
from datetime import datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "margin_test.py")).read()
     .split('print(f"{\'window\'')[0])
from isochrone_firepower import resolution_factor, median_gap_minutes

DEFAULT_WA = 1.49


def res_factor(slug):
    r = iso.get(slug)
    if not r or not r["rates"]:
        return 1.0
    gaps = [(datetime.fromisoformat(b) - datetime.fromisoformat(a)).total_seconds() / 60
            for a, b, _ in r["rates"]]
    return resolution_factor(float(np.median(gaps)))


def score(ms, obs, lbl):
    ms = np.array(ms); o = np.array(obs); N = len(ms); npos = int((o >= 4).sum())
    tp = int(((ms > 0) & (o >= 4)).sum()); fp = int(((ms > 0) & (o < 4)).sum()); nd = tp + fp
    p = (sum(comb(npos, i) * comb(N - npos, nd - i) for i in range(tp, min(npos, nd) + 1))
         / comb(N, nd) if nd else 1.0)
    print(f"{lbl:50s} TP {tp}/{npos}  FP {fp:2d}/{N-npos}  Fisher p {p:.4f}")


base = []
for r in rows:
    g = max_growth(r["slug"], r["when"], None)
    gh = g * 3600 / 1e4 if np.isfinite(g) else np.nan
    if np.isfinite(gh) and r["crit_ha_h"] > 0:
        base.append([r["slug"], r["obs"], np.log10(gh / r["crit_ha_h"])])
for r in json.load(open(os.path.join(valdata.DATA, "oos_margin.json"))):
    if np.isfinite(r["margin"]):
        base.append([r["slug"], r["obs"], r["margin"]])
for r in base:                       # sec.18: SCQ's pyroCb was on 25 Jul, not the sonde day
    if r[0] == "santa-coloma-de-queralt":
        r[2] = -0.74

fe = valdata.load("fuel_effis")
obs = [r[1] for r in base]
final = [m + np.log10(res_factor(s)) + np.log10(((fe.get(s) or {}).get("w_a_kg_m2")
                                                 or DEFAULT_WA) / DEFAULT_WA)
         for s, o, m in base]
print(f"n={len(base)} columns, {sum(1 for o in obs if o >= 4)} pyroCb, "
      f"{sum(1 for s, _, _ in base if fe.get(s))} with EFFIS fuel\n")
score([r[2] for r in base], obs, "w_a 1.49 uniform, native rate (baseline)")
score([r[2] + np.log10(res_factor(r[0])) for r in base], obs, "+ resolution correction")
score(final, obs, "+ per-fire EFFIS w_a (agri = NFFL 3)")

print(f"\n{'fire':26s} {'obs':>3s} {'w_a':>6s} {'unclass':>8s} {'margin':>7s}")
seen = set()
for (s, o, m), f in sorted(zip(base, final), key=lambda x: -x[1]):
    if s in seen:
        continue
    seen.add(s)
    e = fe.get(s) or {}
    print(f"{s[:26]:26s} {o:3d} {(e.get('w_a_kg_m2') or DEFAULT_WA):6.2f} "
          f"{(e.get('unclassified_fraction', float('nan'))):8.2f} {f:+7.2f}")
