"""Apply the measured resolution correction to every labelled column and rescore.

The correction is real and physically motivated (scripts/coarsen.py), but applied across the 30
labelled columns it crosses no decision boundary: the gap it closes is near-uniform across
observed classes, so it is calibration and not discrimination. Adding a w_a correction on top
makes the score worse. This is the record of that negative result --
docs/graf_vs_pyflam_2026-07-26.md sec. 21.11.

Reads the vendored corpus via scripts/valdata.py; runs from a clone with no external inputs.
"""
import warnings; warnings.simplefilter("ignore")
import json

import valdata, os, sys, numpy as np
from math import comb
from datetime import datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "margin_test.py")).read()
     .split('print(f"{\'window\'')[0])
BETA, WREF, CAP = 0.44, 10.0, 3.0

def res_factor(slug):
    r = iso.get(slug)
    if not r or not r['rates']: return 1.0, np.nan
    g = [(datetime.fromisoformat(b)-datetime.fromisoformat(a)).total_seconds()/60
         for a, b, _ in r['rates']]
    nat = float(np.median(g))
    return min((max(nat, WREF)/WREF)**BETA, CAP), nat

cr = []
for r in rows:
    g = max_growth(r["slug"], r["when"], None); gh = g*3600/1e4 if np.isfinite(g) else np.nan
    if np.isfinite(gh) and r["crit_ha_h"] > 0:
        f, nat = res_factor(r["slug"])
        cr.append([r["slug"], r["obs"], np.log10(gh/r["crit_ha_h"]), f, nat])
for r in json.load(open(os.path.join(valdata.DATA,"oos_margin.json"))):
    if np.isfinite(r["margin"]):
        f, nat = res_factor(r["slug"]); cr.append([r["slug"], r["obs"], r["margin"], f, nat])
for r in cr:
    if r[0] == "santa-coloma-de-queralt": r[2] = -0.74      # sec.18: pyroCb on 25 Jul

def score(ms, os_, lbl):
    ms = np.array(ms); os_ = np.array(os_); N = len(ms); npos = int((os_ >= 4).sum())
    tp = int(((ms > 0) & (os_ >= 4)).sum()); fp = int(((ms > 0) & (os_ < 4)).sum()); nd = tp+fp
    p = (sum(comb(npos, i)*comb(N-npos, nd-i) for i in range(tp, min(npos, nd)+1))/comb(N, nd)
         if nd else 1.0)
    print(f"{lbl:46s} TP {tp}/{npos}  FP {fp:2d}/{N-npos}  Fisher p {p:.4f}")

obs = [r[1] for r in cr]
print(f"n={len(cr)} columns, {sum(1 for o in obs if o>=4)} pyroCb\n")
score([r[2] for r in cr], obs, "baseline (w_a 1.49, native isochrone rate)")
score([r[2]+np.log10(r[3]) for r in cr], obs, "+ measured resolution correction (beta 0.44)")
score([r[2]+np.log10(r[3])+np.log10(4.0/1.49) for r in cr], obs, "+ resolution + w_a 4.0 kg/m2")
print(f"\n{'fire':26s} {'obs':>3s} {'gap min':>8s} {'xfac':>5s} {'base':>6s} {'+res':>6s} {'+both':>6s}")
seen = set()
for s, o, m, f, nat in sorted(cr, key=lambda x: -(x[2]+np.log10(x[3]))):
    if s in seen: continue
    seen.add(s)
    print(f"{s[:26]:26s} {o:3d} {nat:8.0f} {f:5.2f} {m:+6.2f} {m+np.log10(f):+6.2f} "
          f"{m+np.log10(f)+np.log10(4.0/1.49):+6.2f}")
