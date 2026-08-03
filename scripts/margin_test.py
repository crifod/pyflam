"""Fire-level capability margin: observed max growth vs the growth this atmosphere demands.

Both sides fire-level, so the published fire-level labels are commensurable. Continuous, so
all sondes contribute rather than only where a binary gate fires.
"""
import warnings; warnings.simplefilter("ignore")
import sys, glob, json, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import valdata
from datetime import datetime, timedelta
from ladder_val import load, diagnose, OBS, G
from splice_pft import pft_spliced
from match import best
from pyflam.atmosphere import pyroconvection_type, PYROCONVECTION_TYPE_LEVEL

S = os.path.dirname(os.path.abspath(__file__))
W_A, ALPHA, HEAT = 1.49, 0.7, 15.0e6
LADDER = "adaptive"
times  = json.load(open(os.path.join(valdata.DATA,"sonde_times.json")))
labels = {f['slug']: f['obs_class'] for f in json.load(open(os.path.join(valdata.DATA,"fire_labels.json")))}
port   = {p['slug']: p for p in json.load(open(os.path.join(valdata.DATA,"portal_fires.json")))}
iso    = json.load(open('docs/isochrone_firepower.json'))

def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
    n = x.size
    rx = np.argsort(np.argsort(x)).astype(float); ry = np.argsort(np.argsort(y)).astype(float)
    r = np.corrcoef(rx, ry)[0, 1]
    t = r * np.sqrt(max(n - 2, 1) / max(1 - r**2, 1e-12))       # t approx, df = n-2
    from math import erf, sqrt
    p = 2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2))))
    return r, p, n

def pft_of(path, slug, when):
    nc = valdata.era5_for(slug, when) or ''
    if not os.path.exists(nc):
        cand = sorted(glob.glob(os.path.join(valdata.ERA5, f"era5*_{slug}_{when:%Y%m%d}_*.nc")))
        if not cand: return float("nan")
        nc = min(cand, key=lambda c: abs(int(c[-5:-3]) - when.hour))
    try:
        r, _, _, _ = pft_spliced(path, nc, port[slug]["lat"], port[slug]["lon"])
    except Exception: return float("nan")
    return float(np.asarray(r["pft_gw"]).ravel()[0])

def max_growth(slug, when, hours):
    """Peak dA/dt (m2/s) within +/-hours of launch; None window = whole fire."""
    if slug not in iso: return float("nan")
    best_r = float("nan")
    for t0, t1, r in iso[slug]['rates']:
        a, b = datetime.fromisoformat(t0), datetime.fromisoformat(t1)
        if hours is not None:
            if b < when - timedelta(hours=hours) or a > when + timedelta(hours=hours):
                continue
        best_r = r if not np.isfinite(best_r) else max(best_r, r)
    return best_r

rows = []
for f in sorted(glob.glob(valdata.SONDES+'/*.csv')):
    b = os.path.basename(f)
    slug = best(re.match(r'\d+-(WF|PF)_(.+?)(_Ambient.*)?\.csv$', b).group(2), list(port))
    if not slug or slug not in labels: continue
    k = re.match(r'(\d+)-', b).group(1)
    if k not in times: continue
    a = load(f)
    if a is None: continue
    d = diagnose(a, "rib")
    if d is None: continue
    t = times[k]; dd, mo, y = t['date'].split('/')
    when = datetime(int(y), int(mo), int(dd), int(t['time'][:2]), int(t['time'][3:5]))
    pft = pft_of(f, slug, when)
    kw = dict(lcl_abl_ratio=d["ratio"], ml_theta_gradient=d["ml"], gamma_theta=d["gam"],
              ladder=LADDER, rh_top_abl=d["rht"] if np.isfinite(d["rht"]) else None)
    rows.append(dict(slug=slug, obs=OBS[labels[slug]], pft=pft, when=when, path=f,
                     pred=PYROCONVECTION_TYPE_LEVEL[pyroconvection_type(**kw)],
                     crit_ha_h=pft*1e9/(ALPHA*HEAT*W_A)*3600/1e4))

print(f"{'window':16s} {'n':>3s} {'rho':>7s} {'p':>9s}")
for lbl, hrs in (("+/- 2 h", 2), ("+/- 4 h", 4), ("+/- 8 h", 8), ("whole fire", None)):
    marg, obs = [], []
    for r in rows:
        g = max_growth(r["slug"], r["when"], hrs)
        m = (np.log10(max(g, 1e-9)*3600/1e4 / r["crit_ha_h"])
             if np.isfinite(g) and np.isfinite(r["crit_ha_h"]) and r["crit_ha_h"] > 0 else np.nan)
        marg.append(m); obs.append(r["obs"])
    rho, p, n = spearman(marg, obs)
    print(f"{lbl:16s} {n:3d} {rho:+7.3f} {p:9.4f}")

# reference: how well does the ladder itself rank the observed class?
rho, p, n = spearman([r["pred"] for r in rows], [r["obs"] for r in rows])
print(f"\n{'ladder prediction':16s} {n:3d} {rho:+7.3f} {p:9.4f}   (reference)")
print(f"\n{'fire':26s} {'obs':>3s} {'crit ha/h':>10s} {'max ha/h (whole)':>17s} {'margin':>7s}")
for r in sorted(rows, key=lambda x: -x["obs"]):
    g = max_growth(r["slug"], r["when"], None)
    gh = g*3600/1e4 if np.isfinite(g) else np.nan
    m = np.log10(gh/r["crit_ha_h"]) if np.isfinite(gh) and r["crit_ha_h"]>0 else np.nan
    print(f"{r['slug'][:26]:26s} {r['obs']:3d} {r['crit_ha_h']:10.0f} {gh:17.0f} {m:+7.2f}")

from collections import Counter
obs = [r["obs"] for r in rows]; pred = [r["pred"] for r in rows]
c = Counter(obs)
print(f"\nclass balance (obs): {dict(sorted(c.items()))}")
maj, majn = c.most_common(1)[0]
print(f"majority-class baseline: always predict {maj} -> {majn}/{len(obs)} = {majn/len(obs):.0%}")
print(f"ladder exact           : {sum(p==o for p,o in zip(pred,obs))}/{len(obs)} = "
      f"{sum(p==o for p,o in zip(pred,obs))/len(obs):.0%}")
print(f"ladder predictions     : {dict(sorted(Counter(pred).items()))}")

# binary detection of pyroCb (class 4) by the capability margin, whole-fire window
m = []
for r in rows:
    g = max_growth(r["slug"], r["when"], None)
    gh = g*3600/1e4 if np.isfinite(g) else np.nan
    m.append(np.log10(gh/r["crit_ha_h"]) if np.isfinite(gh) and r["crit_ha_h"]>0 else np.nan)
m = np.array(m); o = np.array(obs); ok = np.isfinite(m)
for thr_cls in (4, 3):
    tp = int(((m[ok] > 0) & (o[ok] >= thr_cls)).sum()); fp = int(((m[ok] > 0) & (o[ok] < thr_cls)).sum())
    fn = int(((m[ok] <= 0) & (o[ok] >= thr_cls)).sum()); tn = int(((m[ok] <= 0) & (o[ok] < thr_cls)).sum())
    from math import comb
    npos = tp + fn; ndraw = tp + fp; N = tp+fp+fn+tn
    p = sum(comb(npos,i)*comb(N-npos,ndraw-i) for i in range(tp, min(npos,ndraw)+1))/comb(N,ndraw)
    print(f"\nmargin>0 detects class>={thr_cls}: TP={tp} FP={fp} FN={fn} TN={tn}  "
          f"(n={N}, one-tailed Fisher p={p:.3f})")
