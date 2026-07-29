"""How much of the peak growth rate does a mapping window destroy? Measure it, don't assume it.

Sec. 21.10 charged part of the fire-side gap to a dA/dt averaged over a ~60 min mapping
interval, with an assumed 1.5-2x penalty. This measures it: take the fires mapped at <=30 min,
coarsen their own cumulative-area curves, and fit the decay. Median beta = 0.44, and a 60 min
window costs 2.21x against a 10 min reference. See docs/graf_vs_pyflam_2026-07-26.md sec. 21.11.

Reads docs/isochrone_firepower.json; run from the repo root.
"""
import json, numpy as np
from datetime import datetime

iso = json.load(open('docs/isochrone_firepower.json'))

def cumulative(rates):
    """[(t0,t1,dA/dt)] -> (t_minutes, cumulative_area_m2), piecewise-linear."""
    t = [0.0]; a = [0.0]
    t0 = datetime.fromisoformat(rates[0][0])
    for s0, s1, r in rates:
        d0 = (datetime.fromisoformat(s0) - t0).total_seconds()
        d1 = (datetime.fromisoformat(s1) - t0).total_seconds()
        if d0 > t[-1] + 1:                      # gap in mapping: carry area flat
            t.append(d0); a.append(a[-1])
        t.append(d1); a.append(a[-1] + r * (d1 - d0))
    return np.array(t) / 60.0, np.array(a)

def peak_at(tm, ar, w_min):
    """Max mean growth (ha/h) over any window of w_min minutes."""
    if tm[-1] < w_min: return np.nan
    grid = np.arange(tm[0], tm[-1] - w_min + 1e-9, 1.0)
    if grid.size == 0: return np.nan
    A = np.interp(grid, tm, ar); B = np.interp(grid + w_min, tm, ar)
    return float(np.max(B - A) / (w_min / 60.0) / 1e4)

print(f"{'fire':22s} {'native':>7s} {'peak@native':>12s} {'@30':>9s} {'@60':>9s} "
      f"{'@120':>9s} {'nat/60':>7s}")
fac = []
for s, r in sorted(iso.items()):
    rates = r['rates']
    if len(rates) < 4: continue
    gaps = [(datetime.fromisoformat(b) - datetime.fromisoformat(a)).total_seconds()/60
            for a, b, _ in rates]
    nat = float(np.median(gaps))
    if nat > 30: continue                       # too coarse to say anything about 60 min
    tm, ar = cumulative(rates)
    pn, p30, p60, p120 = (peak_at(tm, ar, nat), peak_at(tm, ar, 30),
                          peak_at(tm, ar, 60), peak_at(tm, ar, 120))
    ratio = pn / p60 if (p60 and np.isfinite(p60) and p60 > 0) else np.nan
    if np.isfinite(ratio): fac.append(ratio)
    print(f"{s[:22]:22s} {nat:7.0f} {pn:12.0f} {p30:9.0f} {p60:9.0f} {p120:9.0f} {ratio:7.2f}")
fac = np.array(fac)
print(f"\nnative-resolution peak / 60-min peak: median {np.median(fac):.2f}, "
      f"range {fac.min():.2f}-{fac.max():.2f}  (n={fac.size} fires)")
print("linear interpolation inside each interval, so sub-interval bursts are invisible:")
print("this is a LOWER BOUND on the true instantaneous enhancement.")
