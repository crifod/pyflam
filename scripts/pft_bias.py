"""Is the PFT systematically high? Sweep a uniform scale factor over all 30 fire-observations.

Only 1 of 30 columns exceeds its PFT, which reads as a threshold biased high. Sweeping a scale
factor separates the two readings: if a modest uniform correction captures both observed pyroCb
without admitting many false positives, the bar is misplaced; if it does not, the gap is
elsewhere. Combined with the independent anchoring of PFT magnitude, the answer is that the gap
sits on the firepower side -- docs/graf_vs_pyflam_2026-07-26.md sec. 21.10.

Depends on the campaign-side harness (margin_test.py) and its scratchpad inputs; it is the
record of how the numbers in sec. 21.10 were produced, not a reusable library.
"""
import warnings; warnings.simplefilter("ignore")
import json, os, sys, numpy as np
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"margin_test.py")).read()
     .split('print(f"{\'window\'')[0])

comb_rows = []
for r in rows:
    g = max_growth(r["slug"], r["when"], None)
    gh = g*3600/1e4 if np.isfinite(g) else np.nan
    if np.isfinite(gh) and r["crit_ha_h"] > 0:
        comb_rows.append([r["slug"], r["obs"], np.log10(gh/r["crit_ha_h"])])
for r in json.load(open("/tmp/oos_margin.json")):
    if np.isfinite(r["margin"]): comb_rows.append([r["slug"], r["obs"], r["margin"]])

# SCQ: sec.18 establishes the pyroCb was on 25 Jul, not the 24 Jul sonde day. Use the 25th's
# most favourable hour, consistent with taking the fire's max growth over its life.
SCQ_CORRECTED = -0.74
for r in comb_rows:
    if r[0] == "santa-coloma-de-queralt": r[2] = SCQ_CORRECTED

m = np.array([r[2] for r in comb_rows]); o = np.array([r[1] for r in comb_rows])
npos = int((o >= 4).sum()); N = len(m)
print(f"n={N}, pyroCb={npos} (SCQ corrected to its 25 Jul atmosphere)\n")
print(f"{'PFT scale':>9s} {'threshold':>9s} {'TP':>3s} {'FP':>3s} {'FN':>3s} {'TN':>3s} "
      f"{'precision':>9s} {'Fisher p':>9s}")
for k in (1.0, 2.0, 3.0, 5.5, 8.0, 12.0, 20.0, 50.0):
    t = -np.log10(k)                        # scaling PFT by 1/k moves the bar to -log10(k)
    tp = int(((m > t) & (o >= 4)).sum()); fp = int(((m > t) & (o < 4)).sum())
    fn = npos - tp; tn = N - npos - fp
    nd = tp + fp
    p = (sum(comb(npos,i)*comb(N-npos,nd-i) for i in range(tp, min(npos,nd)+1))/comb(N,nd)
         if nd else 1.0)
    print(f"{k:9.1f} {t:+9.2f} {tp:3d} {fp:3d} {fn:3d} {tn:3d} "
          f"{(tp/nd if nd else float('nan')):9.2f} {p:9.4f}")
print("\nfires admitted as the bar drops (non-pyroCb first to enter):")
for s, ob, mm in sorted(comb_rows, key=lambda x: -x[2])[:8]:
    print(f"   {s[:26]:26s} obs={ob} margin={mm:+.2f}")
