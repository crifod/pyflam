"""Do C-Haines and the PFT order the same columns the same way?

The two diagnostics in routine use weight the atmosphere differently, and the disagreement
has never been measured on identical columns:

* **C-Haines** = CA + CB, where CA is the 850->700 hPa lapse. A *strong* lapse scores high.
* **PFT** = 0.3 z_fc^2 U_ML dtheta_fc. A *weak* lid means small dtheta_fc and a low threshold,
  i.e. pyroCb is cheap.

So a weakly-capped column is unfavourable to one and cheap to the other. Until now the two
could not be compared, because the daily runner exported the PFT fields but not C-Haines.
It does now (``diag_chaines_*.tif``), so this script scores both on the same cells.

What it reports, over the classifiable land of a 3-day run:

1. rank correlation between C-Haines and the PFT (and its inverse, the "ease" of pyroCb);
2. the same restricted to the **degenerate** population -- cells whose PFT has collapsed --
   which is where the frameworks should part company most sharply;
3. how the two rank the cells the fuel-gated product actually flags (class >= 3), i.e.
   whether either diagnostic would have picked out the operationally interesting cells.

Usage:  PYTHONPATH=src python scripts/chaines_vs_pft.py [forecast-folder] [rundate]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import rasterio
from scipy.stats import spearmanr, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    REPO, "docs", "forecast_2026-07-30_3day")
RUNDATE = sys.argv[2] if len(sys.argv) > 2 else "2026-07-30"
DAYS = ["2026-07-30", "2026-07-31", "2026-08-01"]
HOURS = [9, 12, 15, 18]                       # daytime only: the classifier is blank at night
ABL_MIN_M = 150.0
# Below this the PFT has degenerated -- the population identified in the 2026-07-30 run, where
# threshold cells carried PFTs of a few GW against a domain median near 90.
DEGENERATE_PFT_GW = 10.0


def band(day, name, hour):
    p = os.path.join(OUT, f"rasters_hybrid_{day}", f"diag_{name}_{hour:02d}Z.tif")
    with rasterio.open(p) as d:
        return d.read(1).astype("float64")


def cls(day, kind, hour):
    p = os.path.join(OUT, f"rasters_hybrid_{day}", f"pyroconv_{kind}_{hour:02d}Z.tif")
    with rasterio.open(p) as d:
        return d.read(1)


def collect():
    """Every classifiable daytime land column of the run, as flat arrays."""
    ch, pft, marg, abl, gate, zfc, dth = [], [], [], [], [], [], []
    for day in DAYS:
        for h in HOURS:
            a = band(day, "abl", h)
            ok = np.isfinite(a) & (a >= ABL_MIN_M)
            for dst, name in ((ch, "chaines"), (pft, "pft_gw"), (marg, "pft_margin"),
                              (zfc, "z_fc"), (dth, "delta_theta_fc")):
                dst.append(band(day, name, h)[ok])
            abl.append(a[ok])
            gate.append(cls(day, "gated", h)[ok])
    f = lambda x: np.concatenate(x)
    return dict(ch=f(ch), pft=f(pft), margin=f(marg), abl=f(abl),
                gate=f(gate), z_fc=f(zfc), dtheta=f(dth))


def main():
    d = collect()
    fin = np.isfinite(d["ch"]) & np.isfinite(d["pft"]) & (d["pft"] > 0)
    ch, pft = d["ch"][fin], d["pft"][fin]
    gate, zfc, dth = d["gate"][fin], d["z_fc"][fin], d["dtheta"][fin]
    print(f"columns compared: {fin.sum():,} (classifiable daytime land, {len(DAYS)} days x "
          f"{len(HOURS)} hours)\n")

    # --- 1. do they agree at all? -----------------------------------------------------
    r, p = spearmanr(ch, pft)
    print("1. C-Haines vs PFT, all columns")
    print(f"   Spearman rho = {r:+.3f}  (p = {p:.2e})")
    print("   PFT low = pyroCb cheap, so AGREEMENT would be rho < 0 "
          "(high C-Haines <-> low threshold).")
    print(f"   -> {'agree' if r < -0.1 else 'DISAGREE' if r > 0.1 else 'no relationship'}\n")

    # which term of the PFT drives that? ----------------------------------------------
    for nm, v in (("z_fc", zfc), ("delta_theta_fc", dth)):
        rr, _ = spearmanr(ch, v)
        print(f"   C-Haines vs {nm:15s} rho = {rr:+.3f}")
    print()

    # --- 2. the degenerate population -------------------------------------------------
    deg = pft < DEGENERATE_PFT_GW
    print(f"2. Degenerate columns (PFT < {DEGENERATE_PFT_GW:g} GW): {deg.sum():,} "
          f"({100*deg.mean():.1f}% of columns)")
    if deg.sum() > 30:
        u, pu = mannwhitneyu(ch[deg], ch[~deg])
        print(f"   median C-Haines  degenerate {np.median(ch[deg]):.2f}  vs  "
              f"rest {np.median(ch[~deg]):.2f}   (Mann-Whitney p = {pu:.2e})")
        verdict = ("LOWER -- C-Haines calls unfavourable what the PFT calls cheap"
                   if np.median(ch[deg]) < np.median(ch[~deg]) else
                   "higher -- the two agree on this population")
        print(f"   -> the cells the PFT finds cheapest score {verdict}\n")

    # --- 3. would either diagnostic have found the flagged cells? ---------------------
    hit = gate >= 3
    print(f"3. Cells the fuel-gated product flags class >= 3: {hit.sum()}")
    if hit.sum() > 0:
        for nm, v, hi_is_bad in (("C-Haines", ch, False), ("PFT [GW]", pft, True)):
            pct = 100.0 * (v[~hit] < v[hit].mean()).mean() if hi_is_bad else \
                  100.0 * (v[~hit] > v[hit].mean()).mean()
            print(f"   {nm:10s} flagged mean {v[hit].mean():8.2f}   "
                  f"all-column mean {v.mean():8.2f}   "
                  f"non-flagged columns beyond that mean: {pct:.1f}%")
        print("   (a diagnostic with skill would leave few non-flagged columns beyond its "
              "own flagged mean)")


if __name__ == "__main__":
    main()
