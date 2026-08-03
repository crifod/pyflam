"""ICON-EU pressure-levels vs model-levels: which fixes the ML gradient?

Degrades IGRA soundings onto ICON-2I (5 p-lev), ICON-EU (20 p-lev) and ICON-EU model-level
heights, and compares the recovered mixed-layer dtheta/dz to full-resolution truth. Run with
PYTHONPATH=../../src after fetching soundings (see README).
"""

import numpy as np, sys
from sounding_degradation import parse_igra, prep, parcel_abl, ml_grad_inside, ELEV
from pyflam.atmosphere import theta_kelvin
THR = 1.1e-3

# Approximate ICON-EU model-level full heights AGL (SLEVE, 74 levels): fine near the
# surface, stretching aloft. Lowest ~10 m; ~20 levels below 2 km.
ICON_EU_ML = np.array([10,30,55,85,120,160,205,255,310,375,450,535,630,735,850,975,
                       1110,1260,1425,1605,1800,2010,2235,2475,2730,3000,3300,3620,
                       3960,4320,4700,5100,5520], float)
# ICON-2I pressure levels, as heights (what we have today)
ICON_2I_Z = np.array([110, 770, 1457, 3012, 5574], float)
# ICON-EU pressure levels, as approximate heights
ICON_EU_P = np.array([110, 540, 770, 1010, 1140, 1460, 1620, 1930, 2270, 3010, 4200, 5570], float)

def sample_on_heights(z, t, p, zt):
    keep = zt[(zt >= z.min()) & (zt <= z.max())]
    if keep.size < 2: return None
    tt = np.interp(keep, z, t)
    pp = np.interp(keep, z, p)
    return keep, tt, pp

def fit_in_ml(p, z, t, ml_top, zmin=80.0):
    th = theta_kelvin(t - 273.15, p)
    ins = np.where((z <= ml_top) & (z >= zmin))[0]
    if ins.size < 2: return np.nan, ins.size
    return float(np.polyfit(z[ins], th[ins], 1)[0]), ins.size

for stn, name in [("ITM00016144","INLAND (Capofiume)"), ("ITM00016245","COAST (Pratica)")]:
    rows=[]
    for s in parse_igra(f"por_{stn}/{stn}-data.txt", months=(6,7,8), hour=12):
        try: p,z,t,rh,u,v = prep(s, ELEV[stn])
        except Exception: continue
        if len(p)<30 or z.max()<4000: continue
        ml = parcel_abl(p,z,t)
        if ml < 300: continue
        truth = ml_grad_inside(p,z,t,ml)
        if not np.isfinite(truth): continue
        rec=[truth]
        for grid in (ICON_2I_Z, ICON_EU_P, ICON_EU_ML):
            out = sample_on_heights(z,t,p,grid)
            if out is None: rec += [np.nan, 0]; continue
            zz, tt, pp = out
            mlc = parcel_abl(pp, zz, tt)
            g, n = fit_in_ml(pp, zz, tt, mlc)
            rec += [g, n]
        rows.append(rec)
    r=np.array(rows,float); tr = r[:,0] <= THR
    print(f"\n{'='*88}\n{name}   n={len(r)}   truth median = {np.median(r[:,0]):.2e}   capable base rate {100*tr.mean():.0f}%\n{'='*88}")
    print(f"  {'grid':<24} {'lev in ML':>9} {'bias':>10} {'MAE':>10} {'r':>6} {'sens':>5} {'spec':>5} {'J':>6}")
    for lbl, gi, ni in [("ICON-2I  (5 p-lev)",1,2), ("ICON-EU  (20 p-lev)",3,4), ("ICON-EU  (model levels)",5,6)]:
        g=r[:,gi]; nl=r[:,ni]; m=np.isfinite(g)&np.isfinite(r[:,0])
        if m.sum()<20: print(f"  {lbl:<24} too few"); continue
        bias=np.median(g[m]-r[m,0]); mae=np.median(np.abs(g[m]-r[m,0]))
        cc=np.corrcoef(r[m,0],g[m])[0,1]
        pred=g[m]<=THR; t2=tr[m]
        tp=int((t2&pred).sum()); tn=int((~t2&~pred).sum()); fp=int((~t2&pred).sum()); fn=int((t2&~pred).sum())
        se=tp/max(tp+fn,1); sp=tn/max(tn+fp,1)
        print(f"  {lbl:<24} {np.nanmean(nl):9.1f} {bias:10.2e} {mae:10.2e} {cc:6.2f} {se:5.2f} {sp:5.2f} {se+sp-1:6.2f}")
