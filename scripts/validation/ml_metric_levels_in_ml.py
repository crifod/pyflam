"""Candidate D: gradient between two REAL levels both inside the mixed layer.
Plus: how many model levels are actually inside the ML? That number decides whether
the mixed-layer gradient is measurable at all from this archive.
"""
import numpy as np, sys
from sounding_degradation import (parse_igra, prep, abl_rib, parcel_abl, ml_grad_inside,
                     ml_grad_surface_to_abl, degrade, ELEV)
from pyflam.atmosphere import theta_kelvin
THR=1.1e-3

def cand_D(p,z,t,ml_top,zmin=80.0):
    th=theta_kelvin(t-273.15,p)
    ins=np.where((z<=ml_top)&(z>=zmin))[0]
    if ins.size<2: return np.nan
    return float(np.polyfit(z[ins],th[ins],1)[0])

def rates(tr,pr):
    tp=int((tr&pr).sum()); tn=int((~tr&~pr).sum()); fp=int((~tr&pr).sum()); fn=int((tr&~pr).sum())
    se=tp/max(tp+fn,1); sp=tn/max(tn+fp,1); return se,sp,(tp+tn)/max(len(tr),1),se+sp-1

for stn,name in [("ITM00016144","INLAND (Capofiume)"),("ITM00016245","COAST (Pratica)")]:
    R=[]; nin=[]
    for s in parse_igra(f"por_{stn}/{stn}-data.txt",months=(6,7,8),hour=12):
        try: p,z,t,rh,u,v=prep(s,ELEV[stn])
        except Exception: continue
        if len(p)<30 or z.max()<4000: continue
        ml=parcel_abl(p,z,t)
        if ml<300: continue
        truth=ml_grad_inside(p,z,t,ml)
        try: pc,zc,tc,rhc,uc,vc=degrade(p,z,t,rh,u,v)
        except Exception: continue
        if len(pc)<3: continue
        mlc=parcel_abl(pc,zc,tc)
        nin.append(int(((zc<=mlc)&(zc>=80)).sum()))
        B=ml_grad_surface_to_abl(pc,zc,tc,mlc)
        D=cand_D(pc,zc,tc,mlc)
        if not np.isfinite(truth): continue
        R.append((truth,B,D))
    r=np.array(R); nin=np.array(nin); tr=r[:,0]<=THR
    print(f"\n{name}  n={len(r)}  base rate capable={100*tr.mean():.1f}%")
    print(f"  MODEL LEVELS INSIDE THE MIXED LAYER (5-level column):")
    for k in range(0,4):
        print(f"     {k} level(s): {100*np.mean(nin==k):5.1f}%")
    print(f"  {'metric':<30} {'flagged':>8} {'sens':>5} {'spec':>5} {'J':>5}  corr   n_usable")
    for lbl,c in [("B sfc->parcel top (=0.5/depth)",1),("D fit across in-ML levels",2)]:
        ok=np.isfinite(r[:,c])
        if ok.sum()<20:
            print(f"  {lbl:<30} unusable: only {ok.sum()} columns have >=2 in-ML levels"); continue
        pred=r[ok,c]<=THR; t2=tr[ok]
        se,sp,ac,J=rates(t2,pred); cc=np.corrcoef(r[ok,0],r[ok,c])[0,1]
        print(f"  {lbl:<30} {100*pred.mean():7.1f}% {se:5.2f} {sp:5.2f} {J:5.2f}  {cc:+.2f}   {ok.sum()}")
