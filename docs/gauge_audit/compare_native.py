"""Old (every cell at 1000 hPa) vs new (surface pressure) classifier on the same real ICON-2I field."""
import os, sys, importlib.util, numpy as np
R = sys.argv[1]
spec = importlib.util.spec_from_file_location("m", R + "/tests/pyroconv_icon2i_native.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
from pyflam.atmosphere import (theta_kelvin, lcl_height_m, pyroconvection_type,
                               PYROCONVECTION_TYPE_LEVEL, relative_humidity_from_dewpoint)

def legacy(T2m, RH2m, Tlev):                 # verbatim pre-fix body of classify()
    nlat, nlon = T2m.shape; out = np.zeros((nlat, nlon), np.int16); z = m.PLEV_BARO
    th = theta_kelvin(Tlev, np.array(m.LEVELS)[:, None, None])
    i850, i700, i500 = m.LEVELS.index(850), m.LEVELS.index(700), m.LEVELS.index(500)
    th_sfc = theta_kelvin(T2m, 1000.0)
    ml_grad = (th[i850] - th_sfc) / (z[i850] - 110.0)
    gamma = (th[i500] - th[i700]) / (z[i500] - z[i700])
    lcl = lcl_height_m(T2m, RH2m)
    abl = np.full((nlat, nlon), z[-1])
    for k in range(1, len(z)):
        crossed = (th[k] > th_sfc + 0.5) & (abl == z[-1])
        denom = np.where(th[k] != th[k-1], th[k] - th[k-1], 1e9)
        frac = np.clip((th_sfc + 0.5 - th[k-1]) / denom, 0.0, 1.0)
        abl = np.where(crossed, z[k-1] + frac*(z[k]-z[k-1]), abl)
    abl = np.maximum(abl, 100.0); ratio = lcl/abl
    for a in range(nlat):
        for b in range(nlon):
            if abl[a, b] < 600.0: continue
            out[a, b] = PYROCONVECTION_TYPE_LEVEL[pyroconvection_type(
                lcl_abl_ratio=float(ratio[a, b]), ml_theta_gradient=float(ml_grad[a, b]),
                gamma_theta=float(gamma[a, b]))]
    return out, ml_grad

d = m.load()
RH = relative_humidity_from_dewpoint(d["T2m"], d["Td2m"])
bins = [(1100, 975), (975, 925), (925, 875), (875, 0)]
for hour in (0, 12):
    si = d["idx"][m.HOURS.index(hour)]
    T2m, Tlev, PS = d["T2m"][si], d["Tlev"][si], d["PS"][si]
    old, g_old = legacy(T2m, RH[si], Tlev)
    new = m.classify(T2m, RH[si], Tlev, PS)
    print(f"--- {hour:02d}Z  (classes 0 surface, 1 conv plume, 2 overshoot, 3 resilient, 4 deep; -1 nodata)")
    print(f"{'surface p':>12} {'cells':>6} {'old mean':>9} {'new mean':>9} {'changed':>8} {'old ML grad K/m':>16}")
    for hi_, lo_ in bins:
        sel = (PS <= hi_) & (PS > lo_)
        ok = sel & (new >= 0)
        if not sel.any(): continue
        print(f"{lo_:>5}-{hi_:<5} hPa {sel.sum():6d} {old[sel].mean():9.2f} "
              f"{(new[ok].mean() if ok.any() else float('nan')):9.2f} {(old[ok] != new[ok]).sum():8d} "
              f"{np.nanmedian(g_old[sel]):16.2e}   nodata={int((sel & (new < 0)).sum())}")
    for c in range(5):
        print(f"   class {c}: old {int((old == c).sum()):6d}  new {int((new == c).sum()):6d}")
