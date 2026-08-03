"""Resolution-degradation experiment.

Same radiosonde column, twice: at full sounding resolution (ground truth) and
subsampled onto ICON-2I's 5 usable pressure levels. Isolates the vertical-resolution
bias in (a) the bulk-Richardson ABL and (b) the mixed-layer dtheta/dz, with model
error held out entirely. Uses pyflam's shipped functions, so it tests the real code.

IGRA v2 note: only mandatory levels carry a geopotential height, and moisture is
reported as a dewpoint depression, not RH. Heights are therefore reconstructed
hypsometrically from p and virtual temperature (anchored at the station elevation),
which keeps every significant level -- discarding them would destroy exactly the
vertical resolution this experiment is about.
"""
import sys, numpy as np

from pyflam.atmosphere import (
    bulk_richardson_abl_grid, theta_kelvin, specific_humidity_from_rh,
    virtual_potential_temperature, saturation_vapour_pressure_pa,
    DEFAULT_PYROCONV_THRESHOLDS, _G, _RD)

ICON_LEVELS = [1000, 925, 850, 700, 500]
MISS = (-9999, -8888, -999)
ELEV = {"ITM00016245": 12.0, "ITM00016144": 11.0}     # station elevation (m)


def _num(s):
    try:
        v = int(s)
    except ValueError:
        return np.nan
    return np.nan if v in MISS else float(v)


def parse_igra(path, months=(5, 6, 7, 8, 9), hour=12, min_levels=30):
    out, cur = [], None
    for line in open(path):
        if line.startswith("#"):
            if cur is not None and len(cur["p"]) >= min_levels:
                out.append(cur)
            mo, hr = int(line[18:20]), int(line[24:26])
            cur = dict(date=f"{line[13:17]}-{line[18:20]}-{line[21:23]} {hr:02d}Z",
                       p=[], t=[], dpdp=[], wd=[], ws=[]) \
                if (mo in months and hr == hour) else None
            continue
        if cur is None:
            continue
        p, t = _num(line[9:15]), _num(line[22:27])
        if not np.isfinite(p) or not np.isfinite(t):
            continue
        cur["p"].append(p / 100.0)                       # Pa -> hPa
        cur["t"].append(t / 10.0 + 273.15)              # tenths C -> K
        cur["dpdp"].append(_num(line[34:39]) / 10.0)
        cur["wd"].append(_num(line[40:45]))
        cur["ws"].append(_num(line[46:51]) / 10.0)
    if cur is not None and len(cur["p"]) >= min_levels:
        out.append(cur)
    return out


def _fill(x, p):
    """Interpolate missing values against log-pressure (monotone coordinate)."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    if ok.sum() < 2:
        return np.full_like(x, np.nan)
    lp = np.log(p)
    return np.interp(lp, lp[ok][::-1], x[ok][::-1]) if lp[0] > lp[-1] else \
        np.interp(lp, lp[ok], x[ok])


def prep(s, elev):
    p = np.asarray(s["p"], float)
    o = np.argsort(-p)                                   # descending p = ascending height
    p = p[o]
    t = np.asarray(s["t"], float)[o]
    dpdp = _fill(np.asarray(s["dpdp"], float)[o], p)
    wd = _fill(np.asarray(s["wd"], float)[o], p)
    ws = _fill(np.asarray(s["ws"], float)[o], p)

    td = t - np.clip(dpdp, 0, None)
    rh = np.clip(100.0 * saturation_vapour_pressure_pa(td)
                 / saturation_vapour_pressure_pa(t), 1.0, 100.0)
    q = specific_humidity_from_rh(rh, t, p * 100.0)
    tv = t * (1.0 + 0.61 * q)

    # Hypsometric integration from the station elevation.
    z = np.empty_like(p)
    z[0] = elev
    for i in range(1, len(p)):
        tvb = 0.5 * (tv[i] + tv[i - 1])
        z[i] = z[i - 1] + _RD * tvb / _G * np.log(p[i - 1] / p[i])

    u = -ws * np.sin(np.deg2rad(wd))
    v = -ws * np.cos(np.deg2rad(wd))
    good = np.isfinite(p) & np.isfinite(t) & np.isfinite(rh) & np.isfinite(u) & np.isfinite(v)
    p, z, t, rh, u, v = p[good], z[good], t[good], rh[good], u[good], v[good]
    return p, z - z[0], t, rh, u, v                      # z now AGL


def abl_rib(p, z, t, rh, u, v):
    """pyflam's shipped Rib ABL, driven as a single column."""
    th = theta_kelvin(t - 273.15, p)
    q = specific_humidity_from_rh(rh, t, p * 100.0)
    thv = virtual_potential_temperature(th, q)
    g = lambda a: np.asarray(a, float)[:, None, None]
    return float(bulk_richardson_abl_grid(
        g(z), g(thv), g(u), g(v),
        theta_v_surface=np.array([[thv[0]]]),
        wind_u_surface=np.array([[u[0]]]), wind_v_surface=np.array([[v[0]]]))[0, 0])


def parcel_abl(p, z, t, excess=0.5):
    """Holzworth parcel method: first height where theta exceeds surface theta + excess.
    The classic CBL mixing depth -- an independent check on the Rib estimate."""
    th = theta_kelvin(t - 273.15, p)
    tgt = th[0] + excess
    for i in range(1, len(z)):
        if th[i] >= tgt:
            f = (tgt - th[i - 1]) / max(th[i] - th[i - 1], 1e-9)
            return float(z[i - 1] + f * (z[i] - z[i - 1]))
    return float(z[-1])


def ml_grad_surface_to_abl(p, z, t, abl):
    th = theta_kelvin(t - 273.15, p)
    return (float(np.interp(abl, z, th)) - th[0]) / max(abl - 2.0, 1.0)


def ml_grad_inside(p, z, t, top):
    """dtheta/dz *within* the mixed layer (least-squares, 50 m -> top).
    A well-mixed CBL should give ~0. This is the physical target."""
    th = theta_kelvin(t - 273.15, p)
    m = (z >= 50) & (z <= top)
    return float(np.polyfit(z[m], th[m], 1)[0]) if m.sum() >= 3 else np.nan


def degrade(p, z, t, rh, u, v, levels=ICON_LEVELS):
    lp = np.log(p)
    tl = np.array(sorted([L for L in levels if p.min() <= L <= p.max()], reverse=True), float)
    x = np.log(tl)
    o = np.argsort(lp)
    f = lambda a: np.interp(x, lp[o], np.asarray(a, float)[o])
    return tl, f(z), f(t), np.clip(f(rh), 1, 100), f(u), f(v)


if __name__ == "__main__":
    thr = DEFAULT_PYROCONV_THRESHOLDS.ml_stable
    for stn, name in [("ITM00016245", "Pratica di Mare (Rome)"),
                      ("ITM00016144", "S. Pietro Capofiume (Po)")]:
        snds = parse_igra(f"{stn}-data.txt")
        rows = []
        print(f"\n{'='*100}\n{name} -- {len(snds)} warm-season 12Z soundings\n{'='*100}")
        print(f"{'date':<15} {'ABL_rib':>7} {'ABL_crs':>7} {'ABL_par':>7} | "
              f"{'MLg_IN':>9} {'MLg_full':>9} {'MLg_crs':>9} | {'lev':>4}")
        for s in snds:
            try:
                p, z, t, rh, u, v = prep(s, ELEV[stn])
            except Exception:
                continue
            if len(p) < 30 or z.max() < 4000:
                continue
            a_full, a_par = abl_rib(p, z, t, rh, u, v), parcel_abl(p, z, t)
            pc, zc, tc, rhc, uc, vc = degrade(p, z, t, rh, u, v)
            if len(pc) < 3:
                continue
            a_crs = abl_rib(pc, zc, tc, rhc, uc, vc)
            g_in = ml_grad_inside(p, z, t, a_par)
            g_full = ml_grad_surface_to_abl(p, z, t, a_full)
            g_crs = ml_grad_surface_to_abl(pc, zc, tc, a_crs)
            rows.append((a_full, a_crs, a_par, g_in, g_full, g_crs))
            print(f"{s['date']:<15} {a_full:7.0f} {a_crs:7.0f} {a_par:7.0f} | "
                  f"{g_in:9.2e} {g_full:9.2e} {g_crs:9.2e} | {len(p):4d}")
        if not rows:
            continue
        r = np.array(rows)
        med = lambda c: np.nanmedian(r[:, c])
        print("-" * 100)
        print(f"{'MEDIAN':<15} {med(0):7.0f} {med(1):7.0f} {med(2):7.0f} | "
              f"{med(3):9.2e} {med(4):9.2e} {med(5):9.2e} |  n={len(r)}")
        print(f"\n  ABL bias, coarse - full      : {np.nanmedian(r[:,1]-r[:,0]):+7.0f} m")
        print(f"  ABL, Rib(full) - parcel      : {np.nanmedian(r[:,0]-r[:,2]):+7.0f} m")
        print(f"  ML dtheta/dz INSIDE ML (truth): {med(3):.2e} K/m   <- should be ~0")
        print(f"  ML dtheta/dz sfc->ABL, full   : {med(4):.2e} K/m")
        print(f"  ML dtheta/dz sfc->ABL, 5 lev  : {med(5):.2e} K/m")
        print(f"  pyroCu-capable fraction (dtheta/dz <= {thr:.0e}):")
        for lbl, c in [("inside-ML (truth)", 3), ("sfc->ABL full-res", 4), ("sfc->ABL ICON 5-lev", 5)]:
            print(f"    {lbl:22s}: {100*np.nanmean(r[:,c] <= thr):5.1f}%")
