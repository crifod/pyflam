"""Ingest the GRAF-2025 in-plume radiosonde dataset (Zenodo 15264835) into a
canonical, quality-controlled fireABL calibration table.

Each fire has an ENVIRONMENT sonde (descent = ambient profile) and one or more
IN-PLUME sondes (ascent, launched at the fire, rising through the plume). Every
(fire, in-plume sonde) pair is one calibration observation, forced by the MEASURED
theta' excess (in-plume minus environment near the surface) -- so no fire intensity
or flux parameterisation is needed. QC steps:

  1. Ascent/descent -- sort each flight by altitude (both give valid theta(z)).
  2. Ambient diagnostics from the environment sonde: bulk-Richardson ABL, Bolton LCL,
     free-troposphere cap gamma-theta.
  3. Observed fireABL -- highest height where the median-smoothed in-plume warm
     anomaly (theta_ip - theta_env) exceeds a noise floor; CENSORED (lower bound) if
     still warm at the shallower apex.
  4. Measured theta' -- mean in-plume-minus-environment excess over the lowest 200 m.
  5. Predicted fireABL -- fire_induced_abl_grid encroachment with the measured theta'.
  6. Usability -- near-simultaneous pair (<=30 min), uncensored, physical ABL.

Writes docs/inplume_fireabl_calibration.csv (all observations + quality flags).
"""
import glob, os, sys, numpy as np, pandas as pd, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pyflam.atmosphere import (
    theta_kelvin, dewpoint_from_rh, lcl_height_bolton_m, specific_humidity_from_rh,
    virtual_potential_temperature, bulk_richardson_abl_grid, fire_induced_abl_grid)

DATA = os.environ.get("INPLUME_DIR", "/private/tmp/claude-501/-Users-cristianofoderi--softEST-firelab-flammap6-install-0828-2025-pyflam/a1f9bab7-f43b-4aae-86d4-35a99443b113/scratchpad/zenodo_15264835")
_G1 = lambda a: np.asarray(a, float)[:, None, None]
# Martorell already has ERA5+sonde+fireABL in the existing GRAF set; the 5 new fires:
NEW_FIRES = {"CasablancaIII08", "CasablancaIII10", "GranjaEscarp", "Rojals",
             "SantaAna", "SelvadelCamp"}


def launch_time(path):
    return str(pd.read_csv(path, nrows=1).iloc[0, 0]).strip()


def _minutes_apart(a, b):
    try:
        return abs((pd.to_datetime(a) - pd.to_datetime(b)).total_seconds()) / 60.0
    except Exception:
        return np.nan


def load(path):
    d = pd.read_csv(path); d.columns = [c.strip() for c in d.columns]
    z = pd.to_numeric(d["Altitude (m AGL)"], errors="coerce")
    T = pd.to_numeric(d["Temperature (C)"], errors="coerce")
    P = pd.to_numeric(d["Pressure (Pascal)"], errors="coerce")
    RH = pd.to_numeric(d["Relative humidity (%)"], errors="coerce").clip(1, 100)
    s = pd.to_numeric(d["Speed (m/s)"], errors="coerce")
    hd = pd.to_numeric(d["Heading (degrees)"], errors="coerce")
    ok = np.isfinite(z) & np.isfinite(T) & np.isfinite(P) & (P > 1e4) & (z >= 0)
    z, T, P, RH, s, hd = (v[ok].to_numpy() for v in (z, T, P, RH, s, hd))
    o = np.argsort(z); z, T, P, RH, s, hd = z[o], T[o], P[o], RH[o], s[o], hd[o]
    zc, Tc, Pc, RHc, uc, vc = [], [], [], [], [], []
    u = np.where(np.isfinite(s), s*np.sin(np.deg2rad(hd)), 0.0)
    v = np.where(np.isfinite(s), s*np.cos(np.deg2rad(hd)), 0.0)
    for uu in np.unique(z):
        m = z == uu
        zc.append(uu); Tc.append(T[m].mean()); Pc.append(P[m].mean())
        RHc.append(np.nanmean(RH[m])); uc.append(u[m].mean()); vc.append(v[m].mean())
    z, T, P, RH, u, v = map(np.array, (zc, Tc, Pc, RHc, uc, vc))
    th = theta_kelvin(T, P/100.0)
    q = specific_humidity_from_rh(RH, T+273.15, P)
    thv = virtual_potential_temperature(th, q)
    return dict(z=z, T=T, P=P, RH=RH, th=th, thv=thv, u=u, v=v)


def _smooth(x, w=5):
    return pd.Series(x).rolling(w, center=True, min_periods=1).median().to_numpy() if len(x) >= w else x


def ambient(p):
    z, th, thv, u, v = p["z"], p["th"], p["thv"], p["u"], p["v"]
    sfc = z <= z.min() + 150
    Ts, RHs, Ps = p["T"][sfc].mean(), p["RH"][sfc].mean(), p["P"][sfc].mean()
    lcl = float(lcl_height_bolton_m(Ts+273.15, dewpoint_from_rh(Ts, RHs)+273.15, Ps))
    abl = float(bulk_richardson_abl_grid(
        _G1(z), _G1(thv), _G1(u), _G1(v), theta_v_surface=np.array([[thv[sfc].mean()]]),
        wind_u_surface=np.array([[u[sfc].mean()]]), wind_v_surface=np.array([[v[sfc].mean()]]))[0, 0])
    cm = (z >= abl+200) & (z <= abl+1200)
    cap = float(np.polyfit(z[cm], th[cm], 1)[0]) if cm.sum() >= 2 else np.nan
    return lcl, abl, cap


def observed_fireabl(ze, the, zi, thi):
    ztop = min(ze.max(), zi.max())
    grid = np.arange(max(ze.min(), zi.min())+10, ztop, 20.0)
    if grid.size < 5:
        return np.nan, False
    anom = _smooth(np.interp(grid, zi, thi) - np.interp(grid, ze, the))
    floor = max(0.5, 0.25*np.median(anom[grid <= grid.min()+200]))
    warm = grid[anom > floor]
    if warm.size == 0:
        return float(grid.min()), False
    return float(warm.max()), bool(anom[-1] > floor)


def encroach(z, th, theta_mean_below, thp, blh):
    target = theta_mean_below + thp
    for i in range(1, len(z)):
        if z[i] >= blh and th[i] >= target:
            f = np.clip((target - th[i-1]) / max(th[i]-th[i-1], 1e-6), 0, 1)
            return float(max(z[i-1] + f*(z[i]-z[i-1]), blh))
    return float(z.max())


def main():
    fires = {}
    for f in glob.glob(f"{DATA}/*_Environment*.raw_flight_history.csv"):
        fires.setdefault(os.path.basename(f).split("_")[0], {})["env"] = f
    for f in glob.glob(f"{DATA}/*_In-plume*.raw_flight_history.csv"):
        fires.setdefault(os.path.basename(f).split("_")[0], {}).setdefault("ip", []).append(f)

    rows = []
    for name in sorted(fires):
        e = fires[name].get("env"); ips = fires[name].get("ip")
        if not e or not ips:
            continue
        pe = load(e); ze, the = pe["z"], pe["th"]
        if ze.size < 8:
            continue
        lcl, abl, cap = ambient(pe)
        th_below = the[ze <= abl].mean()
        for ip in sorted(ips):
            pi = load(ip); zi, thi = pi["z"], pi["th"]
            if zi.size < 8:
                continue
            gap = _minutes_apart(launch_time(e), launch_time(ip))
            lo = max(ze.min(), zi.min())
            thp = max(thi[zi <= lo+200].mean() - the[ze <= lo+200].mean(), 0.0)
            obs, cens = observed_fireabl(ze, the, zi, thi)
            pred = encroach(ze, the, th_below, thp, abl)
            usable = (np.isfinite(gap) and gap <= 30) and (not cens) and (abl >= 100)
            rows.append(dict(
                fire=name, new_fire=name in NEW_FIRES,
                sonde=os.path.basename(ip).replace(".raw_flight_history.csv", ""),
                ambient_ABL_m=round(abl), LCL_m=round(lcl),
                cap_gamma_Km=round(cap, 6) if np.isfinite(cap) else np.nan,
                theta_excess_K=round(thp, 2), obs_fireABL_m=round(obs),
                pred_fireABL_m=round(pred), abs_err_m=round(abs(pred-obs)),
                decoupling_ratio=round(obs/abl, 2) if abl else np.nan,
                dt_env_min=round(gap) if np.isfinite(gap) else np.nan,
                censored=cens, calibration_usable=usable))

    df = pd.DataFrame(rows).sort_values(["new_fire", "fire", "sonde"], ascending=[False, True, True])
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "inplume_fireabl_calibration.csv")
    df.to_csv(out, index=False)

    new = df[df["new_fire"]]
    good = new[new["calibration_usable"]]
    print(f"Ingested {len(df)} in-plume observations across {df['fire'].nunique()} fires "
          f"({new['fire'].nunique()} new).")
    print(f"\n{'sonde':34} {'ABL':>5} {'LCL':>5} {'θ_exc':>6} {'obs':>5} {'pred':>5} {'|Δ|':>5} "
          f"{'ratio':>5} {'Δt':>4} {'use'}")
    for _, r in new.iterrows():
        u = "USE" if r["calibration_usable"] else ("cens" if r["censored"] else f"{r['dt_env_min']:.0f}m")
        print(f"{r['sonde']:34} {r['ambient_ABL_m']:5.0f} {r['LCL_m']:5.0f} {r['theta_excess_K']:6.1f} "
              f"{r['obs_fireABL_m']:5.0f} {r['pred_fireABL_m']:5.0f} {r['abs_err_m']:5.0f} "
              f"{r['decoupling_ratio']:5.1f} {r['dt_env_min'] if np.isfinite(r['dt_env_min']) else -1:4.0f} {u}")
    print(f"\n5 new fires: {len(new)} observations, {len(good)} calibration-usable "
          f"(near-simultaneous + uncensored).")
    print(f"measured theta' (new fires): {new['theta_excess_K'].min():.1f}–{new['theta_excess_K'].max():.1f} K "
          f"(median {new['theta_excess_K'].median():.1f}) — vs GRAF demo 14–25 K.")
    if len(good):
        print(f"usable encroachment |pred-obs|: median {good['abs_err_m'].median():.0f} m")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
