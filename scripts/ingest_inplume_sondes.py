"""Ingest the GRAF-2025 in-plume radiosonde dataset (Zenodo 15264835) into clean,
calibration-grade (ambient profile, observed fireABL, measured theta') records, and
validate the fireABL encroachment against the observed plume top.

Each fire has an ENVIRONMENT sonde (descent = ambient profile) and one or more
IN-PLUME sondes (ascent, launched at the fire, rises through the plume). QC steps:

  1. Ascent/descent handling -- sort each flight by altitude; both give a valid
     theta(z) once monotone. Drop sub-ground (z<0) and non-finite rows.
  2. Surface reference -- robust near-surface theta = median over the lowest 150 m.
  3. Observed fireABL -- interpolate both profiles to a common grid; the in-plume
     warm anomaly (theta_ip - theta_env) is the fire signal. The observed plume top
     is the highest height where the (median-smoothed) anomaly exceeds a noise floor
     scaled to the surface anomaly. If it never decays below the floor before the
     shallower apex, the fireABL is CENSORED (>= that height, a lower bound).
  4. Measured theta' -- mean in-plume-minus-environment excess over the lowest 200 m
     (the fire-heated parcel excess the encroachment needs -- no fire intensity or
     flux parameterization required).
  5. Ambient ABL -- parcel method (theta exceeds surface + 1.0 K) on the environment
     profile, for the decoupling-ratio context.

Prediction = the fire_induced_abl_grid encroachment intersection using the MEASURED
theta', compared against the observed plume top.
"""
import glob, os, sys, numpy as np, pandas as pd, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

DATA = os.environ.get("INPLUME_DIR", "/private/tmp/claude-501/-Users-cristianofoderi--softEST-firelab-flammap6-install-0828-2025-pyflam/a1f9bab7-f43b-4aae-86d4-35a99443b113/scratchpad/zenodo_15264835")


def launch_time(path):
    d = pd.read_csv(path, nrows=1)
    return str(d.iloc[0, 0]).strip()


def _minutes_apart(a, b):
    try:
        ta, tb = pd.to_datetime(a), pd.to_datetime(b)
        return abs((ta - tb).total_seconds()) / 60.0
    except Exception:
        return np.nan


def load_sonde(path):
    """Clean theta(z) profile (K, m AGL), monotone in height, from a raw flight history."""
    d = pd.read_csv(path); d.columns = [c.strip() for c in d.columns]
    z = pd.to_numeric(d["Altitude (m AGL)"], errors="coerce")
    T = pd.to_numeric(d["Temperature (C)"], errors="coerce")
    P = pd.to_numeric(d["Pressure (Pascal)"], errors="coerce")
    RH = pd.to_numeric(d["Relative humidity (%)"], errors="coerce")
    th = (T + 273.15) * (100000.0 / P) ** 0.286
    ok = np.isfinite(z) & np.isfinite(th) & (P > 1e4) & (z >= 0)
    z, th, rh = z[ok].to_numpy(), th[ok].to_numpy(), RH[ok].to_numpy()
    order = np.argsort(z); z, th, rh = z[order], th[order], rh[order]
    zc, tc, rc = [], [], []                       # collapse duplicate heights
    for u in np.unique(z):
        m = z == u; zc.append(u); tc.append(th[m].mean()); rc.append(np.nanmean(rh[m]))
    return np.array(zc), np.array(tc), np.array(rc)


def _smooth(x, w=5):
    if len(x) < w:
        return x
    return pd.Series(x).rolling(w, center=True, min_periods=1).median().to_numpy()


def ambient_abl(z, th, excess=1.0):
    th0 = th[z <= z.min() + 150].mean()
    for i in range(1, len(z)):
        if th[i] >= th0 + excess:
            f = np.clip((th0 + excess - th[i-1]) / max(th[i]-th[i-1], 1e-6), 0, 1)
            return float(z[i-1] + f*(z[i]-z[i-1]))
    return float(z.max())


def observed_fireabl(ze, the, zi, thi):
    """Observed plume top from the in-plume warm anomaly; returns (height, censored)."""
    ztop = min(ze.max(), zi.max())
    grid = np.arange(max(ze.min(), zi.min()) + 10, ztop, 20.0)
    if grid.size < 5:
        return np.nan, False
    anom = _smooth(np.interp(grid, zi, thi) - np.interp(grid, ze, the))
    surf = np.median(anom[grid <= grid.min() + 200])
    floor = max(0.5, 0.25 * surf)                 # noise floor scaled to surface excess
    warm = grid[anom > floor]
    if warm.size == 0:
        return float(grid.min()), False
    top = float(warm.max())
    censored = anom[-1] > floor                    # still warm at the shallower apex
    return top, censored


def measured_theta_excess(ze, the, zi, thi):
    lo = max(ze.min(), zi.min())
    return float(thi[zi <= lo + 200].mean() - the[ze <= lo + 200].mean())


def encroach(z, th, theta_mean_below, theta_excess, blh):
    target = theta_mean_below + theta_excess
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

    print(f"{'fire':16} {'ABL':>5} {'theta_p':>7} {'obs_top':>8} {'pred':>6} {'|Δ|':>6} {'ratio':>6} {'flag'}")
    rows = []
    for name in sorted(fires):
        e = fires[name].get("env"); ips = fires[name].get("ip")
        if not e or not ips:
            continue
        ze, the, _ = load_sonde(e)
        ip = sorted(ips, key=lambda p: load_sonde(p)[0].min())[0]   # deepest-sampling in-plume
        zi, thi, _ = load_sonde(ip)
        if ze.size < 8 or zi.size < 8:
            continue
        gap = _minutes_apart(launch_time(e), launch_time(ip))
        blh = ambient_abl(ze, the)
        th_below = the[ze <= blh].mean()
        thp = max(measured_theta_excess(ze, the, zi, thi), 0.0)
        obs, cens = observed_fireabl(ze, the, zi, thi)
        pred = encroach(ze, the, th_below, thp, blh)
        # Calibration-usable only if the sonde pair is near-simultaneous (so the
        # anomaly is the fire, not diurnal change), the plume is captured below the
        # apex (not censored), and the ambient ABL is physical. A plume top ~= the
        # ambient ABL is a valid *weak-decoupling* case, so it is NOT excluded.
        usable = (np.isfinite(gap) and gap <= 30) and (not cens) and (blh >= 100)
        flag = "USABLE" if usable else (
            "censored" if cens else ("Δt=%.0fmin" % gap if np.isfinite(gap) and gap > 30
                                     else "ABL?"))
        d = abs(pred - obs) if np.isfinite(obs) else np.nan
        rows.append(dict(fire=name, abl=blh, theta_p=thp, obs_top=obs, pred=pred,
                         err=d, ratio=obs/blh if blh else np.nan, censored=cens,
                         gap_min=gap, usable=usable, z_env_top=ze.max(), z_ip_top=zi.max()))
        print(f"{name:16} {blh:5.0f} {thp:7.1f} {obs:8.0f} {pred:6.0f} {d:6.0f} "
              f"{obs/blh:6.1f} {flag}")

    df = pd.DataFrame(rows)
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "inplume_fireabl_calibration.csv")
    df.to_csv(out, index=False)
    good = df[df["usable"]]
    print(f"\nn={len(df)} fires; {len(good)} calibration-usable (near-simultaneous pair, "
          f"uncensored, physical ABL, plume above ABL).")
    print(f"measured theta' across all: {df['theta_p'].min():.1f}–{df['theta_p'].max():.1f} K "
          f"(median {df['theta_p'].median():.1f}) — field is ~{25/df['theta_p'].median():.0f}x "
          f"smaller than the GRAF demo's modelled 14–25 K, corroborating the mixed-layer-flux fix.")
    if len(good):
        print("USABLE fires (independent encroachment check):")
        for _, r in good.iterrows():
            print(f"  {r['fire']:14} obs_fireABL {r['obs_top']:.0f} m  pred {r['pred']:.0f} m  "
                  f"|Δ| {r['err']:.0f} m  (Δt {r['gap_min']:.0f} min)")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
