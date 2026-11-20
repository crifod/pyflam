"""Exercise the WHOLE pyflam atmosphere pipeline on the real GRAF-2025 environment
sondes (Zenodo 15264835) and sanity-check the physics on independent field data.

For each ambient sounding this runs, end to end, the same atmosphere functions the
gridded pyroconvection product uses, and asserts physical plausibility:

  theta/theta_v/q  -> lcl_height_bolton_m  -> bulk_richardson_abl_grid
  -> parcel_mixing_depth_grid -> ML dtheta/dz + cap gamma-theta -> shear height
  -> pyroconvection_type (the Castellnou ladder) -> fire_induced_abl_grid
     (forced by mixed_layer_fire_flux for a reference fire).

These are the first tests of the atmosphere module against real high-resolution
soundings (not synthetic fixtures, not GRAF-derived tables).
"""
import glob, os, sys, numpy as np, pandas as pd, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pyflam.atmosphere import (
    theta_kelvin, dewpoint_from_rh, lcl_height_bolton_m, specific_humidity_from_rh,
    virtual_potential_temperature, bulk_richardson_abl_grid, parcel_mixing_depth_grid,
    shear_height_window, fire_induced_abl_grid, mixed_layer_fire_flux,
    pyroconvection_type, PYROCONVECTION_TYPES)

DATA = os.environ.get("INPLUME_DIR", "/private/tmp/claude-501/-Users-cristianofoderi--softEST-firelab-flammap6-install-0828-2025-pyflam/a1f9bab7-f43b-4aae-86d4-35a99443b113/scratchpad/zenodo_15264835")
REF_FLI_KW = 1.0e4          # reference 10 MW/m fire for the fireABL forcing
G1 = lambda a: np.asarray(a, float)[:, None, None]        # 1-D profile -> (nlev,1,1) grid


def load(path):
    d = pd.read_csv(path); d.columns = [c.strip() for c in d.columns]
    z = pd.to_numeric(d["Altitude (m AGL)"], errors="coerce")
    T = pd.to_numeric(d["Temperature (C)"], errors="coerce")
    RH = pd.to_numeric(d["Relative humidity (%)"], errors="coerce").clip(1, 100)
    P = pd.to_numeric(d["Pressure (Pascal)"], errors="coerce")
    s = pd.to_numeric(d["Speed (m/s)"], errors="coerce")
    hd = pd.to_numeric(d["Heading (degrees)"], errors="coerce")
    ok = np.isfinite(z) & np.isfinite(T) & np.isfinite(P) & (P > 1e4) & (z >= 0)
    z, T, RH, P, s, hd = (v[ok].to_numpy() for v in (z, T, RH, P, s, hd))
    o = np.argsort(z)
    z, T, RH, P, s, hd = z[o], T[o], RH[o], P[o], s[o], hd[o]
    u = np.where(np.isfinite(s), s * np.sin(np.deg2rad(hd)), 0.0)
    v = np.where(np.isfinite(s), s * np.cos(np.deg2rad(hd)), 0.0)
    th = theta_kelvin(T, P / 100.0)
    q = specific_humidity_from_rh(RH, T + 273.15, P)
    thv = virtual_potential_temperature(th, q)
    return dict(z=z, T=T, RH=RH, P=P, u=u, v=v, th=th, thv=thv, q=q)


def layer_grad(z, th, zlo, zhi):
    m = (z >= zlo) & (z <= zhi)
    if m.sum() < 2:
        return np.nan
    return float(np.polyfit(z[m], th[m], 1)[0])


def run(p, fire):
    z, th, thv, u, v = p["z"], p["th"], p["thv"], p["u"], p["v"]
    sfc = z <= z.min() + 120
    Ts, RHs, Ps = p["T"][sfc].mean(), p["RH"][sfc].mean(), p["P"][sfc].mean()
    Tds = dewpoint_from_rh(Ts, RHs)
    # --- the pipeline ---
    lcl = float(lcl_height_bolton_m(Ts + 273.15, Tds + 273.15, Ps))
    thv_s = thv[sfc].mean()
    abl = float(bulk_richardson_abl_grid(G1(z), G1(thv), G1(u), G1(v),
                theta_v_surface=np.array([[thv_s]]), wind_u_surface=np.array([[u[sfc].mean()]]),
                wind_v_surface=np.array([[v[sfc].mean()]]))[0, 0])
    parcel = float(parcel_mixing_depth_grid(G1(z), G1(th), np.array([[th[sfc].mean()]]))[0, 0])
    ml_grad = layer_grad(z, th, z.min(), parcel)                       # ML stability
    cap = layer_grad(z, th, abl + 200, abl + 1200)                     # free-trop cap
    shear_h = float(shear_height_window(z, u, v))
    rh_top = float(p["RH"][(z >= abl - 150) & (z <= abl + 150)].mean()) if np.isfinite(abl) else np.nan
    th_below = th[z <= abl].mean()
    q_ref = mixed_layer_fire_flux(REF_FLI_KW * 1000.0, abl)            # W/m^2
    fabl = float(fire_induced_abl_grid(G1(z), G1(th), theta_mean_below=np.array([[th_below]]),
                 heat_flux=np.array([[q_ref]]), blh=np.array([[abl]]))[0, 0])
    shear_dist = (min(abs(shear_h - abl), abs(shear_h - lcl)) / abl
                  if np.isfinite(shear_h) and abl > 0 else np.nan)
    try:
        cls = pyroconvection_type(lcl_abl_ratio=lcl / abl, ml_theta_gradient=ml_grad,
                gamma_theta=cap if np.isfinite(cap) else 5e-3, rh_top_abl=rh_top,
                shear_distance=shear_dist if np.isfinite(shear_dist) else 1.0,
                fireline_intensity_kw=REF_FLI_KW, ladder="adaptive")
        cls = cls[0] if isinstance(cls, tuple) else cls
    except Exception as e:
        cls = f"ERR:{str(e)[:20]}"
    # --- physical sanity ---
    checks = {
        "theta_v>=theta": bool(np.all(thv[np.isfinite(thv)] >= th[np.isfinite(thv)] - 1e-6)),
        "ABL in[150,4500]": bool(150 <= abl <= 4500),
        "LCL>0": bool(lcl > 0),
        "cap>=0 (stable aloft)": bool(np.isnan(cap) or cap >= -1e-4),
        "fireABL>=ABL": bool(fabl >= abl - 1e-6),
        "class valid": cls in PYROCONVECTION_TYPES,
    }
    return dict(fire=fire, ztop=z.max(), Ts=Ts, RHs=RHs, LCL=lcl, ABL=abl, parcel=parcel,
                mlgrad=ml_grad, cap=cap, shear=shear_h, rh_top=rh_top, fireABL=fabl,
                lcl_abl=lcl/abl, cls=cls, ok=all(checks.values()), checks=checks)


def main():
    envs = sorted(glob.glob(f"{DATA}/*_Environment*.raw_flight_history.csv"))
    print(f"{'fire':16} {'ztop':>5} {'Ts':>4} {'RH':>3} {'LCL':>5} {'ABL':>5} {'parcel':>6} "
          f"{'MLgrad':>8} {'cap':>8} {'LCL/ABL':>7} {'fABL':>5} {'class':>18} {'phys'}")
    rows = []
    for e in envs:
        fire = os.path.basename(e).split("_")[0]
        try:
            r = run(load(e), fire); rows.append(r)
            print(f"{fire:16} {r['ztop']:5.0f} {r['Ts']:4.0f} {r['RHs']:3.0f} {r['LCL']:5.0f} "
                  f"{r['ABL']:5.0f} {r['parcel']:6.0f} {r['mlgrad']:8.1e} {r['cap']:8.1e} "
                  f"{r['lcl_abl']:7.1f} {r['fireABL']:5.0f} {r['cls']:>18} {'OK' if r['ok'] else 'FAIL'}")
        except Exception as ex:
            print(f"{fire:16} ERROR {str(ex)[:50]}")
    print(f"\n{sum(r['ok'] for r in rows)}/{len(rows)} sondes pass all physical-sanity checks.")
    # which checks failed anywhere
    allchecks = {}
    for r in rows:
        for k, val in r["checks"].items():
            allchecks.setdefault(k, []).append(val)
    print("physical checks (pass/total):")
    for k, vals in allchecks.items():
        print(f"  {k:22} {sum(vals)}/{len(vals)}")
    print("\nPipeline functions exercised on real sondes: theta/theta_v/q, lcl_height_bolton_m,")
    print("bulk_richardson_abl_grid, parcel_mixing_depth_grid, ML+cap dtheta/dz, shear_height_window,")
    print("pyroconvection_type (Castellnou ladder), mixed_layer_fire_flux, fire_induced_abl_grid.")


if __name__ == "__main__":
    main()
