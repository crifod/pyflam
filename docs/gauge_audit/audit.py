"""Gauge-invariance audit of pyflam's PFT: perturb reference constants, measure response."""
import sys, numpy as np
sys.path.insert(0, "tests"); sys.path.insert(0, "src")
import pyflam.atmosphere as A
from test_atmosphere import _synthetic_column

orig = dict(theta=A.theta_kelvin, ml=A._mixed_layer_state, sp=A._saturation_point)

def install(p0_pa=A._P0_PA, kappa=None, kappa_fwd=A._KAPPA_DRY, dtheta_mode="theta"):
    """Re-gauge: theta reference pressure p0 and Poisson exponent (forward and inverse)."""
    kinv = A._KAPPA_DRY if kappa is None else kappa
    kf = kappa_fwd if kappa is None else kappa
    def theta_kelvin(temp_c, p_hpa):
        return (np.asarray(temp_c, float) + 273.15) * (p0_pa / 100.0 / np.asarray(p_hpa, float)) ** kf
    def ml(z, theta, q, ps, usable, shape, *, iterations=6):
        z_base = np.min(np.where(usable, z, np.inf), axis=0); z_base = np.where(np.isfinite(z_base), z_base, 0.0)
        depth = np.full(shape, 500.0)
        for _ in range(iterations):
            w = usable & (z <= (z_base + depth)[None, ...])
            wt = np.where(w, np.maximum(z, 1.0), 0.0); tot = np.maximum(wt.sum(0), 1e-9)
            th_ml = (wt * np.where(w, theta, 0.0)).sum(0) / tot
            q_ml = (wt * np.where(w, q, 0.0)).sum(0) / tot
            t_sfc = th_ml * (ps / p0_pa) ** kinv
            e = q_ml * ps / (A._EPSILON + (1 - A._EPSILON) * q_ml)
            le = np.log(np.maximum(e, 1.0) / 611.2); td = 243.5 * le / (17.67 - le)
            depth = np.clip(A.lcl_height_bolton_m(t_sfc, td + 273.15, ps), 100.0, 5000.0)
        return th_ml, q_ml, depth
    def sp(beta, th_ml, q_ml, ps, h2m):
        t_ml = th_ml * (ps / p0_pa) ** kinv
        th_p = (1 + beta) * th_ml
        # moisture increment: per K of theta (as coded) or per K of surface temperature (invariant)
        dT = beta * (th_ml if dtheta_mode == "theta" else t_ml)
        q_p = q_ml + dT / h2m * 1e-3
        t_p = th_p * (ps / p0_pa) ** kinv
        e = q_p * ps / (A._EPSILON + (1 - A._EPSILON) * q_p)
        le = np.log(np.maximum(e, 1.0) / 611.2); td = 243.5 * le / (17.67 - le)
        z_sp = A.lcl_height_bolton_m(t_p, td + 273.15, ps)
        t_sp = t_p - A._G / A._CP_DRY * z_sp
        p_sp = ps * np.maximum(t_sp / np.maximum(t_p, 1.0), 1e-6) ** (1 / kinv)
        return z_sp, t_sp, p_sp, t_p
    A.theta_kelvin, A._mixed_layer_state, A._saturation_point = theta_kelvin, ml, sp

def pft(col, **kw):
    out = A.pyrocb_firepower_threshold_grid(col["height_agl_m"], col["temperature_k"], col["pressure_pa"],
            col["spec_humidity"], col["wind_u"], col["wind_v"], surface_pressure_pa=col["surface_pressure_pa"],
            beta_steps=kw.pop("beta_steps", 60), **kw)
    f = lambda k: float(out[k][0, 0])
    return f("pft_gw"), f("z_fc_m"), f("delta_theta_fc_k"), f("theta_ml_k")

def elevated(col, ps_new):
    """Same column physics shifted to a surface pressure ps_new (elevated terrain)."""
    c = dict(col); s = ps_new / float(col["surface_pressure_pa"][0, 0])
    c["pressure_pa"] = col["pressure_pa"] * s
    th = col["temperature_k"] * (A._P0_PA / col["pressure_pa"]) ** A._KAPPA_DRY      # keep theta profile
    c["temperature_k"] = th * (c["pressure_pa"] / A._P0_PA) ** A._KAPPA_DRY
    es = 611.2 * np.exp(17.67 * (c["temperature_k"] - 273.15) / (c["temperature_k"] - 29.65))
    es0 = 611.2 * np.exp(17.67 * (col["temperature_k"] - 273.15) / (col["temperature_k"] - 29.65))
    rh = col["spec_humidity"] * col["pressure_pa"] / 0.622 / es0                 # approx RH
    c["spec_humidity"] = 0.622 * rh * es / c["pressure_pa"]
    c["surface_pressure_pa"] = np.full((1, 1), ps_new)
    return c

cases = {"hot-dry (rh25, abl1800)": _synthetic_column(),
         "moist (rh45, abl1200)": _synthetic_column(rh=45.0, abl=1200.0, gamma_free=4.5e-3)}
rows = []
for name, c0 in cases.items():
    for ps in (100000.0, 85000.0):
        c = c0 if ps == 100000.0 else elevated(c0, ps)
        tag = f"{name}, ps={ps/100:.0f} hPa"
        for g, kw in [("baseline (as coded)", {}),
                      ("p0=1013.25 hPa", dict(p0_pa=101325.0)),
                      ("p0=local ps", dict(p0_pa=ps)),
                      ("kappa forward=inverse=Rd/cp", dict(kappa=A._KAPPA_DRY)),
                      ("moisture per K of T, not theta", dict(dtheta_mode="T")),
                      ("p0=local ps + per-K-of-T", dict(p0_pa=ps, dtheta_mode="T"))]:
            install(**kw)
            for bs in (60, 600):
                rows.append((tag, g, bs) + pft(c, beta_steps=bs))
        A.theta_kelvin, A._mixed_layer_state, A._saturation_point = orig["theta"], orig["ml"], orig["sp"]
        for h in (12.0, 18.0):  # physical-uncertainty comparator: fire heat-to-moisture ratio
            rows.append((tag, f"h2m={h:.0f} K/(g/kg) [comparator]", 600) + pft(c, beta_steps=600, heat_to_moisture=h))
base = {}
print(f"{'case':38s} {'gauge':36s} {'steps':>5s} {'PFT GW':>8s} {'dPFT%':>7s} {'z_fc m':>7s} {'dth K':>6s} {'thML':>6s}")
for r in rows:
    tag, g, bs, p, z, d, th = r
    if g.startswith("baseline"): base[(tag, bs)] = p
    b = base.get((tag, bs)); dp = 100 * (p / b - 1) if b and np.isfinite(p) else float("nan")
    print(f"{tag:38s} {g:36s} {bs:5d} {p:8.1f} {dp:7.2f} {z:7.0f} {d:6.2f} {th:6.1f}")
