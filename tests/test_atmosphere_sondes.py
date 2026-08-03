"""Regression tests: the atmosphere pipeline on real GRAF field soundings.

Runs the full pyflam atmosphere stack on vendored ambient radiosondes (Zenodo
15264835, CC-BY; see tests/data/graf_sondes/ATTRIBUTION.md) and asserts physical
plausibility. Unlike the synthetic-column tests elsewhere, these exercise the
module against real, high-resolution, dry-fire-weather profiles it has never seen.
"""
from __future__ import annotations

import glob
import os

import numpy as np
import pytest

from pyflam.atmosphere import (
    theta_kelvin, dewpoint_from_rh, lcl_height_bolton_m, specific_humidity_from_rh,
    virtual_potential_temperature, bulk_richardson_abl_grid, parcel_mixing_depth_grid,
    shear_height_window, fire_induced_abl_grid, mixed_layer_fire_flux,
    pyroconvection_type, PYROCONVECTION_TYPES,
)

_DATA = os.path.join(os.path.dirname(__file__), "data", "graf_sondes")
_SONDES = sorted(glob.glob(os.path.join(_DATA, "*_Environment*.raw_flight_history.csv")))
_G1 = lambda a: np.asarray(a, float)[:, None, None]


def _load(path):
    """Real balloon flight history -> monotone-in-height ambient profile arrays."""
    import csv
    z, T, RH, P, spd, hdg = [], [], [], [], [], []
    with open(path) as f:
        rd = csv.reader(f)
        header = [c.strip() for c in next(rd)]
        ix = {name: header.index(col) for name, col in {
            "z": "Altitude (m AGL)", "T": "Temperature (C)", "RH": "Relative humidity (%)",
            "P": "Pressure (Pascal)", "spd": "Speed (m/s)", "hdg": "Heading (degrees)"}.items()}

        def num(row, k):
            try:
                return float(row[ix[k]])
            except (ValueError, IndexError):
                return np.nan
        for row in rd:
            z.append(num(row, "z")); T.append(num(row, "T")); RH.append(num(row, "RH"))
            P.append(num(row, "P")); spd.append(num(row, "spd")); hdg.append(num(row, "hdg"))
    z, T, RH, P, spd, hdg = (np.array(a, float) for a in (z, T, RH, P, spd, hdg))
    ok = np.isfinite(z) & np.isfinite(T) & np.isfinite(P) & (P > 1e4) & (z >= 0)
    z, T, RH, P, spd, hdg = (a[ok] for a in (z, T, RH, P, spd, hdg))
    order = np.argsort(z)
    z, T, RH, P, spd, hdg = (a[order] for a in (z, T, RH, P, spd, hdg))
    RH = np.clip(RH, 1, 100)
    u = np.where(np.isfinite(spd), spd * np.sin(np.deg2rad(hdg)), 0.0)
    v = np.where(np.isfinite(spd), spd * np.cos(np.deg2rad(hdg)), 0.0)
    th = theta_kelvin(T, P / 100.0)
    q = specific_humidity_from_rh(RH, T + 273.15, P)
    thv = virtual_potential_temperature(th, q)
    return dict(z=z, T=T, RH=RH, P=P, u=u, v=v, th=th, thv=thv, q=q)


def _diagnose(p):
    z, th, thv, u, v = p["z"], p["th"], p["thv"], p["u"], p["v"]
    sfc = z <= z.min() + 120
    Ts, RHs, Ps = p["T"][sfc].mean(), p["RH"][sfc].mean(), p["P"][sfc].mean()
    lcl = float(lcl_height_bolton_m(Ts + 273.15, dewpoint_from_rh(Ts, RHs) + 273.15, Ps))
    abl = float(bulk_richardson_abl_grid(
        _G1(z), _G1(thv), _G1(u), _G1(v),
        theta_v_surface=np.array([[thv[sfc].mean()]]),
        wind_u_surface=np.array([[u[sfc].mean()]]),
        wind_v_surface=np.array([[v[sfc].mean()]]))[0, 0])
    parcel = float(parcel_mixing_depth_grid(_G1(z), _G1(th), np.array([[th[sfc].mean()]]))[0, 0])
    cap_m = (z >= abl + 200) & (z <= abl + 1200)
    cap = float(np.polyfit(z[cap_m], th[cap_m], 1)[0]) if cap_m.sum() >= 2 else np.nan
    q_ref = mixed_layer_fire_flux(1.0e7, abl)              # 10 MW/m reference fire
    fabl = float(fire_induced_abl_grid(
        _G1(z), _G1(th), theta_mean_below=np.array([[th[z <= abl].mean()]]),
        heat_flux=np.array([[q_ref]]), blh=np.array([[abl]]))[0, 0])
    return dict(lcl=lcl, abl=abl, parcel=parcel, cap=cap, fabl=fabl)


@pytest.fixture(scope="module")
def sonde_paths():
    assert _SONDES, "vendored GRAF sonde fixtures missing"
    return _SONDES


def test_fixtures_present(sonde_paths):
    assert len(sonde_paths) >= 6                            # 6 fires, 8 sondes vendored


@pytest.mark.parametrize("path", _SONDES, ids=[os.path.basename(p).split("_")[0] for p in _SONDES])
def test_pipeline_physics_on_real_sonde(path):
    """Every core atmosphere diagnostic yields physically plausible values."""
    p = _load(path)
    assert p["z"].size >= 8
    fin = np.isfinite(p["thv"])
    assert np.all(p["thv"][fin] >= p["th"][fin] - 1e-6)    # theta_v >= theta
    d = _diagnose(p)
    assert 150.0 <= d["abl"] <= 4500.0                      # ABL in plausible band
    assert 150.0 <= d["parcel"] <= 4500.0
    assert 0.0 < d["lcl"] < 8000.0                          # LCL positive, sane
    assert np.isnan(d["cap"]) or d["cap"] >= -1.0e-4        # free troposphere ~stable
    assert d["fabl"] >= d["abl"] - 1e-6                     # fireABL cannot sit below ABL


@pytest.mark.parametrize("path", _SONDES, ids=[os.path.basename(p).split("_")[0] for p in _SONDES])
def test_classifier_returns_valid_class(path):
    """The Castellnou ladder returns a defined class on every real column."""
    p = _load(path)
    d = _diagnose(p)
    z, th = p["z"], p["th"]
    ml_m = (z >= z.min()) & (z <= d["parcel"])
    ml_grad = float(np.polyfit(z[ml_m], th[ml_m], 1)[0]) if ml_m.sum() >= 2 else 0.0
    rh_top = float(p["RH"][(z >= d["abl"] - 150) & (z <= d["abl"] + 150)].mean())
    cls = pyroconvection_type(
        lcl_abl_ratio=d["lcl"] / d["abl"], ml_theta_gradient=ml_grad,
        gamma_theta=d["cap"] if np.isfinite(d["cap"]) else 5e-3,
        rh_top_abl=rh_top if np.isfinite(rh_top) else 30.0,
        shear_distance=1.0, fireline_intensity_kw=1.0e4, ladder="adaptive")
    cls = cls[0] if isinstance(cls, tuple) else cls
    assert cls in PYROCONVECTION_TYPES


def test_granjaescarp_reference_values():
    """Golden regression: GranjaEscarp's ambient diagnostics stay in a tight band."""
    path = [p for p in _SONDES if "GranjaEscarp" in p][0]
    d = _diagnose(_load(path))
    assert 400 < d["abl"] < 900                            # ~640 m
    assert 3000 < d["lcl"] < 4800                          # ~3.9 km (very dry, RH ~12%)
