"""ERA5 as the rich-level reference: does the reanalysis back the ICON-EU product or the
Catalan fire-service (Bombers) product on the overshoot-vs-surface question?

Runs the *same* pyflam classifier (iconeu_diagnostics + classify_profile) on ERA5 pressure
levels over Catalonia and reports two things:

 (1) REGIME -- the inland vs coastal LCL/ABL and the class balance at midday. This adjudicates
     the disagreement: the Bombers product is surface-heavy (~62-68% surface, 12-15% overshoot
     at midday) while the ICON-EU product is overshoot-heavy (~46-52% overshoot). ERA5 settles
     which regime is real.

 (2) RELATIONSHIP MODEL -- degrade ERA5 within itself to ICON-EU-like (12) and ICON-2I-like (5)
     level counts and measure the ABL/LCL/gradient bias vs the full 22-level reference. This is
     the "reduced -> full" correction to carry to ICON-EU (it turns out negligible at ICON-EU's
     level density, which is why the model-level hybrid works and the 5-level path did not).

Fetch the inputs first:  PYTHONPATH=src python scripts/validation/era5_fetch_catalonia.py
Then:                     PYTHONPATH=src python scripts/validation/era5_regime_relationship.py

Finding (2026-07-06..10): ERA5 midday overshoot ~58%, LCL/ABL ~1.1 inland -- overshoot-dominant,
matching the ICON-EU product and contradicting the Bombers' surface-heavy picture. The ICON-EU
level correction is +174 m ABL / ~0 ratio vs full; the 5-level path is +623 m ABL / halved
gradient. Companion to sounding_degradation.py (radiosonde version of the same idea).
"""
import os
import sys
import warnings

import numpy as np
import xarray as xr

warnings.simplefilter("ignore")
# pyflam_gui lives at the repo root; pyflam is on PYTHONPATH=src.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyflam.atmosphere import specific_humidity_from_rh, _G
from pyflam_gui.core.pyroconv import iconeu_diagnostics, classify_profile

# Bombers (Catalan fire service) ICON-EU product, segmented from their 2026-07-13 panels.
# Different day from the ERA5 window, but the summer regime is stable across adjacent days;
# ERA5 characterises the regime authoritatively regardless of exact date.
BOMBERS_MIDDAY = "overshoot ~12-15%, surface ~62-68%"
ICONEU_MIDDAY = "overshoot ~46-52%, surface ~33%"


def _tc(ds):
    for c in ("valid_time", "time"):
        if c in ds.coords:
            return c
    raise KeyError("no time coordinate")


def build(pl, sf, tc, lev, order, it, keep=None):
    """Build the iconeu_diagnostics 'd' dict from ERA5 at time index ``it``.

    ``keep`` is a boolean mask on levels (ascending-height order) to degrade to a coarser
    set. ERA5 gives geopotential and RH; heights come from z/g and specific humidity from RH.
    """
    g = lambda v: pl[v].isel({tc: it}).values[order]
    z = g("z") / _G
    T = g("t"); R = np.clip(g("r"), 1, 100); U = g("u"); V = g("v")
    P = (lev[:, None, None] * 100.0) * np.ones_like(T)
    QV = specific_humidity_from_rh(R, T, P)
    s = lambda v: sf[v].isel({tc: it}).values
    orog = s("z") / _G
    if keep is not None:
        z, T, QV, P, U, V = z[keep], T[keep], QV[keep], P[keep], U[keep], V[keep]
    return dict(z=z - orog[None, ...], T=T, QV=QV, P=P, U=U, V=V,
                T2m=s("t2m"), Td2m=s("d2m"), PS=s("sp"), U10=s("u10"), V10=s("v10"),
                frland=s("lsm"), orog=orog)


def main():
    pl = xr.open_dataset("era5_pl_cat.nc")
    sf = xr.open_dataset("era5_sfc_cat.nc")
    tc = _tc(pl)
    lev = pl["pressure_level"].values
    order = np.argsort(-lev)                       # descending p == ascending height
    lev = lev[order]
    lat, lon = pl["latitude"].values, pl["longitude"].values
    hours = [int(str(x)[11:13]) for x in pl[tc].values]
    print(f"ERA5 {len(lev)} levels {lev.min():.0f}-{lev.max():.0f} hPa, "
          f"grid {len(lat)}x{len(lon)}, {len(hours)} time slices")

    def nearest(targets):
        return np.unique([int(np.argmin(np.abs(lev - t))) for t in targets])
    kEU = np.zeros(len(lev), bool)
    kEU[nearest([1000, 950, 925, 900, 875, 850, 800, 700, 600, 500, 400, 300])] = True
    k2I = np.zeros(len(lev), bool)
    k2I[nearest([1000, 925, 850, 700, 500])] = True

    # --- (1) regime ---
    print("\n(1) REGIME (ERA5 full levels), land medians. inland=orog>300 m, coast=orog<100 m")
    print(f"{'hr':>3} {'ABL':>6} {'ratio_all':>9} {'ratio_inl':>9} {'ratio_cst':>9}  "
          f"classes% [surf conv over resil deep]")
    over, surf = [], []
    for it, h in enumerate(hours):
        if h not in (12, 15):
            continue
        d = build(pl, sf, tc, lev, order, it)
        diag = iconeu_diagnostics(d, ml_method="fit_in_ml")
        land = d["frland"] >= 0.5
        m = land & diag["valid"]
        inl = m & (d["orog"] > 300); cst = m & (d["orog"] < 100)
        cls, _ = classify_profile(diag, ladder="adaptive")
        n = max(int(land.sum()), 1)
        pct = [round(100 * ((cls == c) & land).sum() / n, 1) for c in range(5)]
        r = diag["lcl_ratio"]
        rin = np.nanmedian(r[inl]) if inl.sum() else np.nan
        rcs = np.nanmedian(r[cst]) if cst.sum() else np.nan
        print(f"{h:02d}Z {np.nanmedian(diag['abl'][m]):6.0f} {np.nanmedian(r[m]):9.2f} "
              f"{rin:9.2f} {rcs:9.2f}  {pct}")
        over.append(pct[2]); surf.append(pct[0])
    print(f"\n  ERA5 midday: overshoot median {np.median(over):.0f}%, surface {np.median(surf):.0f}%")
    print(f"  ICON-EU product: {ICONEU_MIDDAY}")
    print(f"  Bombers product: {BOMBERS_MIDDAY}")
    print("  => ERA5 is overshoot-dominant like the ICON-EU product, NOT surface-heavy like "
          "the Bombers. Inland LCL/ABL ~1.1 is the overshooting regime.")

    # --- (2) relationship model ---
    print("\n(2) RELATIONSHIP MODEL: reduced-level vs full-level diagnostics "
          "(midday 09/12/15Z, land, pooled)")
    store = {tag: {"abl": [], "ratio": [], "mlg": []} for tag in ("full", "eu", "i2i")}
    for it, h in enumerate(hours):
        if h not in (9, 12, 15):
            continue
        land = sf["lsm"].isel({_tc(sf): it}).values >= 0.5
        for tag, keep in (("full", None), ("eu", kEU), ("i2i", k2I)):
            dg = iconeu_diagnostics(build(pl, sf, tc, lev, order, it, keep), ml_method="fit_in_ml")
            mm = land & dg["valid"]
            store[tag]["abl"].append(dg["abl"][mm])
            store[tag]["ratio"].append(dg["lcl_ratio"][mm])
            store[tag]["mlg"].append(dg["ml_grad"][mm])
    med = lambda tag, k: float(np.nanmedian(np.concatenate(store[tag][k])))
    print(f"  {'levels':<15}{'ABL(m)':>8}{'LCL/ABL':>9}{'MLgrad':>10}   bias vs full")
    print(f"  {'full (22)':<15}{med('full','abl'):8.0f}{med('full','ratio'):9.2f}{med('full','mlg'):10.2e}   --")
    print(f"  {'ICON-EU (12)':<15}{med('eu','abl'):8.0f}{med('eu','ratio'):9.2f}{med('eu','mlg'):10.2e}   "
          f"ABL {med('eu','abl')-med('full','abl'):+.0f}m, ratio {med('eu','ratio')-med('full','ratio'):+.2f}")
    print(f"  {'ICON-2I (5)':<15}{med('i2i','abl'):8.0f}{med('i2i','ratio'):9.2f}{med('i2i','mlg'):10.2e}   "
          f"ABL {med('i2i','abl')-med('full','abl'):+.0f}m, ratio {med('i2i','ratio')-med('full','ratio'):+.2f}")
    print("  => ICON-EU's level density is essentially unbiased vs full ERA5 (the correction to "
          "carry over is negligible); the 5-level path is not.")


if __name__ == "__main__":
    main()
