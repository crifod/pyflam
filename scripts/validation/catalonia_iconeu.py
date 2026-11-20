"""ICON-EU pyroconvection over Catalonia -- POTENTIAL only (no fuel gate).

The pyflam hybrid *atmosphere* (ICON-EU native model levels) applied outside Italy. ICON-2I
2.2 km is an Italy-only domain, so over Catalonia there is no ICON-2I fuel gate, surface field
or render grid: the atmosphere, surface state and grid are all ICON-EU, and only the potential
(atmospheric upper-bound) map is produced. Sea is masked via the land fraction.

The ICON-EU GRIB are full-Europe files, so a run already fetched for another region (e.g.
Tuscany) is reused directly -- only the read-time subset changes.

Usage:  PYTHONPATH=src python scripts/validation/catalonia_iconeu.py [YYYY-MM-DD] [run]
Env:    PYROCONV_OUT (output dir), PYROCONV_EU_CACHE (GRIB cache root).
"""
import os
import sys
import warnings
from datetime import datetime, timezone

import numpy as np

warnings.simplefilter("ignore")
from pyflam.atmosphere import (
    fetch_icon_eu, ICON_EU_MODEL_LEVELS, PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL,
    PYROCONVECTION_TYPE_COLOR, PYROCONVECTION_TYPE_LABEL)
# pyflam_gui is importable when running with PYTHONPATH=src from the repo root.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyflam_gui.core.pyroconv import read_icon_eu, iconeu_diagnostics, classify_profile

HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
BBOX = (42.9, 0.0, 40.4, 3.4)                 # Catalonia (N, W, S, E)

RUNDATE = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%d")
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d")
STAMP = f"{RUNDT:%Y%m%d}{RUN:02d}"
REPO = os.path.join(os.path.dirname(__file__), "..", "..")
OUT = os.environ.get("PYROCONV_OUT") or os.path.join(REPO, "docs", "daily")
EU_CACHE = os.environ.get("PYROCONV_EU_CACHE") or f"/tmp/pyflam_iconeu/{STAMP}"

# Overnight/shallow-ABL suppression. The default measured gradient (fit_in_ml) needs >=3
# model levels inside the mixed layer, which the ~200 m nocturnal ABL does not provide, so
# those cells fall to surface plume; the 600 m ABL floor reinforces it. PYROCONV_RELAX_NIGHT
# relaxes both -- it fills the gradient from the mixing-depth proxy where fit_in_ml is
# unmeasurable, and lowers the ABL floor -- so a stable nocturnal column is shown as a
# convection plume (a fire penetrating a stable layer) instead of a surface plume, matching
# the convention of the Catalan fire-service product.
RELAX = os.environ.get("PYROCONV_RELAX_NIGHT", "0") not in ("0", "", "false", "no")
ABL_MIN = float(os.environ.get("PYROCONV_ABL_MIN", "150" if RELAX else "600"))
# Moisture gate on classes 3-4 (resilient/deep). Tracks the shipped default
# (DEFAULT_PYROCONV_THRESHOLDS.rh_top_moist, now 60% -- lowered from 80% because 80 gated
# pyroCu/pyroCb out on dry Mediterranean afternoons where they occur). PYROCONV_RH_TOP_MOIST
# overrides it for sensitivity tests.
from dataclasses import replace as _replace
from pyflam.atmosphere import DEFAULT_PYROCONV_THRESHOLDS
_SHIPPED_RH = DEFAULT_PYROCONV_THRESHOLDS.rh_top_moist
RH_TOP_MOIST = float(os.environ.get("PYROCONV_RH_TOP_MOIST", str(_SHIPPED_RH)))
THRESH = _replace(DEFAULT_PYROCONV_THRESHOLDS, rh_top_moist=RH_TOP_MOIST)
_relax_tag = "_relaxed" if RELAX else ""
_moist_tag = "" if RH_TOP_MOIST == _SHIPPED_RH else f"_rh{int(RH_TOP_MOIST)}"
SUFFIX = _relax_tag + _moist_tag


def diagnostics_for_hour(step):
    files = fetch_icon_eu(RUNDT, run=RUN, step=step, cache_dir=os.path.join(EU_CACHE, f"{step:03d}"))
    d = read_icon_eu(files, BBOX, ICON_EU_MODEL_LEVELS)
    diag = iconeu_diagnostics(d, ml_method="fit_in_ml")
    if RELAX:
        # Fill the gradient from the always-computable mixing-depth proxy where the
        # measured fit is unavailable (too few levels in a shallow ML), and rebuild valid.
        proxy = iconeu_diagnostics(d, ml_method="surface_to_parcel")
        mlg = np.where(np.isfinite(diag["ml_grad"]), diag["ml_grad"], proxy["ml_grad"])
        diag["ml_grad"] = mlg
        diag["valid"] = (np.isfinite(diag["abl"]) & np.isfinite(diag["lcl_ratio"])
                         & np.isfinite(mlg))
    cls, used = classify_profile(diag, ladder="adaptive", abl_min_m=ABL_MIN,
                                 thresholds=THRESH)               # no fli
    sea = d["frland"] < 0.5
    return d["lat"], d["lon"], np.where(sea, -1, cls), diag, sea, sorted(used)


def render(lat, lon, stack, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    flip = lat[0] > lat[-1]
    ext = [lon.min(), lon.max(), lat.min(), lat.max()]
    colors = ["#cfe4ef"] + [PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-1.5, 5.5, 1), cmap.N)
    fig, ax = plt.subplots(1, len(HOURS), figsize=(2.2 * len(HOURS), 3.0),
                           constrained_layout=True, squeeze=False)
    ax = ax[0]
    for hi, hour in enumerate(HOURS):
        a = stack[hi][::-1] if flip else stack[hi]
        ax[hi].imshow(a, origin="lower", extent=ext, cmap=cmap, norm=norm,
                      aspect="auto", interpolation="nearest")
        ax[hi].set_title(f"{RUNDATE} {hour:02d}Z", fontsize=8)
        ax[hi].set_xticks([]); ax[hi].set_yticks([])
    notes = ("  [overnight relaxed]" if RELAX else "") + \
            (f"  [RH-top gate {int(RH_TOP_MOIST)}%]" if RH_TOP_MOIST != 80 else "")
    fig.suptitle("Pyroconvection type -- POTENTIAL (atmosphere only, no fuel gate)\n"
                 f"ICON-EU model levels -- Catalonia -- VALID {RUNDATE} "
                 f"(run {RUNDATE} {RUN:02d}Z){notes}", fontsize=11)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8.5, frameon=False,
               title="Pyroconvection class (0 = lowest activity -> 4 = highest)",
               bbox_to_anchor=(0.5, -0.12))
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"[catalonia] ICON-EU potential, VALID {RUNDATE}, bbox {BBOX}")
    stack, lat, lon = [], None, None
    for step in HOURS:                          # valid = run day, so step == hour
        lat, lon, cls, diag, sea, used = diagnostics_for_hour(step)
        stack.append(cls)
        land = ~sea
        n = max(int(land.sum()), 1)
        pct = [round(100 * ((cls == c) & land).sum() / n, 1) for c in range(5)]
        print(f"  {step:02d}Z  ~{diag['n_levels']} in ML  ladder={','.join(used) or '-'}  "
              f"classes%={pct}")
    png = os.path.join(OUT, f"pyroconv_catalonia_iconeu_potential{SUFFIX}_{RUNDATE}.png")
    render(lat, lon, np.stack(stack), png)
    print("wrote", png)


if __name__ == "__main__":
    main()
