"""3-day Catalonia pyroconvection forecast -- ICON-EU only (potential map).

Catalonia is outside the ICON-2I 2.2 km domain (Italy-only) and there is no Catalan
fuel .lcp, so this is the ICON-EU-native, potential-only counterpart to the Tuscany
hybrid product: profile diagnostics (incl. the dry-pyrocloud fireABL/decoupling) from
the ICON-EU model levels, classified as the atmospheric upper bound. Reuses the same
ICON-EU cache the Tuscany run downloads (read_icon_eu just subsets a different bbox).

Usage: PYTHONPATH=src python scripts/catalunya_iconeu_3day.py [run-day] [run]
"""
import os, sys
from datetime import datetime, timedelta, timezone
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(_REPO, "src"))      # pyflam
sys.path.insert(0, _REPO)                            # pyflam_gui (repo root)
from pyflam.atmosphere import (
    fetch_icon_eu, ICON_EU_MODEL_LEVELS,
    PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_COLOR,
    PYROCONVECTION_TYPE_LEVEL, PYROCONVECTION_TYPE_LABEL)
from pyflam_gui.core.pyroconv import read_icon_eu, iconeu_diagnostics, classify_profile

RUNDATE = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%d")
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d").replace(hour=RUN)
HOURS = [int(h) for h in os.environ.get("PYROCONV_HOURS", "9,12,15").split(",")]
DAYS = [(RUNDT + timedelta(days=k)).strftime("%Y-%m-%d") for k in (0, 1, 2)]
# Catalonia bounding box (n, w, s, e)
NORTH, WEST, SOUTH, EAST = 42.95, 0.0, 40.4, 3.35
EU_CACHE = os.environ.get("PYROCONV_EU_CACHE") or f"/tmp/pyflam_iconeu/{RUNDT:%Y%m%d}{RUN:02d}"
OUT = os.environ.get("PYROCONV_OUT") or os.path.join(
    os.path.dirname(__file__), "..", "docs", f"catalunya_{RUNDATE}_3day")
os.makedirs(OUT, exist_ok=True)
RASTERDIR = os.path.join(OUT, f"rasters_iconeu_{RUNDATE}")
os.makedirs(RASTERDIR, exist_ok=True)


def diag_for(valid, h):
    vdt = datetime.strptime(valid, "%Y-%m-%d").replace(hour=h)
    step = int((vdt - RUNDT).total_seconds() // 3600)
    files = fetch_icon_eu(RUNDT, run=RUN, step=step, cache_dir=os.path.join(EU_CACHE, f"{step:03d}"))
    eud = read_icon_eu(files, (NORTH, WEST, SOUTH, EAST), ICON_EU_MODEL_LEVELS)
    diag = iconeu_diagnostics(eud, ml_method="fit_in_ml")
    return eud, diag


def main():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import rasterio
    from rasterio.transform import from_origin

    cmap = ListedColormap([PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    nrow, ncol = len(DAYS), len(HOURS)
    fig, ax = plt.subplots(nrow, ncol, figsize=(2.4*ncol, 2.4*nrow),
                           constrained_layout=True, squeeze=False)
    lat = lon = ext = None
    for di, valid in enumerate(DAYS):
        for hi, h in enumerate(HOURS):
            eud, diag = diag_for(valid, h)
            lat, lon = eud["lat"], eud["lon"]
            cls, used = classify_profile(diag, ladder="adaptive")   # potential (no gate)
            flip = lat[0] > lat[-1]
            ext = [lon.min(), lon.max(), lat.min(), lat.max()]
            a = cls[::-1] if flip else cls
            ax[di][hi].imshow(a, origin="lower", extent=ext, cmap=cmap, norm=norm,
                              aspect="auto", interpolation="nearest")
            ax[di][hi].set_title(f"{valid} {h:02d}Z", fontsize=8)
            ax[di][hi].set_xticks([]); ax[di][hi].set_yticks([])
            # rasters: class + the dry-pyrocloud decoupling ratio
            dlon = float(abs(lon[1]-lon[0])); dlat = float(abs(lat[1]-lat[0]))
            tr = from_origin(lon.min()-dlon/2, lat.max()+dlat/2, dlon, dlat)
            arr = cls if flip else cls[::-1]
            for tag, data, dt in [("potential", arr, "int16"),
                                  ("decoupling", (diag["decoupling"] if flip else diag["decoupling"][::-1]), "float32")]:
                fn = os.path.join(RASTERDIR, f"cat_{tag}_{valid}_{h:02d}Z.tif")
                with rasterio.open(fn, "w", driver="GTiff", height=a.shape[0], width=a.shape[1],
                                   count=1, crs="EPSG:4326", transform=tr, dtype=dt,
                                   nodata=(-1 if dt == "int16" else None)) as d:
                    d.write(np.asarray(data, dt), 1)
            sys.stderr.write(f"  {valid} {h:02d}Z  ladder={','.join(used)}\n")
    fig.suptitle("Catalonia pyroconvection type -- POTENTIAL (atmosphere only, upper bound)\n"
                 f"ICON-EU 6.5 km model levels (DWD open data) -- run {RUNDATE} {RUN:02d}Z", fontsize=12)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, -0.04))
    png = os.path.join(OUT, f"pyroconv_catalunya_iconeu_potential_{RUNDATE}.png")
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"OK Catalonia 3-day potential -> {png}")
    return png


if __name__ == "__main__":
    main()
