"""Montale 2017 — maps, observed-vs-simulated comparison, and combined report."""
import os, sys, numpy as np, warnings
warnings.simplefilter("ignore")
import rasterio, geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

BASE = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "case_montale_2017")
PERIM = "/Users/cristianofoderi/Dropbox/BKgen19/Desktop/incendio_Montale/perim_lamma.shp"
IGN_3035 = (4404680.0, 2317385.0)
CELL_HA = 10 * 10 / 1e4


def load(name, mod="crown_fire"):
    with rasterio.open(os.path.join(BASE, mod, name)) as d:
        return d.read(1), d.bounds, d.transform


def main():
    arr, bnds, tr = load("crown_arrival_time_min.tif")
    fli, _, _ = load("crown_fireline_intensity_kW_m.tif")
    ftype, _, _ = load("fire_type.tif")
    ext = [bnds.left, bnds.right, bnds.bottom, bnds.top]
    perim = gpd.read_file(PERIM).to_crs(3035)
    obs_ha = perim.area.sum() / 1e4

    # area vs time; time to reach the observed 286 ha
    finite = np.isfinite(arr)
    times = np.sort(arr[finite])
    area_series = (np.arange(1, times.size + 1)) * CELL_HA
    t_obs = float(times[np.searchsorted(area_series, obs_ha)]) if obs_ha < area_series[-1] else float("inf")
    isochrones = [30, 60, 120, 240]

    # --- Map A: arrival isochrones + observed perimeter ---
    fig, ax = plt.subplots(figsize=(8, 8))
    am = np.where(finite, arr, np.nan)
    im = ax.imshow(am, extent=ext, origin="upper", cmap="turbo", vmin=0, vmax=300)
    cs = ax.contour(np.clip(am, 0, 400), levels=isochrones, extent=ext, origin="upper",
                    colors="k", linewidths=0.8)
    ax.clabel(cs, fmt="%d min", fontsize=7)
    perim.boundary.plot(ax=ax, color="white", linewidth=2.2, label="observed perimeter (286 ha)")
    ax.plot(*IGN_3035, "r*", ms=18, label="ignition")
    ax.set_title("Montale 2017 — simulated arrival time (min) + observed perimeter", fontsize=11)
    ax.legend(loc="upper right"); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, shrink=0.7, label="arrival time (min)")
    fig.savefig(os.path.join(BASE, "reports", "map_arrival_vs_observed.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Map B: crown fire type ---
    fig, ax = plt.subplots(figsize=(7, 7))
    cmap = ListedColormap(["#f7f7f7", "#fdae61", "#d7301f"]); norm = BoundaryNorm([-.5, .5, 1.5, 2.5], 3)
    ax.imshow(np.where(finite, ftype, np.nan), extent=ext, origin="upper", cmap=cmap, norm=norm)
    perim.boundary.plot(ax=ax, color="black", linewidth=1.5)
    ax.plot(*IGN_3035, "b*", ms=14)
    ax.set_title("Crown fire type (grey surface / orange passive / red active)", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(os.path.join(BASE, "reports", "map_crown_fire_type.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Map C: fireline intensity ---
    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(np.where(finite & (fli > 0), fli, np.nan), extent=ext, origin="upper",
                   cmap="inferno", vmin=0, vmax=30000)
    perim.boundary.plot(ax=ax, color="cyan", linewidth=1.5)
    ax.plot(*IGN_3035, "w*", ms=14)
    ax.set_title("Fireline intensity (kW/m, capped)", fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, shrink=0.7, label="kW/m")
    fig.savefig(os.path.join(BASE, "reports", "map_fireline_intensity.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # area-vs-time report
    area_at = {t: float((arr[finite] <= t).sum() * CELL_HA) for t in isochrones}
    open(os.path.join(BASE, "reports", "spread_vs_observed.md"), "w").write(f"""# Montale 2017 — simulated spread vs observed

Observed final perimeter (LAMMA): **{obs_ha:.0f} ha** (reported ~300 ha forest).

## Free-burning simulation (no suppression)
| time | simulated burned area |
|---|---|
| 30 min | {area_at[30]:.0f} ha |
| 60 min | {area_at[60]:.0f} ha |
| 120 min | {area_at[120]:.0f} ha |
| 240 min | {area_at[240]:.0f} ha |

- Simulated time to reach the observed {obs_ha:.0f} ha: **~{t_obs:.0f} min** after ignition.
- The free-burning sim keeps spreading past 286 ha (no suppression modelled); the
  observed fire was held near 286 ha by ~230 fire crews over hours. The **shape and
  direction** (SW spread under the NE wind) and the fireline intensities are the
  physically meaningful comparison, not the free-burn final area.

Maps: map_arrival_vs_observed.png, map_crown_fire_type.png, map_fireline_intensity.png
""")
    print(f"observed {obs_ha:.0f} ha; sim reaches it at ~{t_obs:.0f} min; maps written")
    return obs_ha, t_obs, area_at


if __name__ == "__main__":
    main()
