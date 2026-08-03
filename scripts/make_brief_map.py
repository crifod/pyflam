# SPDX-License-Identifier: AGPL-3.0-or-later
"""Publication-quality isochrone map for the operational brief."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LinearSegmentedColormap
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.environ.setdefault("PYFLAM_NOCFD", "1")
# Respect the caller's AOI config (must match the checkpoint); fall back to the
# production DS=1 / asymmetric-margin framing when unset.
for _k, _v in dict(PYFLAM_DS="1", PYFLAM_MARGIN_W="16", PYFLAM_MARGIN_E="6",
                   PYFLAM_MARGIN_N="9", PYFLAM_MARGIN_S="11").items():
    os.environ.setdefault(_k, _v)
import run_fire_11_09_2pm as R

OUT = R.OUT
ls, ign, xy = R.clip_landscape()
d = np.load(os.path.join(OUT, "march_checkpoint.npz"))
arr = d["arrival"]; cs = d["cellsize"]
times = list(d["times"])

ember = LinearSegmentedColormap.from_list(
    "ember", ["#f7d774", "#f2a63b", "#e2571f", "#a11414", "#5c0a0a"])

fig, ax = plt.subplots(figsize=(9.2, 8.4))
elev = np.asarray(ls.elevation, float)
sh = LightSource(azdeg=315, altdeg=45)
ax.imshow(sh.hillshade(elev, vert_exag=2.2, dx=cs[0], dy=cs[1]),
          cmap="gray", alpha=0.85, zorder=0)
# non-burnable / burnable subtle fuel tint
fuel = np.asarray(ls.fuel_model)
nonburn = np.isin(fuel, [91, 92, 93, 98, 99]) | (fuel <= 0)
tint = np.zeros(fuel.shape + (4,)); tint[..., 1] = 0.45; tint[..., 3] = np.where(nonburn, 0.0, 0.10)
ax.imshow(tint, zorder=1)

am = np.where(np.isfinite(arr), arr, np.nan)
cf = ax.contourf(am, levels=times, cmap=ember, alpha=0.78, zorder=2, extend="neither")
cs_lines = ax.contour(am, levels=times, colors="k", linewidths=0.5, alpha=0.5, zorder=3)
# ignition
ax.plot(ign[1], ign[0], marker="*", color="#111", ms=22, mec="white", mew=1.4, zorder=5)
ax.annotate("ignition\n14:00", (ign[1], ign[0]), textcoords="offset points",
            xytext=(10, 8), fontsize=9, weight="bold", color="#111",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7), zorder=6)

# scale bar (5 km)
px_per_km = 1000.0 / cs[0]
x0 = fuel.shape[1] * 0.06; y0 = fuel.shape[0] * 0.94
ax.plot([x0, x0 + 5 * px_per_km], [y0, y0], color="#111", lw=4, zorder=6, solid_capstyle="butt")
ax.text(x0 + 2.5 * px_per_km, y0 - 6, "5 km", ha="center", va="bottom",
        fontsize=10, weight="bold", color="#111", zorder=6)
# north arrow
nx, ny = fuel.shape[1] * 0.95, fuel.shape[0] * 0.14
ax.annotate("N", (nx, ny - 26), ha="center", fontsize=12, weight="bold", color="#111", zorder=6)
ax.annotate("", (nx, ny - 22), (nx, ny + 18),
            arrowprops=dict(arrowstyle="-|>", color="#111", lw=2.4), zorder=6)

cb = fig.colorbar(cf, ax=ax, fraction=0.035, pad=0.02, ticks=times[::2])
cb.set_label("time of arrival (minutes after 14:00 CEST)", fontsize=10)
cb.ax.tick_params(labelsize=8)

ax.set_title("Calvana wildfire — 30-minute fire-front isochrones\n"
             "ICON-2i · Cruz-2005 crown fire · Tuscany fuel + GEDI canopy · 100 m",
             fontsize=12.5, weight="bold", pad=10)
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_edgecolor("#bbb")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "brief_map.png"), dpi=150, bbox_inches="tight")
print("wrote", os.path.join(OUT, "brief_map.png"))
