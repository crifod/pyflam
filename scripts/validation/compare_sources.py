"""Cell-by-cell comparison of the icon2i and hybrid pyroconvection products for one day.

Reads the per-hour class GeoTIFFs (potential + gated) and the diagnostic rasters both
products write, restricts to land, and emits the tables the comparison report needs:
class distributions per hour, cell-level agreement, and diagnostic medians. Pure read +
arithmetic on the committed rasters, so the report is reproducible.
"""
import glob, os, sys, json
import numpy as np
import rasterio

DATE = sys.argv[1] if len(sys.argv) > 1 else "2026-07-15"
DAILY = os.environ.get("PYROCONV_DAILY",
    os.path.join(os.path.dirname(__file__), "..", "..", "docs", "daily"))
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
CLASSES = ["surface_plume", "convection_plume", "overshooting_pyrocu",
           "resilient_pyrocu", "deep_pyrocu_pyrocb"]
SHORT = ["surface", "convect", "overshoot", "resilient", "deep"]


def rd(path):
    with rasterio.open(path) as d:
        return d.read(1).astype(float)


def load(src, kind):
    """Class stack (nhours, ny, nx) for a source ('icon2i'/'hybrid') and kind ('potential'/'gated')."""
    rdir = os.path.join(DAILY, f"rasters_{src}_{DATE}")
    return np.stack([rd(os.path.join(rdir, f"pyroconv_{kind}_{h:02d}Z.tif")) for h in HOURS])


def load_diag(src, name):
    rdir = os.path.join(DAILY, f"rasters_{src}_{DATE}")
    return np.stack([rd(os.path.join(rdir, f"diag_{name}_{h:02d}Z.tif")) for h in HOURS])


def main():
    out = {"date": DATE, "hours": HOURS, "classes": CLASSES}

    # land mask: cells classifiable (>=0) in either product at midday, plus finite ABL
    a_pot = load("icon2i", "potential")
    h_pot = load("hybrid", "potential")
    a_gate = load("icon2i", "gated")
    h_gate = load("hybrid", "gated")
    # a cell is "land/valid" if it is ever classified >0 in the potential map of either
    land = ((a_pot > 0).any(0) | (h_pot > 0).any(0))
    n = int(land.sum())
    out["n_land_cells"] = n

    def dist(stack):
        return [[round(100 * ((stack[hi] == c) & land).sum() / max(n, 1), 1)
                 for c in range(5)] for hi in range(len(HOURS))]

    out["potential_icon2i"] = dist(a_pot)
    out["potential_hybrid"] = dist(h_pot)
    out["gated_icon2i"] = dist(a_gate)
    out["gated_hybrid"] = dist(h_gate)

    # cell-by-cell agreement (potential) per hour, and confusion at the peak hour
    agree = []
    for hi in range(len(HOURS)):
        m = land
        agree.append(round(100 * ((a_pot[hi] == h_pot[hi]) & m).sum() / max(n, 1), 1))
    out["agreement_potential"] = agree

    peak = HOURS.index(12)
    conf = np.zeros((5, 5), int)     # rows icon2i, cols hybrid
    ai, hj = a_pot[peak][land].astype(int), h_pot[peak][land].astype(int)
    for i in range(5):
        for j in range(5):
            conf[i, j] = int(((ai == i) & (hj == j)).sum())
    out["confusion_12Z_potential"] = conf.tolist()

    # diagnostic medians over land, per hour
    diags = {}
    for name in ("abl", "parcel_ml", "lcl", "lcl_ratio", "ml_grad", "gamma", "rh_top"):
        try:
            ai = load_diag("icon2i", name); hi = load_diag("hybrid", name)
        except Exception:
            continue
        diags[name] = {
            "icon2i": [round(float(np.nanmedian(ai[k][land])), 4) for k in range(len(HOURS))],
            "hybrid": [round(float(np.nanmedian(hi[k][land])), 4) for k in range(len(HOURS))],
        }
    out["diag_medians"] = diags

    # deep + resilient (classes 3-4) fraction, gated, per hour -- the operational signal
    def hi_frac(stack):
        return [round(100 * (((stack[k] == 3) | (stack[k] == 4)) & land).sum() / max(n, 1), 2)
                for k in range(len(HOURS))]
    out["gated_hi_icon2i"] = hi_frac(a_gate)
    out["gated_hi_hybrid"] = hi_frac(h_gate)

    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
