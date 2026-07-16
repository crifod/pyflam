"""3-day Tuscany pyroconvection forecast (hybrid) from a single 00Z run, as one report.

Runs the tested daily pipeline (tests/pyroconv_daily.py) for three valid days -- the run day
and the next two forecast days -- then stitches the three days into a combined 3-row figure
(potential + gated) and a markdown/PDF report with per-day class distributions.

Uses the day+1/day+2 forecast steps available in the same ICON-EU + ICON-2I 00Z run (ICON-EU to
+120 h, ICON-2I to +72 h). Output goes to its own folder.

Usage:  PYTHONPATH=src python scripts/forecast_3day_tuscany.py [YYYY-MM-DD run-day] [run]
Env:    PYROCONV_OUT (folder, default docs/forecast_<rundate>_3day), plus the usual
        PYROCONV_* knobs honoured by tests/pyroconv_daily.py.
"""
import os
import subprocess
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from pyflam.atmosphere import (PYROCONVECTION_TYPES, PYROCONVECTION_TYPE_LEVEL,
                               PYROCONVECTION_TYPE_COLOR, PYROCONVECTION_TYPE_LABEL)

RUNDATE = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%d")
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
RUNDT = datetime.strptime(RUNDATE, "%Y-%m-%d")
DAYS = [(RUNDT + timedelta(days=k)).strftime("%Y-%m-%d") for k in (0, 1, 2)]
HOURS = [0, 3, 6, 9, 12, 15, 18, 21]
OUT = os.environ.get("PYROCONV_OUT") or os.path.join(REPO, "docs", f"forecast_{RUNDATE}_3day")


def run_daily(valid):
    """Invoke the daily pipeline for one valid day into OUT (hybrid default)."""
    env = dict(os.environ, PYROCONV_OUT=OUT, PYROCONV_VALID=valid)
    print(f"\n=== {valid} (forecast day +{(datetime.strptime(valid,'%Y-%m-%d')-RUNDT).days}) ===")
    subprocess.run([sys.executable, os.path.join(REPO, "tests", "pyroconv_daily.py"),
                    RUNDATE, str(RUN), valid], env=env, check=True)


def load_stack(valid, kind):
    rdir = os.path.join(OUT, f"rasters_hybrid_{valid}")
    arrs, transform = [], None
    for h in HOURS:
        with rasterio.open(os.path.join(rdir, f"pyroconv_{kind}_{h:02d}Z.tif")) as d:
            arrs.append(d.read(1)); transform = d.transform
    return np.stack(arrs), transform


def combined_figure(kind, path):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    cmap = ListedColormap([PYROCONVECTION_TYPE_COLOR[t] for t in PYROCONVECTION_TYPES])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    fig, ax = plt.subplots(len(DAYS), len(HOURS), figsize=(2.05 * len(HOURS), 2.7 * len(DAYS)),
                           constrained_layout=True, squeeze=False)
    for r, valid in enumerate(DAYS):
        stack, _ = load_stack(valid, kind)
        for c, hour in enumerate(HOURS):
            ax[r][c].imshow(stack[c], origin="upper", cmap=cmap, norm=norm,
                            aspect="auto", interpolation="nearest")
            ax[r][c].set_xticks([]); ax[r][c].set_yticks([])
            if r == 0:
                ax[r][c].set_title(f"{hour:02d}Z", fontsize=9)
            if c == 0:
                ax[r][c].set_ylabel(f"{valid}\n(+{r}d)", fontsize=9)
    title = ("FUEL-GATED (expected)" if kind == "gated"
             else "POTENTIAL (atmospheric upper bound)")
    fig.suptitle(f"Tuscany pyroconvection -- 3-day forecast -- {title}\n"
                 f"ICON-EU model levels + ICON-2I 2.2 km gate -- run {RUNDATE} {RUN:02d}Z",
                 fontsize=12)
    leg = [Patch(facecolor=PYROCONVECTION_TYPE_COLOR[t], edgecolor="0.4",
                 label=f"{PYROCONVECTION_TYPE_LEVEL[t]}  {PYROCONVECTION_TYPE_LABEL[t]}")
           for t in PYROCONVECTION_TYPES]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, -0.04))
    fig.savefig(path, dpi=135, bbox_inches="tight"); import matplotlib.pyplot as p; p.close(fig)


def dist_table(kind):
    """Per-day, per-hour class % over land (land = cells ever classed > 0 across the 3 days)."""
    stacks = {v: load_stack(v, "potential")[0] for v in DAYS}
    land = np.zeros(next(iter(stacks.values()))[0].shape, bool)
    for s in stacks.values():
        land |= (s > 0).any(0)
    n = max(int(land.sum()), 1)
    rows = []
    for valid in DAYS:
        st, _ = load_stack(valid, kind)
        for hi, h in enumerate(HOURS):
            pct = [round(100 * ((st[hi] == c) & land).sum() / n, 1) for c in range(5)]
            rows.append((valid, h, pct))
    return rows, n


def build_report(png_pot, png_gate):
    md = os.path.join(OUT, f"forecast_3day_{RUNDATE}.md")
    pdf = os.path.join(OUT, f"forecast_3day_{RUNDATE}.pdf")
    rows, n = dist_table("gated")
    lines = ["| Day | Hour | surface | convect | overshoot | resilient | deep |",
             "|:--|:--|--:|--:|--:|--:|--:|"]
    for valid, h, pct in rows:
        if h in (9, 12, 15, 18):        # daytime rows only, to keep the table compact
            lines.append(f"| {valid} | {h:02d}Z | " + " | ".join(str(x) for x in pct) + " |")
    tbl = "\n".join(lines)
    with open(md, "w") as f:
        f.write(f"""---
title: "Tuscany Pyroconvection -- 3-Day Forecast -- run {RUNDATE} {RUN:02d}Z"
subtitle: "Valid {DAYS[0]} .. {DAYS[2]}. Hybrid: ICON-EU model levels + ICON-2I 2.2 km fuel gate. Method after Castellnou et al. (2022)."
geometry: a4paper, landscape, margin=1.1cm
fontsize: 9pt
---

## 3-day forecast

Three valid days from the single {RUNDATE} {RUN:02d}Z run (day 0, +1, +2), using the day+1/+2
forecast steps in the same ICON-EU + ICON-2I datasets. Atmosphere from ICON-EU native model
levels (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta,
ABL-top RH, shear); the 10 MW/m fuel gate uses ICON-2I's 2.2 km surface fields on the Tuscany
.lcp fuels. See scripts/validation/README.md for the diagnostic validation (radiosondes + ERA5).

The **fuel-gated** map is the expected product; the **potential** map is the atmospheric
upper bound (assumes a pyroCu-capable fire in every cell). Read class *counts* as indicative.

## Fuel-gated (expected)

![gated]({os.path.basename(png_gate)}){{width=100%}}

## Potential (atmospheric upper bound)

![potential]({os.path.basename(png_pot)}){{width=100%}}

## Fuel-gated class distribution -- daytime, % of land cells ({n} land cells)

{tbl}

Classes: 0 surface plume, 1 convection plume, 2 overshooting pyroCu, 3 resilient pyroCu,
4 deep pyroCu/pyroCb. Generated by scripts/forecast_3day_tuscany.py.
""")
    try:
        subprocess.run(["pandoc", md, "-o", pdf, "--pdf-engine=tectonic"], check=True,
                       capture_output=True, timeout=300)
        return pdf
    except Exception as e:
        sys.stderr.write(f"PDF build skipped ({e})\n"); return None


def main():
    os.makedirs(OUT, exist_ok=True)
    for valid in DAYS:
        run_daily(valid)
    png_pot = os.path.join(OUT, f"forecast_3day_potential_{RUNDATE}.png")
    png_gate = os.path.join(OUT, f"forecast_3day_gated_{RUNDATE}.png")
    combined_figure("potential", png_pot)
    combined_figure("gated", png_gate)
    pdf = build_report(png_pot, png_gate)
    print(f"\nOK 3-day forecast -> {OUT}")
    print(f"  {png_pot}\n  {png_gate}" + (f"\n  {pdf}" if pdf else ""))


if __name__ == "__main__":
    main()
