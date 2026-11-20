# SPDX-License-Identifier: AGPL-3.0-or-later
"""Atmospheric time series at the ignition point, 6 h before -> 6 h after the run.

Panels: pyroconvection indices (Continuous-Haines, plume loft factor), convective
moisture (LCL, surface dewpoint depression, mid-level RH), dead fuel moisture
(1/10/100 h -- the *lagged* time-lag model the simulation uses, spun up overnight),
10 m wind (speed + direction), and 2 m relative humidity + T. The fuel-moisture
panel also overlays the FFMC-derived fine fuel moisture and annotates the day's
Canadian FWI System codes (fire-danger context). Source: ICON-2i 2.2 km at
43.93561 N, 11.09635 E. Window shaded = the 14:00-20:00 CEST simulation.
"""
import os, sys, csv, datetime as dt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from _icon2i_provider import build_icon2i, profile_at
import run_fire_11_09_2pm as R
from pyflam import (convective_plume_factor, continuous_haines, inverted_v,
                    lcl_height_m, fwi)
from pyflam.atmosphere import DeadFuelMoistureModel

LAT, LON = R.LAT, R.LON
CEST = dt.timezone(dt.timedelta(hours=2))
UTC = dt.timezone.utc
# window: 6 h before start (12Z) .. 6 h after end (18Z) -> 06Z .. 24Z, hourly
t0 = dt.datetime(2026, 7, 3, 6, tzinfo=UTC)
hours = [t0 + dt.timedelta(hours=h) for h in range(0, 19)]
sim0 = dt.datetime(2026, 7, 3, 12, tzinfo=UTC)
sim1 = dt.datetime(2026, 7, 3, 18, tzinfo=UTC)

prov, ds, levels, files = build_icon2i(dt.date(2026, 7, 3), run=0)

# --- lagged dead fuel moisture: spin the time-lag model up from the 00Z run so the
# 10/100-h classes carry the overnight recovery, matching the simulation. Step it
# hourly from 00Z; keep the values on the plotted hours.
spin0 = dt.datetime(2026, 7, 3, 0, tzinfo=UTC)
_fm = DeadFuelMoistureModel.equilibrium(prov.state_at(LAT, LON, spin0))
lagged = {spin0: dict(m_1h=_fm.m_1h, m_10h=_fm.m_10h, m_100h=_fm.m_100h)}
_t = spin0 + dt.timedelta(hours=1)
while _t <= hours[-1]:
    lagged[_t] = _fm.update(prov.state_at(LAT, LON, _t), 60.0)
    _t += dt.timedelta(hours=1)

# --- Canadian FWI System at solar noon (12Z): fire-danger context + an FFMC-based
# fine fuel moisture to cross-check the dead-fuel model. Dry spell -> rain 0; codes
# start from the spring-startup values (indicative -- DMC/DC need a multi-day drive).
_stn = prov.state_at(LAT, LON, sim0)
FWI = fwi.FWISystem().step(
    temperature=_stn.temperature, relative_humidity=_stn.relative_humidity,
    wind_kmh=fwi.wind_ms_to_kmh(_stn.wind_speed), rain_mm=0.0, month=sim0.month)
print(f"FWI 12Z (spring-startup, rain=0): FFMC {FWI.ffmc:.1f}  DMC {FWI.dmc:.1f}  "
      f"DC {FWI.dc:.1f}  ISI {FWI.isi:.1f}  BUI {FWI.bui:.1f}  FWI {FWI.fwi:.1f}  "
      f"(fine fuel moisture {FWI.fine_fuel_moisture:.1f}%)")

rows = []
for t in hours:
    st = prov.state_at(LAT, LON, t)
    prof = profile_at(prov, levels, st, LAT, LON, t)
    m = lagged[t]
    ch = float(continuous_haines(prof))
    iv, depr, midrh = inverted_v(prof)
    rows.append(dict(
        time=t, local=t.astimezone(CEST),
        wind=float(st.wind_speed), wdir=float(st.wind_direction),
        temp=float(st.temperature), rh=float(st.relative_humidity),
        m1h=float(m["m_1h"]) * 100, m10h=float(m["m_10h"]) * 100, m100h=float(m["m_100h"]) * 100,
        chaines=ch, lcl=float(lcl_height_m(st.temperature, st.relative_humidity)),
        depr=float(depr), midrh=float(midrh),
        plume=float(convective_plume_factor(st, profile=prof))))

x = [r["local"] for r in rows]

# ---- style -------------------------------------------------------------------
INK, MUT, LINE = "#1c1712", "#6f675d", "#e6ddd0"
EMBER, AMBER, RED, TEAL, BLUE = "#e2571f", "#f2a63b", "#a11414", "#2a7f7a", "#3b6fa1"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.edgecolor": "#bbb", "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUT, "ytick.color": MUT})

fig, axes = plt.subplots(5, 1, figsize=(10.5, 13.2), sharex=True)
fig.subplots_adjust(hspace=0.16, top=0.945, bottom=0.06, left=0.085, right=0.9)


def shade(ax):
    ax.axvspan(sim0.astimezone(CEST), sim1.astimezone(CEST), color=AMBER, alpha=0.14, lw=0)
    ax.axvline(sim0.astimezone(CEST), color=EMBER, ls="--", lw=1, alpha=0.7)
    ax.axvline(sim1.astimezone(CEST), color=EMBER, ls="--", lw=1, alpha=0.7)
    ax.grid(True, axis="both", color=LINE, lw=0.7)
    ax.set_axisbelow(True)


def g(k): return [r[k] for r in rows]

# P1 pyroconvection: C-Haines + plume factor
ax = axes[0]; shade(ax)
ax.plot(x, g("chaines"), color=RED, lw=2.2, marker="o", ms=3.5, label="Continuous-Haines")
ax.set_ylabel("Continuous-Haines\n(0–13)", color=RED); ax.set_ylim(0, 10)
ax.tick_params(axis="y", colors=RED)
axb = ax.twinx()
axb.plot(x, g("plume"), color=EMBER, lw=1.8, ls="-", marker="s", ms=3, label="plume loft factor")
axb.set_ylabel("plume loft factor (×)", color=EMBER); axb.set_ylim(0.9, 1.5)
axb.tick_params(axis="y", colors=EMBER)
ax.set_title("Pyroconvection assessment", loc="left", weight="bold", color=INK, fontsize=11)

# P2 convective moisture: LCL + dewpoint depression + mid RH
ax = axes[1]; shade(ax)
ax.plot(x, g("lcl"), color=BLUE, lw=2.2, marker="o", ms=3.5, label="LCL height (m)")
ax.set_ylabel("LCL height (m)", color=BLUE); ax.tick_params(axis="y", colors=BLUE)
axb = ax.twinx()
axb.plot(x, g("depr"), color=AMBER, lw=1.8, marker="^", ms=3.5, label="sfc dewpoint depression (°C)")
axb.plot(x, g("midrh"), color=TEAL, lw=1.5, ls="--", marker="v", ms=3, label="mid-level RH (%)")
axb.set_ylabel("dewpt depr (°C) / mid-RH (%)")
ax.set_title("Convective moisture (inverted-V predictors)", loc="left", weight="bold", fontsize=11)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels()
axb.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left", framealpha=0.85)

# P3 fuel moisture (lagged time-lag model) + FFMC fine fuel moisture cross-check
ax = axes[2]; shade(ax)
ax.plot(x, g("m1h"), color=EMBER, lw=2.2, marker="o", ms=3.5, label="1-h (lagged)")
ax.plot(x, g("m10h"), color=AMBER, lw=2, marker="s", ms=3, label="10-h (lagged)")
ax.plot(x, g("m100h"), color="#8a6d3b", lw=2, marker="^", ms=3, label="100-h (lagged)")
ax.axhline(FWI.fine_fuel_moisture, color=RED, ls=":", lw=1.6,
           label=f"FFMC fine fuel ({FWI.fine_fuel_moisture:.1f}%)")
ax.set_ylabel("dead fuel moisture (%)")
ax.set_title("Dead fuel moisture — time-lag model (overnight spin-up)",
             loc="left", weight="bold", fontsize=11)
ax.legend(fontsize=8.5, loc="upper left", ncol=4, framealpha=0.85)
ax.text(0.985, 0.06,
        f"FWI 12Z: FFMC {FWI.ffmc:.0f} · ISI {FWI.isi:.1f} · "
        f"BUI {FWI.bui:.0f} · DC {FWI.dc:.0f} · FWI {FWI.fwi:.1f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
        color=INK, bbox=dict(boxstyle="round,pad=0.3", fc="#fbf7f0", ec="#d8cdbc"))

# P4 wind speed + direction
ax = axes[3]; shade(ax)
ax.plot(x, g("wind"), color=TEAL, lw=2.2, marker="o", ms=3.5, label="10 m wind speed")
ax.set_ylabel("wind speed (m/s)", color=TEAL); ax.tick_params(axis="y", colors=TEAL)
ax.set_ylim(0, max(g("wind")) * 1.3)
axb = ax.twinx()
axb.scatter(x, g("wdir"), color=MUT, s=14, marker="x", label="direction (° from)")
axb.set_ylabel("wind direction (° from)", color=MUT); axb.set_ylim(0, 360)
axb.set_yticks([0, 90, 180, 270, 360]); axb.tick_params(axis="y", colors=MUT)
ax.set_title("10 m wind", loc="left", weight="bold", fontsize=11)

# P5 RH + temperature
ax = axes[4]; shade(ax)
ax.plot(x, g("rh"), color=BLUE, lw=2.2, marker="o", ms=3.5, label="2 m RH")
ax.set_ylabel("relative humidity (%)", color=BLUE); ax.tick_params(axis="y", colors=BLUE)
ax.set_ylim(0, 100)
axb = ax.twinx()
axb.plot(x, g("temp"), color=RED, lw=1.8, ls="--", marker="s", ms=3, label="2 m T")
axb.set_ylabel("temperature (°C)", color=RED); axb.tick_params(axis="y", colors=RED)
ax.set_title("2 m relative humidity & temperature", loc="left", weight="bold", fontsize=11)

axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=CEST))
axes[-1].xaxis.set_major_locator(mdates.HourLocator(interval=2))
axes[-1].set_xlabel("local time (CEST) — 3 → 4 July 2026")
for ax in axes:
    ax.set_xlim(x[0], x[-1])

fig.suptitle("Calvana fire — ICON-2i atmospheric drivers at the ignition point "
             "(43.93561 N, 11.09635 E)\nsimulation window 14:00–20:00 CEST shaded",
             fontsize=12.5, weight="bold", y=0.985)
out_png = os.path.join(R.OUT, "atmo_metrics_timeseries.png")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print("wrote", out_png)

# CSV
out_csv = os.path.join(R.OUT, "atmo_metrics_timeseries.csv")
with open(out_csv, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["time_utc", "time_cest", "wind_ms", "wind_dir_deg", "temp_C", "rh_pct",
                "m1h_lag_pct", "m10h_lag_pct", "m100h_lag_pct", "continuous_haines",
                "lcl_m", "sfc_dewpt_depr_C", "mid_rh_pct", "plume_factor",
                "ffmc", "isi", "bui", "dc", "fwi", "ffmc_fine_fuel_pct"])
    for r in rows:
        w.writerow([r["time"].strftime("%Y-%m-%d %H:%M"), r["local"].strftime("%Y-%m-%d %H:%M"),
                    f"{r['wind']:.2f}", f"{r['wdir']:.0f}", f"{r['temp']:.1f}", f"{r['rh']:.0f}",
                    f"{r['m1h']:.1f}", f"{r['m10h']:.1f}", f"{r['m100h']:.1f}", f"{r['chaines']:.2f}",
                    f"{r['lcl']:.0f}", f"{r['depr']:.1f}", f"{r['midrh']:.0f}", f"{r['plume']:.3f}",
                    f"{FWI.ffmc:.1f}", f"{FWI.isi:.2f}", f"{FWI.bui:.1f}", f"{FWI.dc:.1f}",
                    f"{FWI.fwi:.2f}", f"{FWI.fine_fuel_moisture:.1f}"])
print("wrote", out_csv)
