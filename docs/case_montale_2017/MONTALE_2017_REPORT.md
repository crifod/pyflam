# Montale (Tobbiana) 2017 wildfire — full pyflam pipeline case study

**Fire:** RTInc2017_495 — Incendio boschivo, loc. **Tobbiana, Comune di Montale (PT)**, Tuscany.
**Ignition:** 16 Jul 2017 ~12:40 local, 43°57′39″N 11°02′25″E (EPSG:3035 X=4404680 Y=2317385).
**Observed:** ~286 ha forest (LAMMA perimeter), ~230 crews / 443 personnel over 9 days.
**Grid:** 10 m, EPSG:3035, 600×600 (6×6 km) centred on the ignition.

This runs the whole pyflam chain on the real fire, at high fidelity, storing every
module's outputs and a detailed report in its own folder.

## Modules & headline results

| # | Module | Folder | Headline |
|---|---|---|---|
| 1 | **Landscape (LCP)** | `landscape/` | 10 m LCP: FBFM40 fuel (fuel-tos-10m) + 10 m DTM (elev 69–1041 m, slope med 35 %) + resampled canopy |
| 0 | **Weather** | `weather/` | La-Ferruccia station, peak burn: **T 28 °C, RH 24 %, wind 5.1 m/s @ NE (39°)**, 1-h dead fuel 4.8 % |
| 2 | **Surface fire behavior** | `fire_behavior/` | Rothermel + Eikonal: **ROS med 18 m/min (90th 52)**, FLI 3.5–14 MW/m, flame 3–6 m |
| 3 | **Crown fire** | `crown_fire/` | Cruz 2005: **48 % of cells active crown**, 3 % passive; crown FLI to ~104 MW/m (capped) |
| 4 | **Spotting** | `spotting/` | Firebrand loft/drift: **12 spot ignitions, max 3.0 km** downwind (SW) of the early front |
| 5 | **Pyroconvection** | `pyroconvection/` | **Real ERA5 sounding**: CAPE 0, ABL 1943 m, LCL 2011 m → **overshooting pyroCu**, C-Haines 8.8, decoupling 8.4× |

## The fire–atmosphere story
The ERA5 profile for 16 Jul 2017 15:00 UTC is the **classic dry-pyroconvection
environment**: a deep (~1.9 km) boundary layer, **near-zero CAPE**, and a **very dry
mid-troposphere** (RH 10–14 % at 650–700 hPa) under a strong upper jet. pyflam
classifies it **overshooting pyroCu**, with Continuous Haines 8.8 and a strong
fire-induced-boundary-layer decoupling — pyroconvection here is a *vertical* (LCL/ABL/
dryness-aloft) problem, not a surface-CAPE one. This matches the fire's known extreme,
crew-intensive behaviour.

## Simulated vs observed
The NE wind drove the fire **SW/downslope**, matching the observed perimeter's
orientation. The free-burning sim (no suppression) reaches the observed **286 ha at
~56 min**; beyond that it keeps spreading because suppression by ~230 crews is not
modelled. The trustworthy comparison is **direction, fire type, and intensity**, not
the free-burn final area. See `reports/`:
- `map_arrival_vs_observed.png` — arrival isochrones + observed perimeter + ignition
- `map_crown_fire_type.png`, `map_fireline_intensity.png`
- `spread_vs_observed.md`

## Method notes & caveats
- **Effective-wind-limit cap:** pyflam has no Rothermel effective-wind-speed limit;
  ~1.7 % of cells (SH/TU/TL fuels) drove the wind factor to a singularity and were
  capped at physical maxima (ROS 122 m/min, FLI 104 MW/m) — the same cap FlamMap
  applies. Median/90th behaviour is unaffected.
- **fireABL magnitude** is a diagnostic (read the decoupling ratio, not the absolute
  16 km — the known over-prediction); ABL from ERA5's assimilated BLH (the coarse
  21-level profile under-resolves the mixed layer for bulk-Ri).
- Surface weather from the local station; profile/CAPE/BLH from ERA5. No suppression.

Reproduce: `build_lcp.py` → `fetch_era5.py` → `run_pipeline.py` → `render_report.py`
(in `scripts/case_montale/`).
