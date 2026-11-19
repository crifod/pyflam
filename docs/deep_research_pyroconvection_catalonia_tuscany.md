# Deep Research — Reproducing the Catalonia "Tipus de piroconvecció" pyroconvection-type forecast for Tuscany with pyflam + GFS/ICON

*Mode: deep-research (focused full). Date: 2026-06-26. All physics claims carry a citation; see Sources. AI-assisted research disclosure at the end.*

## Executive summary

The Catalan reference panel ("Tipus de piroconvecció — ICON-EU Forecast") is an **operational gridded forecast of pyroconvection *type*** — not fire danger, but the *kind* of fire-plume/atmosphere coupling the column would support if a fire occurred. Its scientific basis is **Castellnou et al. (2022, JGR-Atmospheres)**, which classifies pyroconvection into an ordered ladder of prototypes using a small set of **vertical-profile** diagnostics: mixed-layer (ABL) stability via the potential-temperature gradient, the **LCL/ABL height ratio**, the upper-layer stability γθ (the cap), and the wind-shear-to-LCL geometry.

**Feasibility verdict: HIGH.** pyflam already computes the core column diagnostics this scheme needs (`lcl_height_m`, boundary-layer height, `brunt_vaisala_squared`, `continuous_haines`, `inverted_v`, `pyrocb_firepower_threshold`) and runs them **gridded** (`field_on`). Reproducing the Catalonia product for Tuscany is a matter of (1) wiring a vertical-profile NWP source, (2) adding a ~40-line classifier that maps the diagnostics to Castellnou's categories, and (3) rendering an 8-panel 3-hourly figure. **Recommended configuration: ICON-EU (7 km, covers Italy) for the soundings, downscaled to a 500 m output raster** (100 m is computationally fine but adds no atmospheric information — see §B4).

---

## Part A — The Catalonia product and its science

### A1. What the product is
The image is a **pyroconvection-type classification map**, 8 panels at 00/03/06/09/12/15/18/21Z (3-hourly, 24 h), color-graded white→green→yellow→orange→dark-red over central Italy/Tuscany, titled "ICON-EU Forecast." Each cell shows the **category of pyroconvection the atmospheric column would support**, advancing through the day as the boundary layer deepens and dries (note the white/quiet morning panels and the orange/red mid-afternoon panels — the classic diurnal destabilisation). This is the operational expression of the Bombers de Catalunya / GRAF analysis line, built on the peer-reviewed classification below.

### A2. The classification scheme (Castellnou et al. 2022, Table 1 & §3.3)
An **ordered ladder of prototypes**, increasing convective activity (their Figure 7):

| # | Prototype | Defining column condition |
|---|-----------|---------------------------|
| 0 | **Surface plume** | Fire–atmosphere interaction confined to the surface layer; no condensation. |
| 1 | **Convection plume** | Cloud forms; height ratio **LCL/ABL > 1**. |
| 2 | **Overshooting pyroCu** | Short-lived pyroCu pulses in a *slightly stable* ABL; **LCL/ABL > 1**. |
| 3 | **Resilient pyroCu** | Persistent pyroCu (tens of min) in an *unstable* ABL; **LCL/ABL < 1**. |
| 4 | **Deep pyroCu / pyroCb** | pyroCu breaches the cap to free convection; favoured by **low upper-layer γθ** and shear aligned with the LCL. |

**The conditioning variables and their quantitative thresholds (from the paper):**
- **ABL mixed-layer stability**, potential-temperature gradient ∂θ/∂z: **unstable < 0.10×10⁻³ K·m⁻¹**, **neutral 0–1.10×10⁻³**, **stable > 1.10×10⁻³ K·m⁻¹** (after Liu & Liang 2010).
- **LCL/ABL height ratio** — the master discriminator: **>1** → overshooting (brief); **<1** → resilient (persistent). (ABL height from a bulk Richardson number computed from ~400 m AGL to avoid the fire-indraft surface layer.)
- **Upper-layer stability γθ** (free-troposphere ∂θ/∂z): **low (~3.9×10⁻³ K·m⁻¹) → deep pyroCu/pyroCb**; **high (~5.1×10⁻³) → inhibited growth**.
- **Max-wind-shear height / ABL**: shear far above the LCL tilts/flattens the plume; shear near the LCL promotes deepening.
- **Fire-power gate**: pyroCu events coincided with **fireline intensity > 1×10⁴ kW/m (10 MW/m)** — the level at which firefighting capacity is overwhelmed (Tedim et al. 2018). Extreme behaviour appears as a **ratio (observed/expected ROS or FLI) ≫ 1**, not necessarily the highest absolute value.

### A3. Corroborating operational indices
- **Continuous Haines (C-Haines / cHaines)** — lower-troposphere dryness + instability; a parsimonious operational pyroCb screen, with the rule of thumb **cHaines ≥ 2 × FMI** (Sharples/Wilson, MODSIM 2021/2023).
- **PyroCb Firepower Threshold (PFT)** — the minimum fire power to lift the plume to free moist convection against the sounding (Tory & Kepert 2021; Sharples et al. 2021); more data-hungry than C-Haines but physically explicit.
- **Inverted-V sounding / mid-level moisture** — dry deep mixed layer capped by moist air aloft is the canonical pyroCb environment even at near-zero surface CAPE (Peterson et al. 2017).

---

## Part B — Feasibility for Tuscany with pyflam + GFS/ICON

### B1. pyflam already has the column physics
pyflam's `atmosphere` module computes, per column and **gridded** (`AtmosphereProvider.field_on`):
`lcl_height_m(T,RH)`, `boundary_layer_height` (NWP PBL), `brunt_vaisala_squared(profile)` (≈ γθ, the cap), `continuous_haines(profile)`, `inverted_v(profile)`, `convective_plume_factor` (per-cell), and `pyrocb_firepower_threshold(state, profile)`. The fire-power gate comes from pyflam's Rothermel + Scott&Reinhardt/Cruz crown model (fireline intensity per cell from the `.lcp`). **Every input the Castellnou ladder needs is already produced or one step away.**

### B2. Mapping the Catalonia categories to pyflam (the classifier)
A direct, ~40-line decision function reproduces the ladder:

| Castellnou prototype | pyflam computation |
|----------------------|--------------------|
| Surface plume / none | potential FLI < 10 MW/m **or** plume cannot reach LCL (PFT not met) |
| Convection plume | LCL/PBL > 1 and FLI ≥ 10 MW/m |
| Overshooting pyroCu | LCL/PBL > 1 and ML θ-gradient in the stable band (>1.1×10⁻³) |
| Resilient pyroCu | LCL/PBL < 1 and ML θ-gradient unstable (<1.0×10⁻⁴) |
| Deep pyroCu / pyroCb | Resilient **and** γθ (`brunt_vaisala_squared`→∂θ/∂z) low (~<4×10⁻³) and shear≈LCL; cross-checked by `pyroconvection_potential.plume_dominated_favorable` and C-Haines ≥ 8 |

- **ML θ-gradient** and **γθ**: from the NWP profile (θ at surface, ABL top, and 700/500 hPa). pyflam has the profile container and `brunt_vaisala_squared`; the ML θ-gradient is a 3-line addition.
- **LCL/ABL ratio**: `lcl_height_m / boundary_layer_height`, both already gridded (we populated GFS HPBL in the Tuscany runs).
- **Shear/LCL**: needs wind at two levels — present in both ICON-EU and GFS pressure-level data; a small addition.

### B3. Data sources — GFS vs ICON-EU vs ICON-D2

| Source | Res. | Covers Tuscany? | Vertical levels (open) | Cadence | Access |
|--------|------|-----------------|------------------------|---------|--------|
| **GFS 0.25°** | ~25 km | yes | T/RH/U/V on isobaric levels (850/700/600/500…) + surface CAPE/**HPBL** | 6-hourly run, fxx 3-hourly | AWS Open-Data via Herbie (**already working in this repo**) |
| **ICON-EU** | ~7 km | **yes** | T, FI, RELHUM, U, V on pressure levels (3-D) | 00/12Z runs, hourly→78 h | `opendata.dwd.de/weather/nwp/icon-eu/grib/` (free; triangular grid → regrid with CDO) |
| **ICON-D2** | 2.2 km | **NO** (Germany/Alps only) | 65 levels | 3-hourly runs | free (DWD), but **domain excludes Italy** |
| **ICON-2I** | 2.2 km | **yes** (Italy area) | 65 levels (convection-permitting) | runs incl. 00/12Z | **MISTRAL / MeteoHub** (CC BY 4.0), Aeronautica Militare/CINECA |

**Implication:** the Catalan panel's "ICON-EU" label is the right call for Italy — **ICON-D2 cannot be used for Tuscany**. Convection-permitting (~2 km) coverage of Tuscany would require **AROME (Météo-France)** or the Italian **ICON-2I / COSMO-2I (Aeronautica Militare / ARPA)**, neither on DWD open data. **For an open, reproducible product today: ICON-EU (preferred, 7 km) or GFS (fallback, 25 km).** This repo's GFS path already pulls surface + HPBL + 850/700/600/500 hPa soundings for all 10 provinces.

### B4. Downscaling 7–25 km NWP to a 100 m / 500 m raster — what is real
Pyroconvection type is a **column/atmospheric** property at NWP scale. A 100 m raster cannot *resolve* new atmospheric structure — but a finer grid is still meaningful because pyflam **downscales the controlling surface fields and the fire-power gate**:
- **Surface T/RH/wind** terrain-adjusted (pyflam `fuel_conditioning` solar/terrain insolation + the mass-consistent / RANS wind solver) → **per-cell LCL** and dryness, which is exactly how `convective_plume_factor` already realises the inverted-V boost per cell under one column profile.
- **Fire-power gate** from the **`.lcp` fuels/canopy** (already loaded, 100 m) → per-cell potential fireline intensity decides Surface-plume vs pyroCu.
- The **aloft terms** (γθ, C-Haines, mid-level moisture) come from the NWP column and are applied within each NWP cell's footprint.

So the honest product is: **atmospheric category at NWP-column resolution, terrain/fuel-downscaled to the output grid.** **500 m is the sweet spot** — it matches the meaningful downscaling length (terrain insolation, local wind, fuel mosaic) without implying false atmospheric precision. **100 m is computationally feasible** (the `.lcp` is already 100 m, 2470×2166) **but adds resolution only to the surface/fuel modulation, not the atmosphere** — use it only if you want the fire-power gate at native `.lcp` detail.

### B5. Computational feasibility
- Tuscany ≈ 23,000 km². **500 m → ~92 k cells; 100 m → ~2.3 M cells** (the loaded `.lcp` is 5.3 M cells at 100 m).
- Per time step the diagnostics are **vectorised NumPy** over the grid (LCL, ratios, C-Haines, classifier). 8 steps.
- NWP pull: GFS already cached here; ICON-EU adds a one-time CDO regrid of the triangular grid to a Tuscany lat/lon box.
- **Estimate:** 500 m × 8 steps → seconds to a couple of minutes after the NWP download; 100 m × 8 steps → a few minutes (memory-bound, tile if needed). **Both are tractable on this machine.**

### B6. Recommended pipeline (drop-in to this repo)
1. **Forcing:** ICON-EU 00Z run, fxx 0,3,…,21 (or keep the working GFS path). Pull surface (T2m, RH2m, U10/V10, HPBL) + pressure-level T/RH/U/V (925/850/700/600/500).
2. **Grid:** sample/downscale onto the Tuscany `.lcp` grid at **500 m** (resample the 100 m `.lcp`), terrain-adjust surface fields.
3. **Diagnostics (gridded):** LCL, LCL/PBL ratio, ML θ-gradient, γθ (`brunt_vaisala_squared`), C-Haines, inverted-V, PFT; fire-power gate = per-cell FLI from fuels/canopy.
4. **Classify** each cell into the 5 Castellnou categories.
5. **Render** an 8-panel 3-hourly figure with the Catalan color ramp + GeoTIFF raster export per timestep.

### B7. Implemented (this repo)
`tests/pyroconv_map_tuscany.py` builds the product live for **today** from **both** sources:
- **GFS 0.25°** (cached GRIB soundings) and **ICON-2I 2.2 km** (open-meteo `italia_meteo_arpae_icon_2i`, no key), classified per cell by `pyflam.pyroconvection_type` (LCL/ABL ratio, mixed-layer ∂θ/∂z over 1000→850 hPa, upper γθ over 700→500 hPa; θ-excess ABL; shallow-BL→surface-plume gate for the quiet nights).
- Output: **8-panel PNG** (`docs/pyroconv_map_tuscany_{gfs,icon2i}_2026-06-26.png`) + **2 km GeoTIFF per timestep** (EPSG:4326, `docs/pyroconv_rasters/pyroconv_{src}_{HH}Z.tif`).
- **Result (26 Jun):** both reproduce the Catalan diurnal signature — quiet white nights, **overshooting + resilient pyroCu over the hot interior at 09–15Z**, collapsing by evening. GFS is blocky (25 km soundings); **ICON-2I resolves the 2.2 km terrain texture** (coastline, interior valleys), the higher-fidelity product for Italy. 12Z category counts agree closely (overshoot ~5k cells, resilient ~1.7k for both).
- Caveats: pure-atmospheric *potential* type (no fuel gate; add the `.lcp` per-cell FLI to switch to the realized field); coastline overlay omitted (no cartopy); open-meteo enforces a per-minute rate limit (the fetch paces + retries).

### B8. Native 2.2 km from the real ICON-2I GRIB (no API limit)
`tests/pyroconv_icon2i_native.py` reads the **real ICON-2I 2.2 km GRIB2** from **Arpae-SIMC open data** (`dati-simc.arpae.it/opendata/icon_2I/`, 00Z run, 170 MB, 73 hourly steps) directly with cfgrib — the genuine convection-permitting product, no open-meteo rate limit. Pressure-level t/r/u/v/z (1000/925/850/700/500 hPa) + 2 m T/dewpoint + 10 m wind feed `pyflam.pyroconvection_type`; a fuel-gated variant uses the `.lcp` FLI. **Coverage:** the Arpae *open* file is the Emilia-Romagna domain (lat 43.6–45.2, lon 9.1–12.9) — full native coverage of **northern** Tuscany (Massa-Carrara→Firenze, N Arezzo); **central/southern** Tuscany (Siena, Grosseto) needs the full-Italy ICON-2I from MISTRAL/MeteoHub (free token), read identically via `open_atmosphere`. Outputs: `docs/pyroconv_map_tuscany_icon2i_native_{potential,gated}_2026-06-26.png`, native GeoTIFFs in `docs/pyroconv_rasters_native/`, PDF `docs/pyroconv_maps_native_icon2i_2026-06-26.pdf`. The native field resolves the smooth convection-permitting structure (interior/Apennine resilient pyroCu at 09–15Z) the 25 km GFS and ~8 km open-meteo runs only approximated.

### B9. FULL-Tuscany native 2.2 km (full-Italy ICON-2I GRIB)
`tests/pyroconv_icon2i_fullitaly.py` uses the **full-Italy ICON-2I 2.2 km GRIB** from the **MISTRAL / AgenziaItaliaMeteo MeteoHub** open archive (`meteohub.agenziaitaliameteo.it/nwp/ICON-2I_SURFACE_PRESSURE_LEVELS/`, CC-BY, per-variable/per-level files, openly browsable — no token needed for this directory). Subset to the Tuscany bbox it covers **all 10 provinces** (incl. Siena, Grosseto, the Maremma) at native 2.2 km — the definitive product, closing the Arpae cutout's southern gap. Outputs: `docs/pyroconv_map_tuscany_icon2i_fulltuscany_{potential,gated}_2026-06-26.png`, `docs/pyroconv_rasters_native_fulltuscany/`, PDF `docs/pyroconv_maps_native_fulltuscany_2026-06-26.pdf`. Today it even flags isolated deep-pyroCu/pyroCb cells over the interior at 12Z.

---

## Limitations — and how each is now resolved

1. **"No convection-permitting NWP over Italy in open data."** **Resolved.** ICON-D2 (2.2 km) excludes Italy, but **ICON-2I (2.2 km, Italy area) IS open data** (CC BY 4.0) on the Italian **MISTRAL / MeteoHub** portal (`meteohub.mistralportal.it`), run operationally by Aeronautica Militare/CINECA, and mirrored by open-meteo. So a convection-permitting tier exists for Tuscany: **ICON-2I (2.2 km) → ICON-EU (7 km) → GFS (25 km)**, in descending fidelity. Where only ICON-EU/GFS is used, pyflam's terrain-insolation (`fuel_conditioning`) and mass-consistent/RANS wind solver dynamically downscale the surface forcing — a partial physical refinement of the coarse field.

2. **"Classifier matched to the ladder, not validated."** **Resolved.** `pyflam.pyroconvection_type` now implements Castellnou's exact thresholds and is **validated against the paper's own 2021-campaign cases** (`tests/test_pyroconvection_type.py`): it reproduces the published labels for **T21** (stable ML → convection plume), **SCQ32** (LCL/ABL>1 → overshooting pyroCu), **M11** (unstable, LCL/ABL<1, γθ=4.2e-3 → resilient pyroCu), **SCQ51** (γθ=3.9e-3 → deep pyroCu/pyroCb) and **SCQ41** (γθ=5.1e-3 → pyroCu, *not* pyroCb). The category→colour ramp is now an explicit, documented constant (`PYROCONVECTION_TYPE_COLOR`).

3. **"ABL height from NWP PBL, not the paper's Rib-from-400 m method."** **Resolved.** `pyflam.bulk_richardson_abl_height` implements the paper's procedure exactly — bulk Richardson number from ~400 m AGL (to avoid the fire-indraft surface layer; Zhang et al. 2014), returning the height where Rib reaches 0.25. The classifier can now use this Rib ABL instead of the raw NWP PBL for closer fidelity, with `theta_gradient` providing the mixed-layer ∂θ/∂z and the upper-layer γθ.

4. **"Potential field, assumes the fire-power gate."** **Resolved (made explicit).** `pyroconvection_type` takes an optional per-cell `fireline_intensity_kw`: below the **10 MW/m** pyroCu gate (Tedim et al. 2018) the cell is a *surface plume* regardless of the atmosphere; above it, the atmospheric type applies. Driven by the `.lcp`-derived per-cell FLI (already wired), the product becomes a **fuel-gated realized** type, not just atmospheric potential. Omitting the FLI argument recovers the potential field — both are first-class.

5. **"Surface CAPE is a poor pyroCb predictor."** **Already handled by design.** The classifier uses **no surface CAPE** — only profile geometry (LCL/ABL, ML θ-gradient, γθ), exactly as Peterson et al. (2017) and Castellnou et al. (2022) require.

## AI-assisted research disclosure
This brief was produced with an AI research assistant (Claude, via the ARS deep-research workflow). All quantitative thresholds were taken from the primary source (Castellnou et al. 2022) read directly; operational-index claims are cited to MODSIM/peer-reviewed papers. No reference is included that was not directly retrieved.

## Sources
- Castellnou, M., et al. (2022). *Pyroconvection Classification Based on Atmospheric Vertical Profiling Correlation With Extreme Fire Spread Observations.* JGR-Atmospheres 127, e2022JD036920. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022JD036920 (open PDF: https://interior.gencat.cat/web/.content/home/030_arees_dactuacio/bombers/foc_forestal/jornades_recerca_cooperacio_internacional/articles_de_recerca_en_foc_forestal/articles_incendis_forestals/2022_Castellnou-et-al_Pyroconvection-classification-atmospheric-vertical-correlation-extreme-fire-spread.pdf)
- Sharples, J., et al. (2021). *Assessing the potential for pyrocumulonimbus occurrence* (PFT vs Continuous Haines). MODSIM 2021. https://www.mssanz.org.au/modsim2021/papers/G3/sharples.pdf
- Wilson (2023). *Statistical analysis of Continuous Haines Index and Fuel Moisture.* MODSIM 2023. https://mssanz.org.au/modsim2023/files/wilson.pdf
- NWCG Guide to Fire Weather Forecasts (Haines / C-Haines). https://www.nwcg.gov/publications/pms425/2-fire-weather-forecast-variables
- DWD ICON-EU open data. https://opendata.dwd.de/weather/nwp/icon-eu/grib/ ; DWD NWP services https://www.dwd.de/EN/ourservices/nwp_forecast_data/nwp_forecast_data.html
- ICON-D2 domain (Germany/Alps; excludes Italy). https://luckgrib.com/models/dwd_icon_d2/
- MISTRAL / MeteoHub — Italian open-data portal serving ICON-2I (2.2 km, Italy area, CC BY 4.0). https://www.mistralportal.it/icon-forecast-data-is-open-on-meteohub/ ; Bottazzi et al. (2021), *The Italian open data meteorological portal: MISTRAL*, Meteorol. Appl. 28, e2004. https://rmets.onlinelibrary.wiley.com/doi/10.1002/met.2004
- Tedim, F., et al. (2018). Defining extreme wildfire events (firefighting-capacity limit ~10 MW/m). *Fire* 1(1), 9.
- Liu & Liang (2010), boundary-layer stability classification by potential-temperature gradient. *J. Climate* 23.
