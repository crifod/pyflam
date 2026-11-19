# Wildfire in-situ sounding datasets — a survey to augment the fireABL calibration

**Mode:** deep-research, quick brief · **Date:** 2026-07-22
**Goal:** find publicly available datasets worldwide with **real in-situ atmospheric soundings** (radiosondes / dropsondes / mobile balloons) recorded during or near wildfires, to broaden the pyflam fireABL / pyroconvection calibration beyond the GRAF/WUR Iberian set.

**Verification standard:** every dataset below is confirmed against its data-repository page or a peer-reviewed source. Request-only / non-open datasets are flagged; anything I could not confirm to exist is excluded rather than listed as "uncertain."

---

## Executive summary

The **single highest-value, lowest-effort** augmentation already exists and matches our pipeline byte-for-byte: the GRAF group's **2025 in-plume radiosonde release** (Zenodo 15264835) — 6 more Catalan fires, open, same CSV format as `zenodo_6433389`. Beyond Iberia, real fire soundings exist (RaDFIRE, CSU-MAPS, FASMEE, WE-CAN) but the **observed-fireABL/plume-top pairing** that calibration needs is rare — only GRAF and the SJSU (Clements/Lareau) campaigns provide it directly, and the SJSU data is **request-only**. Operational archives (IGRA) offer global scale but not co-located, in-plume soundings.

## Ranked, annotated list

### Tier A — open, fire-paired soundings, directly ingestible

**1. GRAF In-Plume Radiosonde Fireline Observations** — *top pick*
Castellnou Ribau et al. (2025), companion to the AMT paper. **6 fires** (Casablanca III ×2 dates, Granja Escarp, Martorell, Rojals, Santa Ana, Selva del Camp), **21 CSV soundings** (in-plume + environmental `raw_flight_history`), fire perimeters (KML/KMZ). CC-BY-4.0, 152.9 MB, open.
- Zenodo: https://doi.org/10.5281/zenodo.15264835
- **Same group, same format** as our current data → drop-in. 5 new Catalan fires (Martorell overlaps). Effort: hours.

**2. WE-CAN 2018 radiosondes** (NCAR/EOL)
Western-US wildfire-smoke campaign (Boise, Jul–Aug 2018; C-130 sampled 21 wildfire plumes over 16 flights). Includes **NCAR/EOL Quality-Controlled MW41 Radiosonde Data (NetCDF)** — T/RH/wind vs height. Open.
- Project archive: https://data.eol.ucar.edu/project/WE-CAN · MW41 dataset: https://data.eol.ucar.edu/dataset/547.001
- Non-Mediterranean, open. **Gap:** soundings need pairing to fire behaviour (aircraft/fuel assessments exist separately; no ready observed-fireABL column).

### Tier B — real fire soundings, access = request/contact

**3. RaDFIRE** — Rapid Deployments to Wildfires Experiment (SJSU Fire Weather Research Lab, Clements)
**21–26 western-US wildfires, 2013–2018**, radiosondes (potential-temperature profiles, inferred ABL growth) + scanning Doppler lidar + airborne radar + in-plume aircraft (updraft cores to 58 m/s). The premier in-fire sounding campaign.
- BAMS 2018: https://journals.ametsoc.org/bams/article/99/12/2539 · Lab: https://www.fireweather.org/radfire
- **Not openly archived — contact the lab.** Highest-value out-of-region expansion *if* obtainable (has plume observations to pair with soundings).

**4. CSU-MAPS / Lareau & Clements pyroCu-pyroCb soundings**
Mobile radiosondes launched adjacent to fires + Doppler lidar; key finding that the **convective condensation level (CCL), not the LCL**, best predicts pyroCu base height — directly relevant to our in-plume-vs-ambient-LCL question.
- ACP 2016: https://acp.copernicus.org/articles/16/4005/2016/
- **Data "available on request from the authors."**

**5. FASMEE** — Fire and Smoke Model Evaluation Experiment (USFS)
Large prescribed-fire campaigns (Western/Southwest/Southeast) explicitly designed as **model-evaluation datasets**: radiosondes + micrometeorology + plume dynamics (lidar/radar) + fuels/fire-behaviour. The best-instrumented fire–atmosphere pairing in the US.
- Project: https://research.fs.usda.gov/pnw/centers/fasmee · Plan (Atmosphere 2019): https://doi.org/10.3390/atmos10020066
- **No single open DOI found;** access via USFS PNW / FRAMES (https://www.frames.gov/catalog/20954). Prescribed (not wild) fires.

### Tier C — operational soundings near fires (scalable, not co-located)

**6. IGRA v2** (NOAA NCEI) — Integrated Global Radiosonde Archive
Global, quality-controlled radiosonde T/RH/wind + stability parameters, all continents, open. Use the **nearest-station sounding on a fire day** — the same "closest profile" logic as our ERA5 approach (e.g. Black Summer via Wagga Wagga/Melbourne; Fort McMurray 2016).
- https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00975 (also IEM Rawinsonde Archive; University of Wyoming upper-air).

**7. Global merged PBL-height dataset** (ESSD 2024)
~180 M profiles / 370 stations, radiosonde + ERA5 + GLDAS, 2011–2021, 3-hourly, with derived boundary-layer heights. Not fire-specific — useful for ABL-diagnostic climatology/validation.
- https://essd.copernicus.org/articles/16/1/2024/

### Tier D — aircraft profiles (not balloon soundings)

**8. FIREX-AQ 2019** (NASA ASDC / NOAA CSL)
DC-8 / ER-2 / Twin Otter in-plume vertical profiles over North-American fires, open. Rich plume/chemistry/height data but **aircraft, not radiosondes** — different ingestion, weaker as a clean T/RH/wind sounding.
- https://csl.noaa.gov/projects/firex-aq/ · data via NASA Langley ASDC.

## Gaps identified

- **Australia (Black Summer) & Canada (Fort McMurray, 2023–24):** the most spectacular pyroCb on record, but documented mainly by **satellite + operational sondes** — no dedicated open *in-plume* sonde dataset. Must be assembled from IGRA nearest-station + satellite fire-radiative-power/plume-top.
- **The binding constraint is not soundings — it is the *paired observed fireABL / plume-top*.** Outside GRAF, only RaDFIRE/CSU-MAPS provide it directly, and both are request-only. WE-CAN/IGRA give soundings but require separate fire-behaviour linkage.
- FASMEE is prescribed fire (lower intensity), so less useful for the *extreme* pyroconvective end the fireABL targets.

## Recommendations for pyflam

1. **Now (days):** ingest **Zenodo 15264835** — 5 new Catalan fires through the existing GRAF-format pipeline; immediately grows the sonde set and the Stage-A/multi-fire calibration.
2. **Near-term (weeks):** request **RaDFIRE + CSU-MAPS** soundings from the SJSU lab — the highest-value non-Mediterranean expansion, with plume-top observations for the fireABL magnitude check.
3. **Scalable:** build an **IGRA-nearest-station** loader for documented pyroCb days (Black Summer, Canada 2023–24), paired with satellite plume-top — a global, looser validation set that stress-tests the diagnostic outside the Mediterranean.

## Limitations & AI disclosure

Web-search-based reconnaissance, not exhaustive database queries; some datasets' exact sounding counts and access terms are unconfirmed (flagged inline) and request-only holdings were not inspected. Coverage of non-English and pre-2013 campaigns is thin. AI-assisted search; every dataset listed was verified against a data repository or peer-reviewed source before inclusion, per the deep-research verification standard.
