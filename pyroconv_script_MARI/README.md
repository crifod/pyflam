# Pyroconvection reference scripts — by Mari

**Developer: Mari.** The scripts in this folder were developed by Mari, and the folder is
published as Mari provided it, with credit to Mari as the developer.

These are the original scripts that pyflam's pyroconvection classifier was ported from and
validated against. pyflam's code cites them as "the reference (MARI) implementation", for
example in `src/pyflam/atmosphere.py` (`PyroconvThresholds`) and `pyflam_gui/core/pyroconv.py`.
They are kept here unchanged, so every threshold and design choice in the port can be traced
back to its source.

## Contents

| Folder | Script | What it does |
|---|---|---|
| `pyroconvection_era5/` | `pyroconvection_era5_v1_1_regions.py` (v1.1) | Atmospheric predisposition to pyroconvection from **ERA5** reanalysis (hindcast), following Castellnou et al. (2022). |
| `pyroconvection_icon2i/` | `pyroconvection_icon2i_meteohub_v1_2.py` (v1.2) | The same diagnostic on **ICON-2I** forecasts from the MeteoHub / ItaliaMeteo open-data archive. |
| `pyroconvection_wrf15/` | `pyroconvection_wrf15_v1_7_4.py` (v1.7.4) | The same diagnostic on the operational **WRF 1.5 km** ItaliaMeteo-CIMA-LaMMA forecast, with a quality index per cell; `docs/` holds the accompanying technical report (*Rapporto tecnico sulla predisposizione atmosferica alla piroconvezione con WRF 1.5 km*, PDF/DOCX) and a background research report. |

Each script classifies every grid cell into the Castellnou et al. (2022) ladder: 0 surface
plume, 1 convective plume, 2 overshooting pyroCu, 3 resilient pyroCu, 4 deep pyroCu / pyroCb.
To do so it reconstructs AGL profiles, estimates the ABL with the bulk Richardson number,
and computes the LCL, ABL stability, the θ gradient above the ABL, the height of maximum wind
shear and the RH at the ABL top. Each subfolder has its own README (in Italian) and
`requirements_*.txt`.

## Licence

**These files are not covered by pyflam's AGPL-3.0 licence.** Mari has not attached a
licence to them, so all rights stay with Mari. If you want to reuse, modify or redistribute
them beyond reading them alongside pyflam, contact Mari through the pyflam maintainer.

## Citing

When you use or discuss these scripts, credit **Mari** as the developer, e.g.:

> Mari (2026). *Pyroconvection predisposition scripts for ERA5, ICON-2I and WRF 1.5 km*
> (ERA5 v1.1, ICON-2I v1.2, WRF v1.7.4). In: pyflam repository,
> https://github.com/crifod/pyflam/tree/main/pyroconv_script_MARI

The scientific framework is Castellnou, M., et al. (2022), on which the scripts build.
