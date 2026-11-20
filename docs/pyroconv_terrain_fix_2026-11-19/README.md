# Pyroconvection maps after the terrain fix — ICON-2I 2026-11-19 00Z

Real-data check of the surface-pressure fix to the pressure-level classifier
(`pyflam_gui.core.pyroconv.surface_layer_terms`; see CHANGELOG `[Unreleased]` / Fixed).
The June 2026-06-26 inputs of the original maps are no longer in either archive, so the
latest available run was used.

| Files | Script | Input |
|---|---|---|
| `pyroconv_map_tuscany_icon2i_fulltuscany_{potential,gated}_2026-11-19.png`, `pyroconv_rasters_native_fulltuscany/` | `tests/pyroconv_icon2i_fullitaly.py` | MISTRAL full-Italy ICON-2I (`fetch_icon2i_mistral`, incl. `PS`) |
| `pyroconv_map_tuscany_icon2i_native_{potential,gated}_2026-11-19.png`, `pyroconv_rasters_native/` | `tests/pyroconv_icon2i_native.py` | Arpae-SIMC `icon_2I.202611190000.grib` (incl. `sp`) |

Grey = nodata (surface within 250 m of 850 hPa). "Potential" is atmosphere only; "gated"
applies the fuel/FLI gate, which blanks everything in mid-November.

12Z, same data, before vs after the fix (`docs/gauge_audit/compare_*.py`):

| | old logic | with surface pressure |
|---|---|---|
| Full Tuscany, 13,920 cells | 13,919 convection plume | 7,984 surface / 4,362 convection / 161 overshooting / 496 resilient / 588 deep / 329 nodata |
| N Tuscany + E-R native, 12,393 cells | almost all surface plume | 10,682 / 601 / 72 / 148 / 461 / 429 nodata |

Hand-checked cells: `docs/gauge_audit/ridge_cell.py` (deep-class ridge cell, genuinely well
mixed; the old logic read it as stable) and `lowland_cell.py` (85 m noon inversion that the
old sea-level ABL datum inflated to 195 m).

Reproduce (from the repo root, with the inputs downloaded):

```bash
RUN_DATE=2026-11-19 ICON2I_DIR=<mistral dir> OUTDIR=<out> LCP=<canopy_tuscany.lcp> python tests/pyroconv_icon2i_fullitaly.py
RUN_DATE=2026-11-19 ICON2I_GRIB=<arpae grib> OUTDIR=<out> LCP=<canopy_tuscany.lcp> python tests/pyroconv_icon2i_native.py
```
