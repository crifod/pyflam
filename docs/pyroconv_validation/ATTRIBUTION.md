# Pyroconvection validation corpus — provenance

Vendored so the §21 analysis in `../graf_vs_pyflam_2026-07-26.md` reproduces from a clone with
no re-download. Every file here is third-party data, redistributed under the terms below. **None
of it is pyflam's own work**, and it must keep its attribution if reused.

Regenerate anything with `PYTHONPATH=src:scripts python scripts/<script>.py` — the scripts read
this directory through `scripts/valdata.py` and take no external inputs.

## Sources

### `sondes/` — 27 ambient radiosonde profiles
Castellnou Ribau et al. (2025), *Atmospheric Measurement Techniques* **18**, 7805–7831, and the
associated dataset (Zenodo 15264835). Radiosondes released beside active wildfires by GRAF
(Grup de Recolzament d'Actuacions Forestals, Bombers de la Generalitat de Catalunya) and
partners. One CSV per labelled sonde–fire pair; the full archive holds 61.

Launch datetimes in `sonde_times.json` are from that paper's Table S1. They are **UTC**, which is
what establishes that the perimeter timestamps are local — see
`scripts/isochrone_firepower.py`.

### `perims/` — perimeter KMZs
Wildfire Data Portal, `wildfiredataportal.eu` (EWED/ODET, EU co-funded). The source of every
growth rate in §21. Timestamps are **local time**; `isochrone_firepower.py` converts them.

### `portal_fires.json`, `portal_all.json`, `fire_labels.json`
Wildfire Data Portal REST catalogue (`/wp-json/wp/v2/fire`), including the
`fire-classification` taxonomy. The portal's class names match the campaign paper's ladder and
the two agree on every fire they share, so `fire_labels.json` and the portal labels are one
source, not two independent ones — **they must not be treated as mutual corroboration.**

### `era5/` — ERA5 pressure-level profiles
Copernicus Climate Change Service (C3S) Climate Data Store, `reanalysis-era5-pressure-levels`.
Generated using Copernicus Climate Change Service information. Neither the European Commission
nor ECMWF is responsible for any use of this data. Small lat/lon boxes only.

* `era5s_<slug>_<date>_<hh>.nc` — at each **sonde launch** hour. Used to extend sonde profiles
  above their apex, which is required for the PFT: a sonde topping out at 3–5 km cannot reach
  the free-convection height or the −20 °C level.
* `era5o_<slug>_<date>_<hh>.nc` — at each out-of-sample fire's **peak-growth** hour. These
  fires have no sounding, so the atmosphere comes from reanalysis alone.

An earlier batch taken at each fire's *start* hour is deliberately **not** vendored: splicing it
onto sondes launched hours later was a real error (§18), and excluding it stops the file-matching
fallback from silently reintroducing it.

## Derived files

`oos_margin.json` and `era5_oos_index.json` are computed, not source — regenerate with
`scripts/oos_margin.py` and `scripts/era5_oos.py`. `era5_oos.py` needs a configured `cdsapi`
client; nothing else here does.

## Scope

27 sonde-columns over 19 sonded fires, plus 8 fires with a label but no sounding: **30 usable
fire-columns, of which 2 are observed pyroCb.** Every statistic in §21 rests on that, and the
effective sample is ~17 independent fires rather than 30 columns. The section states this at
each result; it should not be quoted without it.
