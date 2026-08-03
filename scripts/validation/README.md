# Pyroconvection: radiosonde validation of the coarse-resolution diagnostics

Why the daily ICON-2I pyroconvection product measures mixed-layer stability the way it
does. These scripts are the evidence; they are reproducible from public data.

## The question

The ICON-2I open-data archive publishes **5 usable pressure levels** below 500 hPa
(1000/925/850/700/500). The Castellnou et al. (2022) classification ladder conditions on
the **mixed-layer dtheta/dz**, with a "stable" threshold of 1.1e-3 K/m derived from
high-resolution radiosondes. Does that threshold — and that measurement — survive the
transfer to a 5-level model column?

## Method

Same air column, twice. Take an IGRA radiosonde (full resolution), compute the truth, then
subsample the identical column onto ICON-2I's exact 5 pressure levels and recompute. Model
error is held out entirely, so what is left is the vertical-resolution effect alone.

* **Truth** = least-squares dtheta/dz *within* the parcel mixed layer, at full resolution.
* Stations: **ITM00016144** S. Pietro Capofiume (inland Po valley — the closest analogue to
  inland Tuscany) and **ITM00016245** Pratica di Mare (coastal Rome — a sea-breeze regime,
  included deliberately as a contrast).
* Sample: all JJA 12Z soundings in the IGRA period of record (n ≈ 2,835 usable).

```bash
# data (IGRA v2, NOAA NCEI, public)
B=https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por
for S in ITM00016144 ITM00016245; do
  curl -sO $B/${S}-data.txt.zip && unzip -oq ${S}-data.txt.zip -d por_$S
done

PYTHONPATH=../../src python sounding_degradation.py      # ABL + gradient bias
PYTHONPATH=../../src python ml_metric_calibration.py     # which metric reproduces truth
PYTHONPATH=../../src python ml_metric_levels_in_ml.py    # how many levels are IN the ML
```

## What we found

**1. The Rib ABL is sound when it is resolved.** At full resolution it agrees with the
parcel mixing depth to within ~100 m. Vertical resolution alone makes it *shallower*, not
deeper (−643 m inland) — so "coarse levels give an over-deep ABL" was never the mechanism.

**2. The mixed-layer dtheta/dz is NOT measurable from this archive.** A 5-level column
contains **0–2 model levels inside the mixed layer**:

| levels inside the ML | inland | coastal |
|:--|--:|--:|
| 0 | 3.4% | 24.8% |
| 1 | 26.5% | 67.7% |
| 2 | 61.6% | 6.0% |
| 3+ | **8.5%** | 1.6% |

and the lowest of them sits in the superadiabatic surface layer. Every direct estimate
fails: fitting the gradient across the in-ML levels scores **Youden J ≈ 0.00** (no skill),
and taking theta at the *Rib ABL top* instead folds in the entrainment jump Delta-theta —
which Castellnou et al. (2022 §2.1.1) treat as a state variable **separate** from the
gradient the ladder conditions on (mixed-layer slab model; Vilà-Guerau de Arellano 2015).

**3. So the shipped diagnostic is an honest proxy, not a measurement.** `surface_to_parcel`
is identically `0.5 / depth` (the parcel top is *defined* as theta_sfc + 0.5 K), i.e. it is
the **parcel mixing depth restated** — at the 1.1e-3 threshold it says exactly
*"well-mixed layer ≥ 455 m"*. It is the default only because it is the sole candidate with
any skill:

| coarse metric | J inland | J coastal | r vs truth |
|:--|--:|--:|--:|
| `surface_to_abl` (gradient to Rib top) | 0.14 | 0.42 | +0.16 |
| **`surface_to_parcel`** (shipped) | **0.29** | **0.55** | **+0.50** |
| fit across in-ML levels | 0.00 | 0.06 | +0.26 |

**4. It over-flags — do not read the class counts quantitatively.** Inland specificity is
**0.29**: of columns that are truly ML-stable it still calls ~71% pyroCu-capable
(sensitivity 0.99, so it rarely *misses* a capable column). The ML-stability gate is a
screening filter, not a calibrated one.

Fixing this needs more vertical levels than the Italian 2.2 km open data publishes — which
is exactly what the **hybrid** path adds (below).

## 5. The hybrid: ICON-EU model levels fix it

ICON-2I open data has no model levels (MeteoHub publishes only `..._SURFACE_PRESSURE_LEVELS`,
and `ICON_2I_RUC` carries the same 6 pressure levels). **ICON-EU** (DWD open data) is coarser
horizontally (6.5 km vs 2.2 km) but publishes the native 74 model levels + `HHL` heights.
Degrading the same soundings onto the real ICON-EU model-level heights (decoded from HHL over
Tuscany: 10, 42, 94, 164, 249, 348, 461, 586, 723, 872, 1033, 1205, 1388 m …):

| grid | levels in ML | ML-grad bias | MAE | Youden J | **Rib ABL bias** |
|:--|--:|--:|--:|--:|--:|
| ICON-2I (5 p-lev) | 2.1 | −4.7e-4 | 4.8e-4 | 0.00 | **−696 m** |
| ICON-EU (20 p-lev) | 4.3 | −3.8e-4 | 4.2e-4 | 0.01 | — |
| **ICON-EU (real model levels)** | **10.4** | **−0.4e-4** | 2.6e-4 | **0.37** | **−73 m** |

Model levels — *not* ICON-EU's pressure levels, which barely help — make the mixed-layer
gradient a genuine measurement (bias 5× smaller) and fix the ABL itself (−696 m → −73 m).

The production **hybrid** product (`PYROCONV_SOURCE=hybrid` in `tests/pyroconv_daily.py`)
takes the atmospheric profile from ICON-EU model levels (`ml_method="fit_in_ml"`, a real
least-squares dtheta/dz), regrids it onto the 2.2 km grid, and keeps ICON-2I's 2.2 km surface
fields for the fuel gate — where fine terrain actually matters. The threshold still sits inside
the ±2.6e-4 error bar, so class *counts* are indicative, not exact; but the diagnostic is now a
measurement rather than a proxy. `scripts/validation/sounding_degradation.py` and the ICON-EU
level lists reproduce the table above.

## 6. ERA5 referees the ICON-EU product against the Catalan fire service

Compared over Catalonia, the ICON-EU hybrid product and the Catalan fire-service (**Bombers**)
ICON-EU product agree overnight and (after relaxing the moisture gate, §7) on pyroCb, but
diverge sharply on the daytime overshoot/surface balance: the Bombers are surface-heavy
(~62-68% surface, 12-15% overshoot at midday) where the pyflam product is overshoot-heavy
(~46-52% overshoot). No same-day inland Catalan sounding or WRF was available to adjudicate, so
we used **ERA5** — an independent model (ECMWF reanalysis, 31 km, 22+ pressure levels), fetched
for the most recent available days (2026-07-06..10) and run through the *same* classifier.

**ERA5 backs the pyflam product, decisively.** Inland midday LCL/ABL ≈ 1.0-1.2 with a deep
~2600-3500 m ABL — the overshooting regime — and the class balance is **~58% overshoot / ~14%
surface**, i.e. *more* overshoot-heavy than the pyflam product and nothing like the Bombers'
surface-heavy picture. (Coastal cells show LCL/ABL ~2-2.8, the sea-breeze shallow-ABL case,
which is why a coastal Barcelona sounding alone was misleading.) ERA5 also produces the
resilient/deep tail on the unstable days (up to ~12% deep on 07-10), independently confirming
the `rh_top_moist = 60` relaxation. So three independent sources — inland radiosondes (Madrid),
ERA5, and the resolution model below — all place the truth with the pyflam product.

**Relationship model (reduced → full levels within ERA5).** Degrading ERA5 to each model's
level count, against the full 22-level reference:

| level set | ABL (m) | LCL/ABL | ML grad | bias vs full |
|:--|--:|--:|--:|:--|
| full (22) | 2585 | 1.16 | 1.15e-3 | — |
| ICON-EU-like (12) | 2759 | 1.15 | 1.05e-3 | **ABL +174 m, ratio −0.01** |
| ICON-2I-like (5) | 3209 | 1.18 | 6.6e-4 | ABL +623 m, gradient halved |

At ICON-EU's level density the diagnostics are **essentially unbiased** vs the full reanalysis,
so the "scale-and-adapt" correction to carry to ICON-EU is negligible — ICON-EU can be used
directly. The 5-level path carries a real +623 m ABL bias and a halved gradient, which is
exactly why the model-level hybrid was needed. This is the gridded ERA5 version of the
radiosonde degradation in §5, over many more columns.

```bash
pip install cdsapi   # + a ~/.cdsapirc key (https://cds.climate.copernicus.eu/how-to-api)
PYTHONPATH=../../src python era5_fetch_catalonia.py          # writes era5_{pl,sfc}_cat.nc
PYTHONPATH=../../src python era5_regime_relationship.py      # the tables above
```

The `.nc` inputs (~1 MB) are not committed — re-fetch them with the script (ERA5T lags ~5 days,
so pick past dates).

## References

- Castellnou, M., et al. (2022). *JGR: Atmospheres, 127*, e2022JD036920. https://doi.org/10.1029/2022JD036920
- Liu, S., & Liang, X.-Z. (2010). *Journal of Climate, 23*(21), 5790–5809. (source of the 1.1e-3 K/m threshold)
- Zhang, Y., et al. (2014). *GMD, 7*(6), 2599–2611. (Rib_c = 0.33, 200 m start)
- Vogelezang & Holtslag (1996). *Boundary-Layer Meteorology, 81*, 245–269.
- Holzworth, G. C. (1964). *Monthly Weather Review, 92*(5), 235–242. (parcel method)
- Vilà-Guerau de Arellano, J., et al. (2015). *Atmospheric Boundary Layer.* Cambridge UP.
- Durre, Vose & Wuertz (2006). *Journal of Climate, 19*(1), 53–68. (IGRA)
- Hersbach, H., et al. (2020). The ERA5 global reanalysis. *QJRMS, 146*(730), 1999–2049. https://doi.org/10.1002/qj.3803
