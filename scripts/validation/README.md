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

Fixing this properly needs more vertical levels — ERA5's 37 pressure levels, or ICON native
model levels. The open-data archive publishes neither.

## References

- Castellnou, M., et al. (2022). *JGR: Atmospheres, 127*, e2022JD036920. https://doi.org/10.1029/2022JD036920
- Liu, S., & Liang, X.-Z. (2010). *Journal of Climate, 23*(21), 5790–5809. (source of the 1.1e-3 K/m threshold)
- Zhang, Y., et al. (2014). *GMD, 7*(6), 2599–2611. (Rib_c = 0.33, 200 m start)
- Vogelezang & Holtslag (1996). *Boundary-Layer Meteorology, 81*, 245–269.
- Holzworth, G. C. (1964). *Monthly Weather Review, 92*(5), 235–242. (parcel method)
- Vilà-Guerau de Arellano, J., et al. (2015). *Atmospheric Boundary Layer.* Cambridge UP.
- Durre, Vose & Wuertz (2006). *Journal of Climate, 19*(1), 53–68. (IGRA)
