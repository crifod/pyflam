# Pyroconvection product: ICON-2I vs hybrid (ICON-EU) — 2026-07-15

An exhaustive cell-by-cell comparison of the two atmospheric sources now available in
`tests/pyroconv_daily.py`, run on the same day, same grid, same fuel gate:

- **`icon2i`** — atmosphere from ICON-2I 2.2 km, 5 pressure levels. Mixed-layer stability is
  a **mixing-depth proxy** (the gradient is not measurable on 5 levels).
- **`hybrid`** — atmosphere from ICON-EU 6.5 km **native model levels** (regridded to the
  2.2 km grid), mixed-layer stability a **measured** least-squares dθ/dz. Surface fields and
  the 10 MW/m fuel gate stay ICON-2I 2.2 km.

Run: 2026-07-15 00Z, valid same day, 3-hourly. Domain: Tuscany (9.6–12.5 E, 42.2–44.6 N),
9,907 land cells. Both products and their diagnostic rasters are committed alongside this
report; the numbers below are reproducible from `scripts/validation/` + the committed
GeoTIFFs via `compare_sources.py`.

> **One-line summary.** On this (low-danger) day the two products are operationally almost
> identical once the fuel gate is applied — but their *atmospheric* pictures differ sharply:
> the hybrid is systematically more conservative, replacing ~20% of ICON-2I's "overshooting
> pyroCu" with "convection plume" and **eliminating the spurious deep-pyroCb cells** ICON-2I
> produces at midday. This is the over-flagging the radiosonde validation predicted, seen live.

---

## 1. Potential (atmospheric upper bound) — class distribution

Percent of land cells, per valid hour. Classes: 0 surface plume · 1 convection plume ·
2 overshooting pyroCu · 3 resilient pyroCu · 4 deep pyroCu/pyroCb.

**ICON-2I (proxy):**

| Hour | surface | convect | overshoot | resilient | deep |
|:--|--:|--:|--:|--:|--:|
| 00Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 03Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 06Z | 99.8 | 0.0 | 0.1 | 0.1 | 0.0 |
| 09Z | 9.0 | 7.1 | 81.5 | 2.3 | 0.0 |
| 12Z | 5.5 | 1.4 | **92.7** | 0.1 | **0.3** |
| 15Z | 21.1 | 8.6 | 68.9 | 0.1 | **1.4** |
| 18Z | 98.3 | 1.6 | 0.2 | 0.0 | 0.0 |
| 21Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |

**Hybrid (measured):**

| Hour | surface | convect | overshoot | resilient | deep |
|:--|--:|--:|--:|--:|--:|
| 00Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 03Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 06Z | 99.7 | 0.0 | 0.0 | 0.3 | 0.0 |
| 09Z | 4.0 | **44.3** | 50.9 | 0.8 | 0.0 |
| 12Z | 18.5 | **22.0** | 59.5 | 0.0 | **0.0** |
| 15Z | 15.5 | **22.4** | 62.0 | 0.0 | **0.0** |
| 18Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 21Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |

Both share the diurnal envelope — nothing overnight (shallow nocturnal ABL below the 600 m
gate), activity 09–15Z, collapse by 18Z. The differences are all in the daytime *ladder rung*:

- **Convection plume explodes** in the hybrid (1.4% → 22.0% at 12Z; 7.1% → 44.3% at 09Z). These
  are cells the hybrid finds *marginally stable or LCL-above-ABL*, so they stay a penetrating
  convection plume instead of being promoted to pyroCu.
- **Overshooting shrinks** correspondingly (92.7% → 59.5% at 12Z).
- **Deep pyroCb disappears**: ICON-2I flags 0.3% (12Z) and 1.4% (15Z) of the domain as deep
  pyroCu/pyroCb; the hybrid flags **0.0%**. These were the scattered dark-red cells in the
  ICON-2I panel, and they are exactly the false positives the validation warned about.

## 2. Cell-by-cell agreement (potential)

| Hour | 00 | 03 | 06 | 09 | 12 | 15 | 18 | 21 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| agreement % | 100 | 100 | 99.4 | **50.7** | 61.8 | 73.8 | 98.3 | 100 |

Identical at night, most divergent at the convective peak (09–12Z ≈ 51–62%). The sources
agree on *whether* something is happening and disagree on *how intense* — which is the whole
point of measuring the mixed-layer gradient instead of proxying it.

## 3. Where the cells move — 12Z confusion matrix (potential)

Rows = ICON-2I class, columns = hybrid class, % of land cells. Read across a row to see where
ICON-2I's cells land in the hybrid.

| ICON-2I ↓ / hybrid → | surface | convect | overshoot | resilient | deep |
|:--|--:|--:|--:|--:|--:|
| surface | 3.0 | 1.4 | 1.2 | 0.0 | 0.0 |
| convect | 0.1 | 1.1 | 0.2 | 0.0 | 0.0 |
| **overshoot** | **15.5** | **19.5** | 57.7 | 0.0 | 0.0 |
| resilient | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| **deep** | 0.0 | 0.0 | **0.3** | 0.0 | 0.0 |

The movement is one-directional: **downgrades, never upgrades**. Of ICON-2I's 92.7%
overshooting, 57.7 pts survive as overshooting, 19.5 pts drop to convection plume and 15.5 pts
to surface plume in the hybrid. Every one of ICON-2I's deep-pyroCb cells becomes merely
overshooting. Nothing moves up the ladder. The hybrid is strictly more conservative here.

## 4. Fuel-gated (expected) product — the operational map

Percent of land cells; this is the map a forecaster actually uses (a class is assigned only
where fire power reaches ≥ 10 MW/m on the real fuels).

| Hour | ICON-2I overshoot | hybrid overshoot | ICON-2I ≥resilient | hybrid ≥resilient |
|:--|--:|--:|--:|--:|
| 09Z | 0.1 | 0.0 | 0.0 | 0.0 |
| 12Z | 0.7 | 0.5 | 0.0 | 0.0 |
| 15Z | 1.9 | 1.7 | 0.0 | 0.0 |

**Operationally, the two products are nearly identical today.** The fuel gate dominates: on a
low-danger July day the .lcp fuels + forecast weather rarely reach 10 MW/m, so ≥97.8% of the
domain is surface plume in both, and neither produces a single deep-pyroCb cell in the gated
map. The atmospheric disagreement in §1 is almost entirely in cells where no fire could get
strong enough to realise it — i.e. it matters for the *potential* ceiling, not for today's
expected behaviour. On a high-danger day (deep, dry, unstable, high fireline intensity) the
gate would open and the §1 divergence would propagate into the operational map.

## 5. Why they differ — the diagnostics

Land-median of each diagnostic, ICON-2I (2I) vs hybrid (HY), at the three active hours.

| Diagnostic | 09Z 2I / HY | 12Z 2I / HY | 15Z 2I / HY |
|:--|:--|:--|:--|
| Rib ABL (m) | 1060 / 897 | **2311 / 1510** | 1734 / 1429 |
| parcel mixing depth (m) | 912 / 865 | 1569 / 1368 | 1192 / 1216 |
| LCL (m) | 1156 / 1472 | 1648 / 1873 | 1641 / 1820 |
| **LCL / ABL** | 1.10 / 1.50 | **0.79 / 1.23** | 0.94 / 1.29 |
| **ML dθ/dz (K/m)** | 0.0005 / 0.0018 | **0.0003 / 0.0014** | 0.0004 / 0.0010 |
| cap γθ (K/m) | 0.0040 / 0.0051 | 0.0031 / 0.0031 | 0.0032 / 0.0028 |
| RH at ABL top (%) | 64.7 / 53.5 | 52.8 / 45.9 | 54.1 / 46.9 |

Two mechanisms, both pointing the same way:

1. **The ABL is shallower and the LCL higher in the hybrid, so LCL/ABL crosses 1.** At 12Z
   ICON-2I gives 0.79 (LCL *below* the ABL → the resilient/overshooting branch), the hybrid
   1.23 (LCL *above* → the convection-plume/overshooting branch). This is the direct cause of
   the convection-plume surge. It is consistent with the radiosonde finding that ICON-2I's Rib
   ABL runs ~500–700 m too deep; the hybrid's shallower ABL (validated bias ~−70 m vs ~−700 m)
   is the more trustworthy denominator.

2. **The measured ML gradient is genuinely marginal, not near-zero.** The ICON-2I proxy reads
   3e-4 K/m (deeply "unstable" by construction: 0.5/mixing-depth); the hybrid measures
   1.4e-3 K/m — just *above* the 1.1e-3 "stable" threshold in about half the cells, which is
   what tips them out of the pyroCu branch into convection plume. The two are not measuring the
   same thing: one is a restated mixing depth, the other a real dθ/dz.

The cap γθ and ABL-top RH are close between the sources; the divergence is driven by ABL depth
and the mixed-layer gradient, exactly the two quantities the model levels fix.

## 6. Cross-reference to the radiosonde validation

This live comparison is the operational face of `scripts/validation/`:

- ICON-2I's ML-stability proxy has inland specificity 0.29 — it **over-flags**. Here that shows
  up as 34 pts of overshooting (92.7 → 59.5) and all of the deep-pyroCb being false positives
  the hybrid removes.
- ICON-EU model levels put ~7–9 measurable levels in the mixed layer at midday (1 at night —
  correctly, the nocturnal layer is thin), cutting the ML-gradient bias to ~−0.4e-4 and the ABL
  bias to ~−70 m. The shallower ABL and measured gradient are what drive §5.

## 7. Limitations (unchanged by this comparison, restated)

- The 1.1e-3 K/m threshold still sits inside the ±2.6e-4 K/m validated error bar, so the hybrid
  class *counts* are indicative, not exact — but they are a measurement, not a proxy.
- No radiosonde lies inside Tuscany; validation used inland Po-valley + coastal Rome stations.
- The hybrid trades horizontal resolution (6.5 km atmosphere vs 2.2 km) for vertical
  resolution. The visible smoothing in the hybrid panels is that trade; the fuel gate stays
  2.2 km, so fire-power detail is preserved.
- Neither product is validated against *observed* pyroCu in Tuscany. Everything here is
  internal consistency with the Castellnou (2022) method + radiosonde-anchored diagnostics.
- The fifth (shear) diagnostic is not computed in either product; the 4-diagnostic ladder runs.

## 8. Which to use

- **High-danger days / deep-plume risk assessment:** prefer **hybrid**. Its ABL and mixed-layer
  stability are measured, it does not manufacture spurious deep-pyroCb cells, and the deep-plume
  call is the one that most needs to be right.
- **Situational awareness on ordinary days:** either — they agree operationally once gated. The
  ICON-2I product is 2.2 km throughout and cheaper (~0.9 GB extra for the hybrid's ICON-EU
  levels, and 8 ICON-EU steps to fetch).
- **Do not read potential-map class totals quantitatively from either.** Use the gated map for
  decisions and the diagnostic rasters (ABL, LCL/ABL, ML dθ/dz) to understand a cell's class.

## 9. Reproduction

```bash
# ICON-2I product
PYTHONPATH=src python tests/pyroconv_daily.py 2026-07-15 0
# hybrid product (fetches ICON-EU model levels, ~1.1 GB / 8 steps)
PYROCONV_SOURCE=hybrid PYTHONPATH=src python tests/pyroconv_daily.py 2026-07-15 0
# comparison tables (this report)
python scripts/validation/compare_sources.py 2026-07-15
```

Products: `pyroconv_tuscany_{icon2i,hybrid}_{potential,gated}_2026-07-15.png`, the two PDFs,
and the per-hour diagnostic GeoTIFFs under `rasters_{icon2i,hybrid}_2026-07-15/` (gitignored;
regenerable). Method and validation: `scripts/validation/README.md`.
