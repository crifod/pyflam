---
title: "Tuscany Pyroconvection-Type Forecast -- ICON-EU model levels + ICON-2I 2.2 km gate -- VALID 2026-07-28"
subtitle: "3-hourly, 24 h. Run 2026-07-26 00Z. Method after Castellnou et al. (2022), JGR-Atmos."
geometry: a4paper, landscape, margin=1.2cm
fontsize: 9pt
---

## How to read this product

Two panels are produced. The **fuel-gated** map is the expected, operationally
comparable product (the equivalent of the Catalan "tipus de piroconveccio" map):
a pyroCu/pyroCb class is assigned **only where a fire could actually reach >= 10
MW/m** of fireline intensity on the real Tuscany fuels. The **potential** map is an
unconditional **upper bound** -- it assumes a pyroCu-capable fire in *every* cell.
Use the gated panel for situational awareness; use the potential panel only to see
the atmospheric ceiling.

## Method and its limits (read this before using the classes)

The column is classified from a **real vertical profile**, not a standard atmosphere:
per-cell heights from the ICON-EU model-level heights (HHL), the mixing depth from the **bulk Richardson number** (first Rib >= 0.33
above 200 m AGL, referenced to the 2 m / 10 m state), the LCL from the exact
**Bolton (1980)** formula, and the humidity at the ABL top from the model's RH. The cap
gamma-theta is taken over ABL+200 m to ABL+1200 m.

The **mixed-layer stability** here is a genuine measurement, not a proxy. The ICON-EU
native model levels put ~10 levels inside the mixed layer (this run: 13 median over land at
12-15Z, when the layer is mature -- the night and evening counts are far lower, but nothing is
classified then). Over the classified land the fit is supported in **100%** of columns at peak of day (i.e. that share holds at least 3 levels inside the layer; the rest return a nan gradient and leave the class map, which is also what thins the evening panels). So the
mixed-layer dtheta/dz is a real least-squares fit across those levels. Validated against
IGRA radiosondes (JJA 12Z, period of record), the model-level fit cuts the gradient bias to
~-0.4e-4 K/m (from ~-4.7e-4 on ICON-2I's 5 pressure levels) and the Rib ABL bias to ~-70 m
(from ~-700 m). This is why the hybrid product exists: the ICON-2I 2.2 km open data cannot
resolve the mixed layer, and ICON-EU can.

The trade is horizontal resolution: the atmosphere is 6.5 km (regridded to the 2.2 km grid),
while the **fuel gate keeps ICON-2I's 2.2 km surface fields**, where fine terrain matters.
The mixed-layer gradient is now measured, but the 1.1e-3 K/m threshold still sits inside the
validated error bar (+/- ~2.6e-4), so treat class *counts* as indicative, not exact.

After sunset the layer the gradient is fitted over is the **residual layer** -- the near-neutral
air the decaying convective layer leaves behind -- rather than the shallow nocturnal stable
layer a surface-referenced parcel would find. By day the two coincide. This is what lets the
evening hours carry a forecast at all; it is a newer diagnostic than the rest of the ladder, so
weight those hours accordingly.

**Ladder actually used for this run: `shear`.** The full 5-diagnostic method ran: the model levels resolved a shear-maximum height, so class 4 additionally required that maximum to sit within 0.30 ABL of the ABL/LCL -- a *necessary* condition the reduced ladders omit. Class 4 here therefore carries its shear clause; it remains an alert to inspect the column rather than a calibrated probability.

The classes express atmospheric predisposition **given a fire of sufficient power**;
in the gated panel that power is computed, not assumed. They are paper-informed
thresholds, not locally validated ones.

## Expected pyroconvection type -- FUEL-GATED

![gated](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/forecast_2026-07-26_3day/pyroconv_tuscany_hybrid_gated_2026-07-28.png){width=100%}

## Atmospheric potential -- UPPER BOUND (assumes a pyroCu-capable fire in every cell)

![potential](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/forecast_2026-07-26_3day/pyroconv_tuscany_hybrid_potential_2026-07-28.png){width=100%}

## Dry-pyrocloud decoupling -- DIAGNOSTIC (no class label)

![decoupling](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/forecast_2026-07-26_3day/pyroconv_tuscany_hybrid_decoupling_2026-07-28.png){{width=100%}}

The **decoupling ratio** fireABL / ABL is how high a reference intense fire (a 10 K plume temperature excess, the intense end of the GRAF in-plume
measurements) would grow its own boundary layer by *sensible heat alone*, divided by the ambient ABL. It is the **dry** counterpart to the moist class map above (Castellnou et al. 2022; Castellnou Ribau et al. 2024): values well above 1 mark deep, hot, dry columns where a fire can punch through and decouple from the surface *even where the moist ladder scores low*. It is a diagnostic, not a calibrated class -- no dry/moist LCL split is applied (the +1 km literature offset is not supported by the GRAF prototype labels; a fit prefers ~ -0.5 km). fireABL from `pyflam.atmosphere.fire_induced_abl_grid` (parcel intersected with the real theta(z) stack).

## Class scale (low -> high pyroconvective activity)

| Level | Colour | Class | Meaning |
|:--:|:--|:--|:--|
| 0 | white | Surface plume | Buoyant smoke plume; no significant cloud development. |
| 1 | green | Convection plume | Plume penetrates a stable mixed layer; condensation possible, no pyroCu. |
| 2 | yellow | Overshooting pyroCu | Brief pyrocumulus; cloud base above the mixing height (LCL/ABL > 1). |
| 3 | orange | Resilient pyroCu | Persistent pyrocumulus in an unstable column (LCL/ABL < 1). |
| 4 | dark red | Deep pyroCu / pyroCb | Deep pyroconvection / pyrocumulonimbus; weak upper cap lets the plume deepen. |

## Classification thresholds (ladder in use: `shear`)

| Diagnostic | Threshold | Effect |
|:--|:--|:--|
| ML dtheta/dz (least-squares fit, in mixed layer) | > 1.1e-3 K/m (stable) | Convection plume only -- no pyroCu |
| ML dtheta/dz (measured on model levels) | <= 1.1e-3 K/m | Column is pyroCu-capable (bias ~-0.4e-4 K/m vs radiosondes) |
| LCL / ABL ratio | 1.0 -- 1.60 | Overshooting pyroCu (brief) |
| LCL / ABL ratio | < 1.0 | Resilient pyroCu (persistent) |
| LCL / ABL ratio | <= 1.10 (+ conditions below) | Admissible for deep pyroCu / pyroCb |
| Cap gamma-theta (ABL+200 m -> ABL+1200 m) | <= 4.2e-3 K/m (weak cap) | Permits deepening to pyroCb |
| Cap gamma-theta | >= 4.8e-3 K/m (strong cap) | Inhibits deepening (resilient at most) |
| RH at the ABL top (mean, ABL +/- 150 m) | >= 60% | Required for classes 3 and 4 (was 80%; relaxed for dry fire weather) |
| Shear-maximum distance / ABL | <= 0.30 | Required for class 4 -- **resolved and applied throughout this run** |
| Fireline intensity (fuel gate, gated panel) | >= 10 MW/m | Minimum fire power for any pyroCu (Tedim et al. 2018) |
| ABL depth | < 150 m | Not classifiable (implausible depth). No 600 m gate: it had no basis in Castellnou et al. (2022) and discarded 5 of the 8 campaign fires -- removed 2026-07-26 |
| Usable pressure levels | < 4 | Cell not classified |

## Reference cases (Castellnou et al. 2022, Table 1)

The paper's labelled events, which anchor the published 3-diagnostic ladder
(`pyflam.pyroconvection_type(ladder="castellnou")`, still the library default and
regression-tested against these rows). The profile ladders used for this map are
stricter at the top: they additionally require a moist ABL top and an LCL close to it.

| Case | Observed type | LCL/ABL | ML dtheta/dz | gamma-theta (700-500) |
|:--|:--|:--:|:--|:--:|
| T21 | Convection plume | -- | stable | -- |
| SCQ32 | Overshooting pyroCu | > 1 | neutral/unstable | -- |
| M11 | Resilient pyroCu | < 1 | unstable | 4.2e-3 (resilient cap) |
| SCQ41 | pyroCu, not pyroCb | < 1 | unstable | 5.1e-3 (strong cap) |
| SCQ51 | Deep pyroCu / pyroCb | < 1 | unstable | 3.9e-3 (weak cap) |

Forcing: atmosphere from ICON-EU 6.5 km native model levels (DWD open data, CC-BY), lowest ~24 levels; ~13 inside the mixed layer here (land, 12-15Z), regridded to the ICON-2I 2.2 km grid. Surface fields and fuel gate from ICON-2I 2.2 km (MISTRAL / AgenziaItaliaMeteo). Classifier: pyflam.pyroconvection_type (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta, ABL-top RH; no surface CAPE).
Gate: Rothermel + Cruz-2005 crown on the .lcp fuels with forecast moisture/wind.
Per-hour diagnostic rasters accompany the classes: ABL, parcel_ml, residual_ml, LCL, LCL/ABL,
ML dtheta/dz, cap gamma-theta, RH-top, fireABL, decoupling, and (hybrid only) delta_theta,
firecape, penetration.
Province borders: ISTAT-derived (openpolis geojson-italy). Generated by tests/pyroconv_daily.py.
