---
title: "Tuscany Pyroconvection-Type Forecast -- ICON-EU model levels + ICON-2I 2.2 km gate -- VALID 2026-08-03"
subtitle: "3-hourly, 24 h. Run 2026-08-03 00Z. Method after Castellnou et al. (2022), JGR-Atmos."
geometry: a4paper, landscape, margin=1.2cm
fontsize: 9pt
lang: en
header-includes: |
  \usepackage{float}
  \floatplacement{figure}{H}
---

## How to read this product

Each map answers one question. Read them in order -- but see the warning below question 3:
they are **not** nested, and a cell can pass question 3 while failing question 2.

| # | Question | Map | The test |
|:--|:--|:--|:--|
| 1 | What is the most intense pyroconvection this **column** could sustain? | POTENTIAL | the ladder, with no fire-side condition at all |
| 2 | Is there enough fire here to make a **plume** at all? | FUEL-GATED | Byram intensity >= 10 MW/m on the .lcp fuels |
| 3 | Is there enough fire here to make a **pyroCb** in this specific column? | PFT MARGIN | total firepower >= that column's own PFT |

Question 1 is a property of the atmosphere alone -- no fuel or land information enters it,
not even in deciding where it is defined. Questions 2 and 3 both condition on fire, and both
descend from the same Byram intensity field, but they are **not** the same test and they do
not agree: (2) compares power *per metre of front* against a fixed 10 MW/m constant, while
(3) converts that to a *total* power through an assumed head-fire length and compares it
against a threshold computed for that column. Cells routinely clear (2) and fail (3) by
orders of magnitude. Question 3 is the operationally decisive one and the hardest to satisfy.

## Method and its limits (read this before using the classes)

The column is classified from a **real vertical profile**, not a standard atmosphere:
per-cell heights from the ICON-EU model-level heights (HHL), the mixing depth from the **bulk Richardson number** (first Rib >= 0.33
above 200 m AGL, referenced to the 2 m / 10 m state), the LCL from the exact
**Bolton (1980)** formula, and the humidity at the ABL top from the model's RH. The cap
gamma-theta is taken over ABL+200 m to ABL+1200 m.

The **mixed-layer stability** here is a genuine measurement, not a proxy. The ICON-EU
native model levels put ~10 levels inside the mixed layer (this run: 15 median over land at
12-15Z, when the layer is mature -- the night and evening counts are far lower, but nothing is
classified then). Over the classified land the fit is supported in **99%** of columns at peak of day (i.e. that share holds at least 3 levels inside the layer; the rest return a nan gradient and leave the class map, which is also what thins the evening panels). So the
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

## 1. Potential -- what could this column sustain? (atmospheric upper bound)

![Question 1 -- POTENTIAL. Pyroconvection class from the atmosphere alone, no fire-side condition. Upper bound, not an expectation.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_potential_2026-08-03.png){width=100%}

The ladder with the fire-side condition switched off entirely, so this is the atmosphere's
answer and nothing else: the most intense pyroconvection the column would support **given a
fire able to exploit it**. That fire is assumed, not computed -- it enters through the ladder's
calibrated thresholds, not through any fuel map. Use this map to see where the atmosphere is
the binding constraint, never as an expectation. Note that the classifier runs on any usable
column, including over water; the sea is masked in the figure but the underlying
`pyroconv_potential_*.tif` rasters are not, so zonal statistics taken straight off them must
be restricted to land.

## 2. Fuel-gated -- is there enough fire for a plume at all?

![Question 2 -- FUEL-GATED. The same ladder with every cell below 10 MW/m of Byram intensity on the Tuscany .lcp fuels forced to class 0. The expected product.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_gated_2026-08-03.png){width=100%}

The same ladder, on the same columns, with each cell forced to class 0 where the Tuscany .lcp
fuels and the forecast surface weather do not support **10 MW/m** of Byram fireline intensity
(Tedim et al. 2018). The gate is binary and one-directional: above the threshold the fuels
have no further influence, below it the cell reads surface plume whatever the atmosphere says.
This map therefore cannot exceed the potential map anywhere -- every difference between the
two is a cell knocked to 0 -- and it is the expected, operational product.

## 3. PFT margin -- is there enough fire for a pyroCb in *this* column?

![Question 3 -- PFT MARGIN. Firepower over each column's own PyroCb Firepower Threshold, log scale pivoting at the criterion. Rings meet the criterion and pass the fuel gate; crosses meet it only because the PFT has degenerated.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_pft_margin_2026-08-03.png){width=100%}

Per-cell **total** firepower divided by that column's own **PyroCb Firepower Threshold**
(Tory & Kepert 2021, eq. 31: `PFT = 0.3 z_fc^2 U_ML dtheta_fc`, in GW). At or above 1 the fire
is powerful enough to force a pyroCb *in that column*; below 1 it is not, however favourable
the atmosphere. The colour pivots exactly at 1 and the scale is logarithmic, because the field
spans four decades.

Byram intensity is power per metre of front and the PFT is a total power, so bridging them
needs an active head-fire length. Rather than assume one it is taken from the fire behaviour
the Wildfire Data Portal publishes -- the median over its 20 fires with both burn ratio and
ROS, 700 m -- with a convective fraction of 0.7 of the released heat entering the plume
(Tory & Kepert app. D). Both are assumptions, and the margin scales linearly with each: a
5 km head fire, which Tory & Kepert give for an extreme case, would multiply every margin here
by about seven. **Read the margin as an order of magnitude, not a number.**

This is why the map matters despite that caveat: over most of the burnable domain the margin
is not near 1 but two to three decades below it, so the conclusion "no pyroCb here" is robust
to the head-fire assumption in a way that a marginal cell would not be. Where the margin does
approach 1, the assumption becomes load-bearing and the cell deserves judgement rather than
the map.

### Questions 2 and 3 are not nested -- read the markers, not the colour

A cell can clear the PFT criterion while **failing** the 10 MW/m fuel gate, and on this run
most of them do. The mechanism is the denominator: the PFT is
`0.3 z_fc^2 U_ML dtheta_fc`, so in a shallow, weakly-capped, light-wind column it collapses --
against a domain median near 90 GW, the cells at margin >= 1 have PFTs of a few GW. A fire of
2-8 MW/m then "passes" a threshold that a fire of that size cannot physically be said to have
beaten, because at that intensity there is no pyroCu to deepen in the first place.

The map marks the two cases apart, and only the first is a forecast:

* **black ring** -- margin >= 1 *and* the fuel gate passed. A genuine pyroCb candidate.
* **grey cross** -- margin >= 1 but the gate failed. The PFT has degenerated, not been beaten.
  Read these as an artefact of the threshold's form at small `z_fc`, not as a signal.

The **also gate-passing** column of the table below is therefore the operational count. Taking
the raw margin count instead would overstate the candidate area several-fold.

| Hour | % grid with firepower | % margin >= 1 | % margin >= 0.1 | max margin | cells >= 1 | **also gate-passing** |
|:--|--:|--:|--:|--:|--:|--:|
| 00Z | 0.0 | 0.0 | 0.0 | nan | 0 | **0** |
| 03Z | 0.0 | 0.0 | 0.0 | nan | 0 | **0** |
| 06Z | 0.0 | 0.0 | 0.0 | 0.02 | 0 | **0** |
| 09Z | 25.4 | 0.1 | 2.3 | 1.49 | 2 | **0** |
| 12Z | 25.4 | 0.8 | 7.9 | 27.46 | 30 | **16** |
| 15Z | 25.6 | 0.3 | 6.1 | 1.58 | 9 | **7** |
| 18Z | 6.4 | 0.0 | 0.5 | 0.25 | 0 | **0** |
| 21Z | 0.1 | 0.0 | 0.0 | 0.01 | 0 | **0** |

*% grid with firepower* is the share of the domain that has both a classifiable column and
burnable fuel; the two margin columns are percentages **of that share**, not of Tuscany. A
percentage of the whole domain would be dominated by cells that carry no fuel and can never
contribute, and would move for reasons unrelated to the forecast. *cells >= 1* is the raw
count behind the first percentage, carried because at these rates a percentage rounds to
something that looks like nothing when it is in fact a handful of specific places.
***also gate-passing*** applies the joint criterion of the previous section and is the count
to act on; the gap between it and *cells >= 1* is the degenerate-PFT population.

## 4. Dry-pyrocloud decoupling -- DIAGNOSTIC (no class label)

![Supporting diagnostic -- DRY-PYROCLOUD DECOUPLING. fireABL/ABL for a reference intense fire. Continuous, no class label.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_decoupling_2026-08-03.png){width=100%}

The **decoupling ratio** fireABL / ABL is how high a reference intense fire
(a 10 K plume temperature excess -- the intense end of the GRAF in-plume measurements,
0.1-13.1 K) would grow its own boundary layer by *sensible heat alone*, divided by the ambient
ABL. It is the **dry** counterpart to the moist class maps above (Castellnou et al. 2022;
Castellnou Ribau et al. 2024): values well above 1 mark deep, hot, dry columns where a fire can
punch through and decouple from the surface *even where the moist ladder scores low*, which is
exactly the situation the gated map under-reports. Like the potential map it assumes a fire
everywhere -- here a fixed reference flux rather than a computed one -- so it is an upper bound,
not an expectation, and no dry/moist LCL split is applied (the +1 km literature offset is not
supported by the GRAF prototype labels). fireABL from
`pyflam.atmosphere.fire_induced_abl_grid`, the parcel intersected with the real theta(z) stack;
cells the classifier rejects (no usable column, or ABL below 150 m) are left blank rather
than painted with a bare quotient.

## 5. Predicted plume-top height -- the only field validated against observation

![PREDICTED PLUME TOP. The cost ladder solved for height instead of firepower, for a declared 5 GW fire. Continuous, no class label.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_plumetop_2026-08-03.png){width=100%}

The three questions above ask *how much fire* a given height costs. Solving the same
inequality the other way -- for the largest height a given fire can afford -- answers the
question a forecaster actually has, and returns a continuous field instead of a class.

It is also the only product here that has been **checked against observation end to end**.
Comparing a plume top against a rung's own target height puts that target on both sides of the
test and measures mostly its own circularity; solving for the height puts a predicted number
against a measured one and nothing else. Scored that way against 1231 plumes digitised from
MISR stereo imagery over 7 regions and 11 biomes, with nothing fitted: Spearman **+0.46**
against the observed top, **+125 m** median bias, 576 m mean absolute error, and 12.8 % better
than predicting a constant on a scale-free error -- ahead of firepower alone (+0.37) and of
boundary-layer depth alone (+0.34), and ordering the regional medians at rho +0.61.

**Read the firepower as a scenario.** The map is drawn for a declared 5 GW *total*
power, stated here because a total power cannot be built from a Byram intensity without
assuming a head-fire length, and that assumption is exactly what this line of work exists to
avoid. Doubling the reference fire does not double the height -- the cost grows with the cube
of the climb -- but the field does shift, and it must not be read as a per-cell estimate of
what will actually burn.

**Two limits worth stating.** The validation sample is MISR's, so it is fixed at roughly 10:30
local solar time: the boundary layer is still growing and fires are smaller than they will be
in the afternoon this product forecasts. And the field is produced only where the vertical
grid resolves the layer the cap is read across, which on this run means the model-level source;
a pressure-level fallback leaves it blank rather than interpolating one.

Unlike the class rasters, `diag_plume_top_*.tif` and `diag_form_*.tif` are written **already
masked** to the classifier's `valid` field, so a zonal statistic taken straight off them cannot
include columns the classifier rejected. The other `diag_*.tif` stay raw, as they always have.

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

Forcing: atmosphere from ICON-EU 6.5 km native model levels (DWD open data, CC-BY), lowest ~24 levels; ~15 inside the mixed layer here (land, 12-15Z), regridded to the ICON-2I 2.2 km grid. Surface fields and fuel gate from ICON-2I 2.2 km (MISTRAL / AgenziaItaliaMeteo). Classifier: pyflam.pyroconvection_type (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta, ABL-top RH; no surface CAPE).
Gate: Rothermel + Cruz-2005 crown on the .lcp fuels with forecast moisture/wind.
Per-hour diagnostic rasters accompany the classes: ABL, parcel_ml, residual_ml, LCL, LCL/ABL,
ML dtheta/dz, cap gamma-theta, RH-top, fireABL, decoupling, and (hybrid only) delta_theta,
firecape, penetration, pft_gw (PyroCb Firepower Threshold, GW), z_fc, delta_theta_fc, u_ml,
firepower_gw, pft_margin.
Province borders: ISTAT-derived (openpolis geojson-italy).
Generated by tests/pyroconv_daily.py at pyflam `894a4b9+local-changes`, 2026-08-03 10:21 UTC -- frozen with that
commit's physics, since DWD drops the run after ~24 h.
