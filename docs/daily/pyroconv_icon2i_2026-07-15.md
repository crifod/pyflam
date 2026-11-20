---
title: "Tuscany Pyroconvection-Type Forecast -- ICON-2I 2.2 km -- VALID 2026-07-15"
subtitle: "3-hourly, 24 h. Run 2026-07-15 00Z. Method after Castellnou et al. (2022), JGR-Atmos."
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
per-cell heights from the model geopotential (FI), the mixing depth from the **bulk Richardson number** (first Rib >= 0.33
above 200 m AGL, referenced to the 2 m / 10 m state), the LCL from the exact
**Bolton (1980)** formula, and the humidity at the ABL top from the model's RH. The cap
gamma-theta is taken over ABL+200 m to ABL+1200 m.

The **mixed-layer stability** diagnostic is the weakest link in this product, and is a
*proxy*, not a measurement. Validated against IGRA radiosondes (JJA 12Z, period of record,
n=2835), a 5-pressure-level column contains only 0-2 model levels inside the mixed layer --
one or none in 92% of coastal columns -- and the lowest sits in the superadiabatic surface
layer. The mixed-layer dtheta/dz is therefore **not measurable from this archive**: fitting
it across the in-mixed-layer levels has no skill (Youden J ~ 0.00), and measuring theta up to
the Rib ABL top folds in the entrainment jump (Delta-theta), which Castellnou et al. (2022,
sec.2.1.1) treat as a variable *separate* from the gradient the ladder conditions on.

What is used instead is the **parcel mixing depth** as a proxy: at the 1.1e-3 K/m threshold
the criterion reduces to "well-mixed layer >= 455 m deep". It is the only candidate with
skill (J = 0.29 inland / 0.55 coastal, r = +0.50 against the radiosonde truth). It
**over-flags**: inland specificity is 0.29, so of the columns that are truly stable it still
calls ~71% pyroCu-capable (sensitivity 0.99 -- it rarely misses a capable column).

**Consequence: read the classes as a screening flag, not as calibrated counts.** The class
*totals* on these maps are not quantitatively trustworthy. The hybrid product
(`PYROCONV_SOURCE=hybrid`, ICON-EU model levels) measures this gradient properly; on this
5-level source, set `PYROCONV_ML_METHOD=surface_to_abl` or `mid_layer` for the superseded
measurements.

**Ladder actually used for this run: `noshear`.** The full method has a fifth
diagnostic -- the distance from the ABL/LCL to the height of maximum wind shear -- which
this product does not yet compute, so the four-diagnostic ladder runs. The shear clause is a
*necessary* condition for the top class, so omitting it makes class 4 somewhat **easier** to
reach here than in the full method. Treat class 4 as an alert to inspect the column, not as a
calibrated probability.

The classes express atmospheric predisposition **given a fire of sufficient power**;
in the gated panel that power is computed, not assumed. They are paper-informed
thresholds, not locally validated ones.

## Expected pyroconvection type -- FUEL-GATED

![gated](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_icon2i_gated_2026-07-15.png){width=100%}

## Atmospheric potential -- UPPER BOUND (assumes a pyroCu-capable fire in every cell)

![potential](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_icon2i_potential_2026-07-15.png){width=100%}

## Class scale (low -> high pyroconvective activity)

| Level | Colour | Class | Meaning |
|:--:|:--|:--|:--|
| 0 | white | Surface plume | Buoyant smoke plume; no significant cloud development. |
| 1 | green | Convection plume | Plume penetrates a stable mixed layer; condensation possible, no pyroCu. |
| 2 | yellow | Overshooting pyroCu | Brief pyrocumulus; cloud base above the mixing height (LCL/ABL > 1). |
| 3 | orange | Resilient pyroCu | Persistent pyrocumulus in an unstable column (LCL/ABL < 1). |
| 4 | dark red | Deep pyroCu / pyroCb | Deep pyroconvection / pyrocumulonimbus; weak upper cap lets the plume deepen. |

## Classification thresholds (ladder in use: `noshear`)

| Diagnostic | Threshold | Effect |
|:--|:--|:--|
| ML stability *proxy*: parcel mixing depth (see Method) | depth < 455 m ("stable") | Convection plume only -- no pyroCu |
| ML stability proxy | depth >= 455 m (= dtheta/dz <= 1.1e-3 K/m) | Column is pyroCu-capable (weak filter: spec. 0.29 inland) |
| LCL / ABL ratio | 1.0 -- 1.60 | Overshooting pyroCu (brief) |
| LCL / ABL ratio | < 1.0 | Resilient pyroCu (persistent) |
| LCL / ABL ratio | <= 1.10 (+ conditions below) | Admissible for deep pyroCu / pyroCb |
| Cap gamma-theta (ABL+200 m -> ABL+1200 m) | <= 4.2e-3 K/m (weak cap) | Permits deepening to pyroCb |
| Cap gamma-theta | >= 4.8e-3 K/m (strong cap) | Inhibits deepening (resilient at most) |
| RH at the ABL top (mean, ABL +/- 150 m) | >= 60% | Required for classes 3 and 4 (was 80%; relaxed for dry fire weather) |
| Shear-maximum distance / ABL | <= 0.30 | Required for class 4 **in the 5-diagnostic ladder only** (not resolvable here -- see Method) |
| Fireline intensity (fuel gate, gated panel) | >= 10 MW/m | Minimum fire power for any pyroCu (Tedim et al. 2018) |
| ABL depth | < 600 m | Held at surface plume (mixing too shallow) |
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

Forcing: ICON-2I 2.2 km full-Italy GRIB (MISTRAL / AgenziaItaliaMeteo, CC-BY), Tuscany subset: geopotential, temperature, RH and wind on 5 pressure levels plus the 2 m / 10 m state, surface pressure and orography. Classifier: pyflam.pyroconvection_type (bulk-Richardson ABL, Bolton LCL, mixed-layer dtheta/dz proxy, cap gamma-theta, ABL-top RH; no surface CAPE).
Gate: Rothermel + Cruz-2005 crown on the .lcp fuels with forecast moisture/wind.
Per-hour diagnostic rasters (ABL, LCL, LCL/ABL, ML dtheta/dz, cap, RH-top) accompany the classes.
Province borders: ISTAT-derived (openpolis geojson-italy). Generated by tests/pyroconv_daily.py.
