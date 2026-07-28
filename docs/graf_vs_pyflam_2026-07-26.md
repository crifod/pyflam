# GRAF vs pyflam hybrid — pyroconvection type, run 2026-07-26 00Z

Comparison of the Catalan fire service (GRAF / Bombers) pyroconvection-type product against
the pyflam hybrid product over central Italy, plus the code change the comparison motivated.

**Status: for evaluation.** Start at **§12**, the standing synthesis of every difference
found; §§1-11 are the working record. Written 2026-07-26; extended the same day with §§7-12 after a
literature and sonde investigation. Two corrections are applied (§5 nodata separation, §9
removal of the 600 m ABL gate). The **extent** divergence (§3.1) is now *explained* but not
resolved; the **evening** divergence (§3.2) is half resolved, with the remaining half
attributed to a distinct cause (§9.3).

---

## 1. What is being compared

| | GRAF (Bombers) | pyflam hybrid |
|---|---|---|
| Product | "Tipus de piroconvecció" | Pyroconvection type — POTENTIAL panel |
| Source | ICON-EU | ICON-EU model levels (+ ICON-2I 2.2 km for the grid and the fuel gate) |
| Run | 26 Jul 2026 00Z | 26 Jul 2026 00Z |
| Valid | 26, 27, 28 Jul — 3-hourly, 00–21Z | identical |
| Domain | Tuscany + Emilia-Romagna, Marche, Umbria, Lazio | Tuscany bbox 9.6–12.5 E, 42.2–44.6 N |
| Classification grid | native ~6.5 km | diagnostics regridded to 2.2 km, thresholded there |
| Fuel gating | none apparent | none in the POTENTIAL panel (the GATED panel is a separate product) |

The POTENTIAL panel is the like-for-like counterpart: it is 100 % ICON-EU physics —
`iconeu_diagnostics` uses ICON-EU's own surface state, and ICON-2I contributes only the target
grid. The GATED panel has no GRAF counterpart and is excluded from this comparison.

### 1.1 Legend inversion — read before comparing anything

The two products number their classes in **opposite** order. Colours and class names
correspond one-to-one; the integers do not. Comparing "class 4" to "class 4" inverts the
entire result.

| Colour | GRAF | pyflam |
|---|---|---|
| dark red | **1** — PyroCb | **4** — Deep pyroCu / pyroCb |
| orange | **2** — Resilient PyroCu | **3** — Resilient pyroCu |
| yellow | **3** — Overshooting PyroCu | **2** — Overshooting pyroCu |
| green | **4** — Convective plumes | **1** — Convection plume |
| white / grey | not in legend | **0** — Surface plume / **n/c** — not classifiable |

All statements below are by **colour and class name**.

---

## 2. Where the two agree

Agreement is strongest on day 0 and is substantive rather than incidental:

- **Night (00–06Z, all three days).** Both blank. pyflam: 99–100 % of land not classifiable
  (ABL below the 600 m floor). Independent of any threshold choice.
- **Day-0 peak hour.** Both peak at **15Z**.
- **Day-0 top class.** Both reach the highest class (dark red / pyroCb).
- **Day-0 feature placement.** Both put the strongest signal on a band along the Apennine
  crest continuing east-southeast toward Marche / San Marino. Most of GRAF's dark red on this
  band lies **outside** the pyflam domain, which is why the two look differently centred: the
  same feature is cropped at our NE boundary.

This is the part of the comparison that supports the method: same hour, same ceiling, same
structure, from an independent implementation of the same published ladder.

---

## 3. Divergences

### 3.1 Extent — pyflam is far more expansive at midday

pyflam POTENTIAL, % of land at class ≥ 2 (overshooting pyroCu or higher):

| Day | 09Z | 12Z | 15Z | 18Z |
|:--|--:|--:|--:|--:|
| 26/07 | 12.7 | 20.3 | **59.2** | 24.9 |
| 27/07 | 34.2 | 65.1 | **89.8** | 0.1 |
| 28/07 | 47.5 | 65.9 | **83.5** | 0.0 |

GRAF over the same hours on 27–28/07 shows only small isolated yellow patches. pyflam places
65–90 % of Tuscany at overshooting-pyroCu or above, much of it contiguous class 3 (resilient).

**Assessment.** Taking GRAF as correct, this is over-flagging by pyflam, and it is the
difference an operational user would notice first. It is consistent with the known weakness
recorded in `pyroconv.profile_diagnostics`: the mixed-layer stability gate has poor
specificity. Note however that on the *hybrid* path the ML gradient is a genuine least-squares
fit (~13 levels inside the layer today, 100 % fit support at peak of day), so the documented
pressure-level proxy weakness does **not** explain it. The cause is unidentified. Candidates,
untested:

1. threshold calibration — the LCL/ABL and cap γθ thresholds are paper values, not locally
   validated, and class 3 requires only LCL/ABL < 1.0 plus a weak cap;
2. the Rib ABL running deep at midday, lowering LCL/ABL across the board;
3. interpolate-then-threshold on the 2.2 km grid spreading marginal cells (see §4);
4. a genuinely different ladder or threshold set in the GRAF implementation.

### 3.2 Evening — direction of the disagreement is inverted

GRAF's **peak** on both forecast days is **18Z**: on 27/07 a broad yellow band across Tuscany
with dark-red cells near Firenze; on 28/07 a dark-red cluster over southern Tuscany / Umbria.
pyflam at those hours went to **zero**.

This one turned out to be a pyflam defect, not a forecast difference — see §5. The corrected
product no longer claims a quiet evening; it declares those hours unclassifiable. **The
underlying disagreement remains**: GRAF resolves a column at 18Z and pyflam does not.

---

## 4. Comparison caveats

- **Domain crop.** GRAF's largest 26/07 feature is mostly outside the pyflam bbox. Extent
  statistics are not comparable outside the Tuscany overlap.
- **Resampling order.** pyflam interpolates the *diagnostics* to 2.2 km and then applies the
  class thresholds; GRAF thresholds at native 6.5 km. Interpolate-then-threshold and
  threshold-then-interpolate differ near a threshold boundary, so class edges will not match
  even under identical physics. This affects edges, and cannot explain a 65-vs-5 % extent gap.
- **Time convention — unresolved.** GRAF's panels are labelled `Z`. If they were in fact local
  (CEST = UTC+2), their 18:00 would be 16Z and the evening offset against the pyflam 15Z peak
  would largely vanish. Against this reading: their day-0 peak at 15Z matches pyflam's day-0
  peak at 15Z, which would not survive a 2 h shift. Treat as an open question to put to the
  GRAF authors, **not** as an explanation.
- **Single case.** One run, one region, three days. Nothing here is a skill estimate.

---

## 5. Correction applied — nodata separated from class 0

### 5.1 The defect

`classify_profile` allocated its output with `np.zeros` and used `continue` for any column it
could not classify — no usable profile, a non-finite diagnostic, or an ABL below `abl_min_m`.
Those cells kept the value `0`, which the legend renders as **"Surface plume"**. The map could
therefore not distinguish *"the ladder ran and found no significant cloud development"* from
*"the ladder could not run"*.

Quantified on the 26 Jul 00Z run, POTENTIAL panel, land cells:

| | painted class 0 | of which unclassifiable | genuine surface-plume |
|:--|--:|--:|--:|
| 27/07 15Z | 6.1 % | 1.4 % | 4.6 % |
| 27/07 18Z | 99.2 % | **99.2 %** | **0.0 %** |
| 28/07 18Z | 97.2 % | **97.2 %** | **0.0 %** |

At 18Z on both forecast days, **not one cell** was a genuine surface-plume diagnosis. The Rib
ABL falls below the 600 m floor across effectively all land — physically plausible, since 18Z
is 20:00 local, near sunset — every cell is dropped, and the product drew the entire domain in
the "no significant convection" colour. It was not forecasting quiet evenings; it was failing
to classify them and rendering that failure as quiet.

This is the mechanism behind §3.2, and it inverted the comparison: exactly where GRAF placed
its peak, pyflam appeared to forecast nothing.

### 5.2 The change

| File | Change |
|---|---|
| `pyflam_gui/core/pyroconv.py` | New `PYROCONV_NODATA = -1`. `classify_profile` now allocates `np.full(..., PYROCONV_NODATA)`, so unclassifiable cells carry the sentinel instead of `0`. |
| `tests/pyroconv_daily.py` | Daily panels render `-1` in grey (`NODATA_COLOR`), with an `n/c — not classifiable` legend entry. |
| `scripts/forecast_3day_tuscany.py` | Same rendering in the 3-day figures. The **sea sentinel moved from `-1` to `-2`** — it previously collided with the new nodata value. `dist_table` gained an `n/c` column so rows sum to 100 across classes 0–4 plus n/c. |
| `scripts/forecast_3day_province_report.py` | Province peak class computed with `initial=-1`, so a province with no classifiable daytime column reports `-1` rather than reading as class 0. |

**Preserved distinction:** a fuel-gated cell below the 10 MW/m threshold stays class **0**.
That *is* a diagnosis — the atmosphere was assessed and the fire power was insufficient — and
must not become nodata. Verified directly:

| Case | Result |
|---|---|
| classifiable column | real class (4) |
| below the 10 MW/m fuel gate | **0** — genuine surface plume |
| ABL below the 600 m floor | **−1** — nodata |
| `valid=False` column | **−1** — nodata |
| non-finite ML gradient | **−1** — nodata |

The GeoTIFF exports already declared `nodata=-1`, so the rasters become self-consistent with
their own metadata for the first time; no consumer change was needed there.

### 5.3 What the correction does and does not do

It does **not** move a single class boundary, and it changes no forecast. Every cell that was
classified before is classified identically now. What changes is that cells which were never
classified stop being drawn as the lowest class. The 18Z panels go from "uniformly quiet" to
"uniformly blank", which is what the pipeline actually knows.

---

## 6. Open items for the follow-up evaluation

*Superseded by §11 — kept as the record of what was open before the investigation in §§7-10.*

1. **The extent gap (§3.1) is unexplained and is the substantive divergence.** The evening
   issue is now closed; this one is not. Priority.
2. **Recover the evening hours.** The 600 m `ABL_MIN_M` floor plus the Rib ABL definition
   removes the entire domain at 18Z, exactly when GRAF still resolves columns and when
   evening plume-driven events occur. Worth testing whether a residual-layer or parcel-based
   depth keeps those columns alive, and what GRAF uses.
3. **Confirm the GRAF time convention** (§4) before drawing any conclusion about the evening.
4. **Threshold-then-interpolate** as an alternative to the current order, to remove the
   resampling artefact from future comparisons.
5. **Repeat over more cases.** One run proves nothing about skill; the day-0 agreement and the
   day-1/2 extent gap could both be coincidental.

---

## 7. Evidence from the GRAF sonde campaign

The repository already vendors the ambient sondes from the campaign behind the method:
Castellnou Ribau et al. (2025), *AMT* 18(24), 7805-7831 (Zenodo 15264835) — eight environment
profiles, seven fires, Catalonia and Chile. Running our own diagnostics on them:

| Fire | launch UTC | ambient ABL | verdict under the 600 m gate |
|:--|--:|--:|:--|
| CasablancaIII08 | 22:27 | **202 m** | discarded, unclassified |
| CasablancaIII10 | 20:56 / 22:30 | **391 / 368 m** | discarded, unclassified |
| Rojals | 11:34 | **205 m** | discarded, unclassified |
| SelvadelCamp | 15:50 | **247 m** | discarded, unclassified |
| GranjaEscarp | 16:58 | 637 m | 4 — deep pyroCu/pyroCb ✓ |
| Martorell | 16:42 | 651 m | 4 — deep pyroCu/pyroCb ✓ |
| SantaAna | 16:41 | 959 m | 4 — deep pyroCu/pyroCb ✓ |

Our ambient ABL reproduces the campaign's own `ambient_ABL_m` to within a few metres
(202/202, 372/391, 651/651, 211/205, 955/959), so this is fidelity, not a computation error.

Two results follow. **The ladder works on the fires it can see** — all three columns clearing
the gate are correctly called deep pyroCu/pyroCb. **The gate discarded five of eight real
pyroconvective wildfires**, whose observed fire-induced boundary layers reached 451-2850 m.
The Catalan launches cluster at 15:50-16:58 UTC (17:50-18:58 local) and the Chilean at
20:56-22:30 UTC: the campaign sampled late afternoon and evening, the window this product
was blanking.

### 7.1 Tested and rejected: the fire-induced depth as the ladder's scale

Hypothesis: substitute fireABL for the ambient ABL in the LCL/ABL ratio, since fireABL/ambient
ran 1.7-20.2 across the sondes. Reclassified from the saved 2026-07-26 rasters (`noshear`
ladder for both baseline and variant, since `shear_dist` is not exported):

| 15Z, % land ≥ pyroCu | ambient depth | fire depth |
|:--|--:|--:|
| 26/07 | 58.9 | 54.4 |
| 27/07 | 89.6 | 70.0 |
| 28/07 | 83.7 | 72.5 |

**Rejected for the extent gap.** A 10-20 point reduction leaves 54-72 % of Tuscany, and class 4
*rises* (9.4 → 12.6 %): fireABL > ambient ABL always, so the ratio falls and cells move *up*
the ladder. The mechanism runs the wrong way. Changing the scale without adding a penetration
condition cannot fix this — which sharpened the question answered in §8.

## 8. The source method, read directly

Castellnou et al. (2022), *JGR-Atmos* 127, e2022JD036920, sections 2.1-2.4.

**§2.4.2 — ABL determination.** Bulk Richardson with **Ri_b > 0.33** (Zhang et al. 2014),
computed *"starting above the surface layer, estimated to be roughly **400 m AGL**"* — an
explicit departure from Zhang's 200 m, because the sonde launch layer is dominated by plume
indraft. **No minimum-depth gate is specified anywhere.**

**§2.4.1 — what the classification actually consumes:**

> *"The classification was based on ABL stability, **plume characteristics** (flattened,
> overshooting, PyroCu, PyroCb), **plume stages** (surface, penetration, or deepening), and
> turbulence position (height) on top of ABL."*

Plume characteristics and stage are *observations of an actual plume*. Our POTENTIAL map has
neither: it evaluates the ambient column and assumes a pyroCu-capable fire everywhere. **We
implemented the atmospheric subset of a classifier whose discriminating inputs are the
fire-side ones.** Saturating on a well-mixed afternoon is the expected behaviour of a
4-variable projection of an ≥8-variable classification applied to every grid cell — not a
mis-tuned threshold. This is the explanation for §3.1.

**§2.1.2 — the missing physical gate:**

> *"The jumps at the entrainment zone between ABL and free atmosphere (**Δq, Δθ**) and lapse
> rates on the free atmosphere (γq, γθ) **assess the ability of a parcel to penetrate above ABL
> and achieve free convection**."*

pyflam uses γθ alone. **Δθ is absent from the codebase entirely** — our comments recognise
Castellnou treats it as a state variable separate from the gradient, but only to justify
excluding it from the gradient; it was never implemented as a diagnostic. Penetration ability
is precisely what should stop an ordinary afternoon from scoring pyroCu.

Two further variables the paper computes and we do not: **FireCAPE** (Potter 2005, their
Eq. 2 — fire-modified CAPE integrating the fire's θ excess; distinct from surface CAPE, which
this product deliberately omits) and **plume updraft vertical speed**. The paper also
distinguishes **fireLCL / fireABL / fireShear** from ambient; we have only fireABL.

### 8.1 Fidelity deviations found

| Item | Paper | pyflam | Status |
|:--|:--|:--|:--|
| Ri_b critical value | 0.33 | 0.33 (gridded) | ✓ |
| Ri_b critical value | 0.33 | **0.25** (`bulk_richardson_abl_height`, 1-D) | ✗ open |
| Ri_b search start | ~400 m AGL | 200 m (gridded) | ✗ open, documented rationale |
| Minimum ABL depth | none | 600 m | ✗ **removed, §9** |

The start height was tested on the cached run: raising 200 → 400 m lifts the median 18Z ABL
from 225 to 404 m (27/07) but changes the share below 600 m by under 0.5 points. It is a real
fidelity issue; it was **not** the cause of the evening blanking.

## 9. Correction 2 — the 600 m ABL gate removed

### 9.1 Rationale

`ABL_MIN_M = 600.0` ("below this ABL depth, held at surface plume") was a pyflam invention. It
has no basis in the source method (§8) and is contradicted by the campaign's own sondes (§7),
which it would discard in five of eight real pyroconvective fires. A gate that rejects the
regime the instrument campaign was built to sample cannot be defended.

### 9.2 The change

`ABL_MIN_M` is now **150 m**, matching the Rib solver's own plausibility clip
(`atmosphere._ABL_MIN_M`). It guards against a near-zero depth turning LCL/ABL into a
meaningless quotient and is explicitly no longer a discriminating gate. The `abl_min_m`
parameter is retained so callers wanting a stricter floor (e.g. `scripts/validation/`) pass
their own. All 66 classifier and sonde tests pass; the campaign sondes go from **3/8 to 8/8**
classifiable.

### 9.3 Effect on the 2026-07-26 Tuscany grid

| | % n/c, 600 m gate | % n/c, after | % ≥ pyroCu, 600 m gate | % ≥ pyroCu, after |
|:--|--:|--:|--:|--:|
| 26/07 09Z | 43.7 | **0.1** | 12.5 | 16.2 |
| 26/07 15Z | 2.1 | 1.2 | 58.9 | 59.3 |
| 26/07 18Z | 68.9 | **31.0** | 25.3 | 37.9 |
| 27/07 15Z | 1.4 | 0.0 | 89.6 | 89.8 |
| 27/07 18Z | 99.2 | **51.9** | 0.1 | 1.8 |
| 28/07 18Z | 97.2 | **58.2** | 0.0 | 0.0 |

The morning recovers almost completely and the evening halves, while midday coverage is
unchanged (89.6 → 89.8 %) — the correction fixes what it should and does not touch the extent
problem. Notably it is also far more conservative than the rejected fireABL substitution,
which put 40-61 % of land at ≥ pyroCu at 18Z.

**The residual evening n/c is a different defect.** It equals the share where the least-squares
mixed-layer dθ/dz is undefined (31.0 / 51.9 / 58.2 %) — there are fewer than `ML_FIT_MIN_PTS`
model levels inside a mixed layer that has stopped existing. That is arguably correct
behaviour, but it means the ladder's second diagnostic is **ill-posed in the evening regime**,
which is exactly the regime the campaign instrumented. A gradient defined on the residual layer
(mixed-layer theory, Vilà-Guerau de Arellano et al. 2015) is the likely resolution.

## 10. Items 2-5 implemented

### 10.1 Item 2 — Richardson parameters reconciled

`bulk_richardson_abl_height` used **Ri_c = 0.25**, matching neither the paper (0.33, after
Zhang et al. 2014) nor pyflam's own gridded path. Fixed to the shared `_RIB_CRITICAL`.

The 200 vs 400 m start height is **resolved as a deliberate divergence, not a defect**. Zhang
et al. recommend 200 m; Castellnou et al. raise it to ~400 m specifically because their sondes
are launched into the indraft converging on the plume base (Charland & Clements 2013). That is
a fire-sonde correction. `bulk_richardson_abl_height` is the in-plume path and keeps 400 m;
`bulk_richardson_abl_grid` runs on ambient forecast columns with no indraft and keeps Zhang's
200 m. Both constants now carry the reasoning and the measured sensitivity.

### 10.2 Items 3-4 — the two missing variables, implemented as diagnostics

New in `pyflam.atmosphere`, exported and carried through `regrid_diagnostics`:

| Function | Quantity | On the campaign sondes | Over Tuscany, 27/07 15Z |
|:--|:--|:--|:--|
| `entrainment_jump_grid` | Δθ across the entrainment zone | 0.07–3.55 K | p10/med/p90 = 1.95 / 2.59 / 3.69 K |
| `fire_cape_grid` | FireCAPE (Potter 2005) | 846–1747 J/kg | 1964 / 2192 / 2361 J/kg |

Both are **diagnostic only** — not wired into the ladder, so no class changes from them. That
restraint is deliberate, and the reason is the next finding.

### 10.3 The penetration gate does not yet discriminate — and why

The obvious use of Δθ is the source method's penetration test: can the fire's θ-excess clear
the jump? Measured over Tuscany at the hour with 89.8 % coverage, the ratio θ′/Δθ has
**p10 = 5.46** — every column clears it 5–10×. A gate at 1 would exclude nothing, so wiring it
in would be a no-op dressed as physics.

The reason is a scaling problem upstream. The reference fire's θ′ is **20.1 K**, while the
campaign *measured* in-plume excesses of **0.1–13.1 K** (median 3.6). Our reference parcel is
hotter than any fire the campaign instrumented. Sensitivity at 27/07 15Z:

| assumed θ′ | source | % land passing θ′/Δθ ≥ 1 |
|--:|:--|--:|
| 20.1 K | current 200 W/m² reference | 100.0 |
| 9.9 K | campaign maximum | 100.0 |
| 3.6 K | campaign median | 88.9 |
| 2.5 K | campaign lower range | **44.4** |

The gate becomes a real discriminator only at θ′ ≈ 2–4 K, i.e. at the *observed* magnitudes.
**Calibrating θ′ is therefore a precondition for using the penetration test at all**, and the
campaign supplies ten measured values to calibrate against. Caveat: the measured `theta_excess_K`
is an in-plume excess while `fire_parcel_theta_excess` returns a mixed-layer forcing for a
cell-averaged flux — they may not be the same quantity, which has to be settled before fitting.

### 10.4 Item 5 — the residual-layer gradient, and it works

`residual_layer_grid` references the parcel to the **coldest level in the lower column** (the
top of the nocturnal stable layer) rather than to the surface, and takes the residual top as
the first height exceeding it by 0.5 K. By day the θ minimum is at the surface and this reduces
exactly to the parcel mixing depth; only in the evening do the two separate.
`iconeu_diagnostics` now fits the mixed-layer gradient over `max(parcel_ml, residual_ml)`.

Classifiable land, the metric the evening divergence turns on:

| | before item 5 | after |
|:--|--:|--:|
| 27/07 15Z | 98.6 % | 100.0 % |
| 27/07 18Z | 48.1 % | **81.0 %** |
| 28/07 18Z | 41.8 % | **80.5 %** |

Together with §9 this closes the evening divergence as a *pyflam defect*: of the 99.2 % of land
unclassifiable at 27/07 18Z before this work, the ABL gate accounted for roughly half and the
ill-posed gradient for the rest. Whether the classes now produced there agree with GRAF is a
separate question, and untested.

### 10.5 Item 8 — the fire θ-excess: a dimensional error, and a disqualified gate

**The definitional question resolved first.** The campaign's `theta_excess_K` is
mean(θ in-plume) − mean(θ environment) over the lowest 200 m — a measured plume-vs-ambient
warm anomaly. `fire_parcel_theta_excess` predicts a temperature scale from a heat flux. They
are the same quantity *by construction of the existing calibration*: step 5 of
`scripts/ingest_inplume_sondes.py` substitutes the measured excess straight into the
encroachment, "so no fire intensity or flux parameterisation is needed". Calibration was
therefore admissible.

**A dimensional error in the conversion.** The formula computed

    theta' = F / (rho * w0)

which carries units of **J/kg, not kelvin** — the specific heat capacity is absent from both
the velocity and the temperature scale — and the result was added directly to a potential
temperature in `fire_induced_abl_grid`. For F = 200 W/m² it returned **20.1 K** where the
dimensionally correct Deardorff value is **0.23 K**, an overstatement of ~88×. The codebase had
already sensed something wrong here — the note on `_FIRE_PLUME_SCALE_M` records that "the
theta' FORM is the problem" — but attributed it to the choice of flux rather than to units.
The end-to-end fireABL validation (bias −15 m, MAE 400 m, r = 0.91) carried "one global
coefficient for the table's unit ambiguity", and r is invariant under constant rescaling, so
the fix does not invalidate it; the production reference path applied no such coefficient.

**Why calibrating the flux cannot work.** Even dimensionally corrected, θ* is the *turbulence*
scale of a convectively mixed layer — a few tenths of a kelvin — while a plume core is a
coherent buoyant structure the campaign measured at 0.1–13.1 K. No cell-averaged flux closes
an order-of-magnitude gap, because it spreads the fire's heat over ground that is mostly not
burning. The reference is therefore now **prescribed in kelvin** (`_REFERENCE_THETA_EXCESS_K
= 10.0`, the intense end of the measurements), and `fire_induced_abl_grid` takes
`theta_excess=` directly — the same route the campaign ingest already took.

**The gate is disqualified, and not by θ′.** Forced with each fire's *own measured* excess —
the best input obtainable — the encroachment reproduces observed fireABL at **151 % relative
error** across the four usable observations:

| Fire | measured θ′ | observed fireABL | predicted | error |
|:--|--:|--:|--:|--:|
| CasablancaIII10 | 4.38 K | 2850 m | 1769 m | 1081 m |
| GranjaEscarp | 4.67 K | 1592 m | 1634 m | 42 m |
| SantaAna | 9.86 K | 371 m | 1628 m | 1257 m |
| Martorell | 2.63 K | 1002 m | 3247 m | 2245 m |

You cannot calibrate the input of a model that does not track its output. **The blocker on the
penetration gate is not an uncalibrated θ′ — it is that the fireABL chain has no demonstrated
skill on the only data that can test it** (n = 4, mean |error| 1156 m on observations of
371–2850 m). Item 3 stays diagnostic for that reason, and the next question is the encroachment
model itself, not its forcing.

Effect on the published diagnostic: θ′ 20.1 → 10.0 K halves it. At 27/07 15Z the decoupling
median goes 2.9 → 2.08 and the penetration ratio 7.8 → 3.86 — still above 1 everywhere, so the
gate would remain a no-op even now.

## 11. Revised open items

| # | Item | Status |
|:--|:--|:--|
| 1 | Remove the unfounded 600 m ABL gate | **done, §9** |
| 2 | `bulk_richardson_abl_height` Ri_c 0.25 → 0.33; reconcile the 200/400 m start | **done, §10.1** |
| 3 | Implement **Δθ** (entrainment jump) | **done as diagnostic, §10.2** — gate blocked on the fireABL model's skill (§10.5), not on item 8 |
| 4 | Implement **FireCAPE** (Potter 2005) | **done as diagnostic, §10.2** |
| 5 | Define the ML gradient on the residual layer | **done, §10.4** — evening classifiable 48 → 81 % |
| 6 | Re-run the GRAF comparison once the above are in | now unblocked |
| 7 | Confirm the GRAF time convention (§4) | open, needs the authors |
| 8 | Calibrate the fire θ-excess against the campaign's measurements | **done, §10.5** — but it closed the question by *disqualifying* the penetration gate, not by enabling it |
| 9 | Regenerate the published products | **2026-07-26 done** (items 1-5 and 8; 48 class rasters unchanged, 96 θ′-dependent diagnostics updated). `docs/forecast_2026-07-25_3day*/` still predates everything, and there is still no provincial report for 26 July |
| 10 | **The decoupling ratio is dominated by its denominator in the evening** — see below | open |

### 11.1 Item 10 — the decoupling ratio's collapsing denominator

Flagged after the item-8 regeneration. The dry-pyrocloud diagnostic is fireABL / ambient ABL,
and removing the 600 m gate (item 1) admitted exactly the columns with the smallest
denominators. The two changes push the ratio in opposite directions, and the denominator wins:

| max decoupling ratio | before items 1-8 | after |
|:--|--:|--:|
| 26/07 18Z | 9.1 | **19.1** |
| 27/07 18Z | 9.7 | **19.4** |
| 28/07 18Z | 11.8 | **18.1** |

The maxima roughly **doubled even though θ′ halved**, because the median 18Z ambient ABL is
now ~230 m rather than the ≥600 m the gate used to enforce. The 18Z rows also read 100 % at
both ratio thresholds — no longer the small-denominator artefact corrected earlier, since the
classifiable share is now 47-82 % rather than 0.6-49 %, but the ratios themselves are being set
by a shallow ambient depth rather than by fire behaviour.

This is the same ambient-versus-in-fire depth problem as the evening class divergence (§12.6
item 2), surfacing in the dry diagnostic. It does not invalidate the θ′ correction, which moved
in the right direction and onto measured ground. It does mean **the evening decoupling panels
must not be read as fire behaviour** until the denominator question is settled — the candidates
being a fire-relevant depth in place of the ambient one, or reporting fireABL directly rather
than as a ratio.

Items 2 and 5 were corrections; 3 and 4 added physics taken from the published method rather
than calibration against GRAF output, so model independence is preserved. Item 8 would
calibrate against the campaign's *published measurements*, not against the GRAF product. **No threshold in
this codebase has been altered to match the GRAF product**, and none should be.

## 12. Synthesis — how pyflam differs from the GRAF model

A standing summary of every difference established so far, organised by kind rather than by
the order it was found. Sections 1-11 are the working record; this is the account to read
first.

### 12.1 The root difference

**They classify a plume; we screen an atmosphere.** Everything consequential follows.

Castellnou et al. (2022) sec. 2.4.1 lists the classifier's inputs: *ABL stability, plume
characteristics (flattened / overshooting / pyroCu / pyroCb), plume stages (surface,
penetration, deepening), and turbulence position on top of ABL.* Two of those four are
observations of a real plume. The pyflam POTENTIAL map has neither -- it evaluates the ambient
column and posits a pyroCu-capable fire in every cell. We implemented the atmospheric
projection of a classifier whose discriminating inputs are the fire-side ones. That is why we
saturate at midday and they do not, and no threshold adjustment reaches it.

### 12.2 Variables they use and we lack

| Variable | Role in the source method | pyflam |
|:--|:--|:--|
| Plume stage / characteristics | classification input | **unobtainable from NWP** — the core gap |
| Plume updraft vertical speed | where updrafts stabilise | unobtainable (they derive it from balloon ascent) |
| **Δθ** entrainment jump | sec. 2.1.2: *"assess the ability of a parcel to penetrate above ABL and achieve free convection"* | computed since 2026-07-26, **diagnostic only** |
| **Δq** moisture jump | same clause | absent |
| **FireCAPE** (Potter 2005) | critical variable, their Eq. 2 | computed since 2026-07-26, **diagnostic only** |
| fireLCL, fireShear | in-fire counterparts of the ambient levels | absent; we have fireABL alone |

### 12.3 Same quantity, different definition

| | GRAF | pyflam | Measured effect |
|:--|:--|:--|:--|
| Ri_b critical value | 0.33 | 0.33 ✓ (1-D path was 0.25 — fixed) | — |
| Ri_b search start | 400 m (in-plume indraft correction) | 200 m (Zhang, ambient column) | 18Z median ABL 404 vs 225 m |
| ABL identification | height of maximum RH, *supplemented* by Ri_b | Ri_b alone | **1.09x on well-mixed columns, 3.8-9.4x on collapsed ones — §14** |
| LCL | MetPy iterative, surface parcel | Bolton, surface parcel | negligible |

### 12.4 Gates pyflam invented

| Gate | Basis in the method | Status |
|:--|:--|:--|
| `ABL_MIN_M` = 600 m | none | **removed** — falsified by 5 of 8 campaign sondes at 202-391 m |
| `rh_top_moist` 80 → 60 % | none; tuned to GRAF's class-4 fraction on 2026-07-15 | **in place**, and applies only to classes 3-4 — class 2 carries no moisture gate at all |
| 10 MW/m fireline-intensity fuel gate | Tedim et al. 2018 | ours alone; no GRAF counterpart |

### 12.5 Pipeline and presentation

- **Order of operations.** They threshold at native 6.5 km; we interpolate the diagnostics to
  2.2 km and threshold there. Affects class edges, not extent.
- **Taxonomy.** Their scale is 1-4 with pyroCb = 1; ours 0-4 with pyroCb = 4 — inverted. They
  have neither a surface-plume class nor a not-classified class, so **white on a GRAF map is
  undefined**: it may mean no significant convection or no usable column. That is exactly the
  conflation removed from pyflam in §5, and it means statements in §2 about both products
  "agreeing on blank nights" rest on an assumption their legend does not license.
- **Palette.** darkred / gold / yellow / lightgreen against ColorBrewer RdYlGn; selectable via
  `PYROCONV_PALETTE=graf`. Cosmetic, but note gold and yellow are adjacent hues, so their
  resilient/overshooting pair cannot be reliably separated by eye in a compressed image --
  which is a limit on any visual reading of their published maps, including the ones in §3.

### 12.6 Ranked by consequence

1. **Missing fire-side inputs** — causes the 65-90 % vs isolated-patches extent gap, now
   measured independently of their maps as a **+1.4 to +1.8 class bias** against labelled
   fires (§17). Unresolved. Δθ and FireCAPE are the candidate proxies; the penetration test
   is a no-op at the reference θ′ (§10.3), and firepower is now computable per fire (§16)
   but has not yet been made to condition the classification.
2. **Ambient vs in-fire levels** — causes the evening disagreement. At 18Z the ambient ABL is
   230 m against a 1105 m LCL, so 98 % of columns exceed the overshooting threshold and fall
   to convection plume, while GRAF reaches pyroCb.
3. **The RH relaxation** — an uncalibrated deviation on classes 3-4, fitted to a single day and
   never revisited.
4. **Ri_b start height** — real, small.
5. **Resampling order, palette, class numbering** — presentational.

### 12.7 What remains unknown

*(Partly superseded: §16 shows the fire behaviour data is public, and §17 measures the extent
gap against labels rather than against their maps.)*

- Whether GRAF applies an additional mask or condition in their gridded product. Their
  documented physics offers FireCAPE and the penetration test as candidates; we have no
  confirmation.
- Their time convention: panels are labelled `Z`, but a local reading would move the evening
  comparison materially (§4).
- Any quantitative GRAF field. The entire extent comparison is ordinal, read from an image.

### 12.8 The one point of genuine corroboration

Day 0: same peak hour (15Z), same class ceiling, same Apennine-crest band running east toward
Marche. Two independent implementations of the same published ladder agreeing on structure and
timing. **This is evidence precisely because pyflam has not been tuned to GRAF output**, and it
is the reason item 8 must calibrate θ′ against the campaign's published *measurements* rather
than against their maps.

## 14. The ABL definition — the evening divergence, explained

Established 2026-07-27 from the full campaign dataset (Zenodo 17886250, 61 ambient sondes,
2021-2025) and the AMT supplement.

### 14.1 They do not define the ABL the way we do

Castellnou Ribau et al. (2025) sec. 2.6:

> *"The height of the maximum RH value is used as a criterion to estimate the height of the
> atmospheric boundary layer. This criterion is based on the observation that specific humidity
> tends to be well mixed in the convective boundary layer... reaching a peak at the inversion
> level. Above this inversion, the air becomes drier and warmer, resulting in a decrease in RH."*

A **moisture** criterion, where `bulk_richardson_abl_grid` is a **dynamic** one. Implemented as
`atmosphere.max_rh_abl_grid` for comparison; the classifier still uses the Rib depth.

### 14.2 The two agree when mixed and diverge when collapsed

Across 59 ambient campaign sondes with both definitions resolvable:

| Rib regime | n | max-RH / Rib |
|:--|--:|--:|
| Rib < 300 m (collapsed, stable) | 36 | **3.83x** |
| Rib 300-800 m | 14 | 1.85x |
| Rib > 800 m (well mixed) | 9 | **1.09x** |

corr(log Rib, log ratio) = -0.41. The same split appears on the ICON-EU Tuscany grid:

| 27/07 | Rib | max-RH | ratio | LCL/Rib | LCL/max-RH |
|:--|--:|--:|--:|--:|--:|
| 15Z (well mixed) | 1766 m | 2163 m | 1.14 | 1.00 | 0.86 |
| 18Z (collapsed) | 225 m | **2198 m** | **9.37** | 4.37 | **0.55** |

### 14.3 What it explains, and what it does not

**The evening divergence (§3.2): explained.** At 18Z the median column moves from LCL/ABL 4.37
-- far above the overshooting band, hence class 1 convection plume -- to **0.55**, inside the
resilient-pyroCu band. That is GRAF's evening peak, produced by their own published definition
with nothing tuned. The mechanism is visible in the numbers: the max-RH ABL is essentially
unchanged from 15Z to 18Z (2163 -> 2198 m) because the *moisture* structure retains the
daytime mixed layer after the thermal and dynamic structure has collapsed. It is, in effect,
detecting the residual layer -- the same physics as item 5, reached by a different route.

**The midday over-extent (§3.1): not explained.** At 15Z the definitions agree to 14 %, and the
change moves the ratio *down* (1.00 -> 0.86), shifting overshooting toward resilient rather
than reducing coverage. Consistent with §12.1: the midday gap is about fire-side conditioning,
not the ABL.

**Open item 10 is also implicated.** The decoupling ratio's collapsing denominator is a Rib
artefact; a max-RH depth would not collapse at 18Z, so switching the denominator would remove
most of that problem too.

### 14.4 Caveat

The source applies this criterion **visually**, with a human rejecting spurious maxima.
`max_rh_abl_grid` is an automated analogue: it searches above a surface-layer floor and
requires RH to fall by a set amount above the peak, without which a monotone profile returns
its top level. On the campaign sondes it still picks implausibly low heights in a minority of
columns. It should not be swapped into the classifier without a hand comparison against
plotted profiles.

## 15. Three definitional gaps -- the pattern to put to GRAF

Each divergence traced to root cause has the same shape: **pyflam implements the atmospheric
skeleton of the method, without the fire-side or observational inputs the method assumes.**

| # | Gap | Consequence | Status |
|:--|:--|:--|:--|
| 1 | **ABL by max-RH** vs bulk Richardson | evening classes: convection plume vs resilient pyroCu | measured, §14 |
| 2 | **In-fire levels** (fireLCL / fireABL / fireShear) vs ambient | our modelled fireABL reproduces the measured in-fire ratio in 5 of 7 sondes, but classes do not move because the cap clause dominates | measured, §12.2 |
| 3 | **Fire-power conditioning** absent | midday coverage 65-90 % against isolated patches | unresolved, the larger gap |

This pattern is a better opening with GRAF than any single number, because it is a question
about *method*, not about their output: we have their published equations and their sondes, and
what we lacked was the fire-side data their classification consumes.

**Correction (2026-07-27): most of that data is already public.** The statement above --
written when neither the AMT supplement nor the Zenodo archive carried fire behaviour -- was
wrong. The **Wildfire Data Portal** (wildfiredataportal.eu, EWED/ODET, EU co-funded) publishes
it per fire; see §16. What genuinely remains unavailable is the *plume-side* observation --
plume stage and plume updraft speed -- which is derived from the in-plume sondes and the
analyst's reading, not tabulated. That is the narrower thing to ask about.

## 16. The Wildfire Data Portal -- the fire-side data, and the firepower distribution

`wildfiredataportal.eu` exposes a REST catalogue (`/wp-json/wp/v2/fire`) of **27 fires**, of
which **23 carry populated fire behaviour**: start hour (UTC), surface affected, rate of spread
mean/max, **burn ratio mean/max in ha/h**, flame lengths head and flanks, torching fraction and
spotting distance, plus a perimeter KMZ and, per fire, a run of the Wageningen **CLASS**
mixed-layer model -- which makes that collaboration concrete rather than inferred.

Burn ratio in ha/h is exactly the `dA/dt` of Tory & Kepert's appendix D, so firepower follows
without any assumption of our own:

    FP = alpha * h * w_a * dA/dt,   alpha = 0.7, h = 15 MJ/kg, w_a = 1.25-4 kg/m2

| Fire | burn ratio | firepower (w_a 1.25-4) |
|:--|--:|--:|
| Guissona | 6000 ha/h | 219-700 GW |
| Varnavas | 2750 ha/h | 100-321 GW |
| Katsimidi | 2000 ha/h | 73-233 GW |
| Pauls, El Valle 2 | ~780 ha/h | 28-91 GW |
| Santa Coloma de Queralt | 450 ha/h | 16-53 GW |
| Martorell | 90 ha/h | 3.3-10.5 GW |
| Rojals (prescribed) | 0.2 ha/h | ~0 GW |

Against published PFTs of 12-1240 GW and our Tuscany field at 17-166 GW, the distribution
straddles the decision boundary rather than sitting to one side -- which is what a useful
conditioning variable should do, and is the quantity the POTENTIAL map lacks entirely.

**Guissona is excluded from any fit**: its 6000 ha/h is internally inconsistent with a reported
ROS max of 12 m/h, and being the largest firepower in the set it would otherwise dominate.
**Rojals** is a natural negative control -- a 0.97 ha prescribed burn at ~0 GW that should clear
no threshold in any column.

`docs/graf_fire_firepower_join.csv` joins this to the sonde archive: 22 sonde-fire pairs over 15
distinct fires, each with date, start hour, area, ROS, burn ratio, firepower bounds and the
sonde-derived Rib ABL / max-RH ABL / LCL. The §14 regime split reappears across these 15
independent fires -- Rib pinned at 200-350 m in most rows while max-RH runs 300-3120 m.

Remaining for the falsifiable PFT test: the column above sonde apex (631-7555 m here) spliced
from ERA5 or ICON to reach the -20 C level. Dates and coordinates are known for every row, so
that step is mechanical.

## 17. Ladder validation against GRAF's observed fire classes

The strongest result here, and the one that no longer depends on reading their maps.

The Wildfire Data Portal labels every fire with GRAF's own observed pyroconvection class.
Running the pyflam ladder on each ambient campaign sonde and comparing to that label gives a
direct test of the classifier against ground truth on real fires -- 26 sondes, ~15 distinct
fires.

| ABL definition | exact | within 1 class | **mean bias** |
|:--|--:|--:|--:|
| bulk Richardson (current default) | 2/26 (8 %) | 11/26 (42 %) | **+1.81** |
| max-RH (theirs, §14) | 2/26 (8 %) | 16/26 (62 %) | **+1.38** |

The error is systematic, not scattered: fires observed as *Convective plume* are repeatedly
classified 4 (Vilanova de Meia, Vega Honda, Junquillos, Granyena, San Patricio). The single
Deep pyroCu/pyroCb fire is classified correctly, which is what an over-predicting classifier
does.

**This converts §3.1 from an eyeball comparison against a JPEG into a measured defect.** The
midday over-extent is not an artefact of comparing the wrong panel or misreading their
colours -- the ladder runs 1.4 to 1.8 classes hot against labelled fires.

**It independently supports §14.** Their ABL definition cuts the bias by 0.43 classes and
lifts within-1 agreement from 42 % to 62 %, tested against labels rather than against their
product. It helps materially; it does not close the gap.

### 17.1 Why the bias is a lower bound

Three caveats, of which the first *strengthens* the finding:

1. **The label is the fire's peak class; the sonde samples one phase.** Predicting high scores
   as correct against a peak label while being wrong about that moment. The measured +1.81 /
   +1.38 therefore **understates** the true over-prediction.
2. **The POTENTIAL framing explains part of it, but not enough.** Our ungated classification
   assumes a pyroCu-capable fire, so over-prediction on a 0.97 ha prescribed burn is by
   construction. But Granyena (6.9-22 GW), Vega Honda (6.2-19.8 GW) and Junquillos
   (1.5-4.7 GW) were substantial fires observed as convective plumes, and we call them class 4.
3. **n = 26 sondes but only ~15 distinct fires** -- San Patricio x3, Pauls x3, Granyena x2 --
   so the effective sample is smaller than it looks.

## 18. PFT against firepower on real fires -- inconclusive, and why

With firepower from §16 and PFT from the spliced sonde + ERA5 column:

| observed class | n | clears its PFT |
|:--|--:|--:|
| Surface plume | 3 | 0 |
| Convective plume | 8 | 0 |
| PyroCu / Overshooting / Resilient | 9 | 0 |
| **Deep pyroCu / pyroCb (Guissona)** | 1 | **YES** -- 139 GW against 219-700 GW |

19 correct negatives, 1 correct positive, 0 false positives. That reads well, but the test is
weak: with one positive case it cannot demonstrate sensitivity, and a threshold that rarely
fires achieves the same score trivially.

Two methodological traps were found and fixed in the course of running it, both of which had
produced wrong conclusions:

* **ERA5 time.** The first run built every request from the *fire's start hour* and spliced it
  onto a sonde launched hours later -- 4 h for Martorell, 6 h for SCQ, 5 h for Manuel
  Rodriguez. Rebuilt at each sonde's own launch datetime from Table S1. This alone removed the
  run's only false positive (Granyena).
* **Mixed-layer window.** `pyrocb_firepower_threshold_grid` averaged over an absolute
  `z <= 500 m`, which selects no levels on a sonde starting at 840 m AGL (Guissona), returning
  a nan PFT that reads as "no firepower suffices". Guissona was reported on that basis as a
  false negative of the method; it is not -- corrected, it clears its threshold, as an observed
  pyroCb should.

**Santa Coloma de Queralt remains unresolved.** Its sonde is from 24 July at 19:17, but
Castellnou Fig. 6d places the pyroCb on the **25th** (SCQ51); the 24th was the pyroCu day. The
portal's fire-level label attributes the 25th's behaviour to a 24th profile. Testing SCQ needs
the 25 July sondes, which are separate Table S1 rows. It is the only fire with per-sonde
types (Castellnou et al. 2022), so it is the one worth doing properly.

## 21. Fuel load for the firepower estimate -- two routes reconciled

Appendix D's firepower needs `w_a`, the fine fuel available for burning. Tory & Kepert give
1.25-4 kg/m2 (Australian eucalypt) and use 9.4 for Chisholm (boreal). A 3.2x band is wider than
the margins it is used to judge -- for Santa Coloma de Queralt it spanned the entire distance
between clearing and not clearing its threshold -- so it had to be pinned.

### 21.1 A route that does not work

Deriving `w_a` from the portal's published flame length via Byram (`FLI = 258 L^2.17`, then
`w_a = FLI / (H ROS)`) **fails**, and fails worst on the fires that matter. The relation is
calibrated for surface flames of roughly 0.5-5 m; the portal reports 15-60 m plume heights for
the plume-dominated fires, and extrapolating gives 92-1863 MW/m against the 4-50 MW/m
Castellnou et al. (2022) actually observed. Inverted, that yields fuel loads of 93-763 kg/m2
against a real range near 2. Where ROS is also small (Guissona, 12 m/h) the division explodes.
Only low-intensity fires land sensibly (Patagual 0.68, Junquillos 0.42 kg/m2).

### 21.2 Two routes that agree

**From observed fire behaviour.** Inverting Byram through *published* FLI rather than
flame-length-derived FLI:

| fire / moment | FLI | ROS | w_a |
|:--|--:|--:|--:|
| Martorell M12 (pyroCu) | 11 MW/m | 1350 m/h | 1.96 kg/m2 |
| SCQ51 (pyroCb) | 18 MW/m | 3300 m/h | 1.31 kg/m2 |
| SCQ31 (pyroCu peak) | 47 MW/m | 3300 m/h | 3.42 kg/m2 |
| campaign mean | 8.5 MW/m | 1500 m/h | 1.36 kg/m2 |

Range 1.31-3.64, median ~2.0 kg/m2.

**From the fuel models.** Scott & Burgan loads for Mediterranean shrub and timber-understory
fuels, taking the **full grass/shrub load** -- 1h + 10h + live herb + live woody, not just the
fine dead and herbaceous components -- and applying an **empirical +30 % Mediterranean-basin
adjustment**:

| model | full load | +30 % |
|:--|--:|--:|
| SH2 | 1.70 | **2.21** |
| SH5 | 1.93 | **2.51** |
| SH7 | 2.73 | **3.56** |
| TU5 | 2.47 | **3.21** |

The two routes land on the same value from independent directions.

> **Provenance of the +30 %.** This is an empirical adjustment from operational use in the
> Mediterranean basin (C. Foderi), not part of Scott & Burgan, applied as an upper bound on the
> model load. It is recorded as a stated assumption rather than folded into a constant, because
> it is a local calibration and should travel with that caveat.

### 21.3 What was wrong in the first cross-check

An earlier version of this analysis reported Scott & Burgan giving 0.2-0.9 kg/m2 and concluded
the models understate Mediterranean fuel by 3x. That was an error on this side, twice over: the
sum omitted the 10h and live-woody components that carry most of the load in shrub fuels, and a
tons/acre conversion was applied to loads already stored in lb/ft2. Corrected, the models do not
conflict with the fire-behaviour route -- they corroborate it.

### 21.4 Effect

`w_a = 2.0 kg/m2` (band 1.5-2.5) is now used in `docs/graf_fire_firepower_join.csv`,
sourced from two agreeing routes rather than an assumed range. Uncertainty falls from 3.2x to
about +/-25 %, and the consequence is to **close** the escape hatch rather than widen it:

* **Guissona** 350 GW (262-438) against a 139 GW threshold -- clears comfortably, as an observed
  pyroCb should.
* **Santa Coloma de Queralt** ~26-33 GW against 66 GW -- still short by ~2x, and the narrowed
  band no longer spans the boundary. The fuel load is no longer a candidate explanation; the
  remaining one is that burn ratio is an hourly average while the pyroCb-generating head fire is
  a short burst, which appendix D acknowledges ("averaged over the time period dt").

## 21. The heat-impulse term: identified, and not yet usable

The single most consequential threshold in the ladder is `lcl_ratio_overshoot_max = 1.60`, an
upper bound on LCL/ABL with no citation in Castellnou et al. (2022). It terminates **16 of the
26** labelled campaign sondes before any later gate is consulted -- more classifications than
every other diagnostic combined.

### 21.1 What it actually is

The MARI reference port, which is where the value comes from, says so in a comment on that
branch:

> *"LCL sopra ABL ma non troppo lontano; **richiede un impulso di calore del fuoco, non
> ricavabile dal solo ERA5**."*
> -- the LCL above the ABL but not too far, which requires a heat impulse from the fire, not
> obtainable from ERA5 alone.

So 1.60 is not a geometric fact about the atmosphere. It is a **stand-in for the fire-power
term**, placed where the fire-side data was missing, by an author who knew that is what he was
doing. The same port heads its threshold block *"soglie operative iniziali"* -- initial
operational thresholds -- and says of the cap pair *"non sono soglie universali"*.

### 21.2 It is uncited and it works

Tested against the observed classes on the adaptive ladder:

| ceiling | exact | within 1 | bias |
|--:|--:|--:|--:|
| **1.60 (current)** | **11/26** | 20/26 | **-0.08** |
| 3.00 | 9/26 | 21/26 | +0.04 |
| 5.00 | 9/26 | 22/26 | +0.08 |
| removed | 4/26 | 21/26 | +0.58 |

(This sensitivity test was run on the 26-sonde sample; §21.3 onward uses 27, after the name
matcher recovered one more fire. The bar re-measures at 12/27 there.)

Removing it collapses exact agreement and pushes the ladder half a class hot. The labels come
from a 2025-26 publication and the value from an earlier port, so this behaves as out-of-sample
corroboration rather than a fit.

It misclassifies individual events -- Santa Coloma de Queralt, an observed pyroCb, sits at ratio
12.8 and is called class 1 -- but raising the ceiling to admit it costs seven other fires. The
ladder is trading errors, not simply making one.

### 21.3 The explicit fire-power term, tested properly

pyflam has every piece needed to replace the proxy with the physics: a PyroCb Firepower
Threshold verified against 5 of 6 published worked cases (§18), fuel load pinned to
~1.5 kg/m2 by three independent routes (§20), and head-fire length measured from 20 fires. The
remaining input was `dA/dt`, and §21.4 below closes it.

An earlier version of this section reported that every explicit variant lost to the proxy. **That
test was broken in two independent ways and its numbers are withdrawn**:

* Firepower resolved for only 11 of 26 sondes, because the isochrone parser missed three whole
  classes of file and compared UTC launch times against local perimeter times (§21.4).
* More seriously, **the PFT resolved for 1 sonde in 27**. It was computed from the sonde alone,
  and a sonde topping out at 3-5 km cannot reach the free-convection height or the -20 degC
  cloud-top level the threshold is defined at, so it returned nan almost everywhere. The gate
  was not losing the comparison, it was **never firing**. Splicing ERA5 above each sonde top
  (`splice_pft.py`) resolves PFT for 27/27.

With both faults fixed, on all 27 labelled sondes, adaptive ladder, `w_a` = 1.49 kg/m2:

| variant | exact | within 1 | bias |
|:--|--:|--:|--:|
| ceiling 1.60 only (the bar) | 12/27 | 20/27 | +0.04 |
| **+ firepower gate, fire-level peak** | **14/27** | **21/27** | -0.07 |
| + firepower gate, per-moment isochrone | 13/27 | 20/27 | -0.15 |
| ceiling removed, per-moment gate only | 6/27 | 21/27 | +0.48 |

The gate now *beats* the bar. **It should still not be adopted**, because of what the firing
pattern shows.

### 21.4 Per-moment firepower from the isochrones

The portal publishes perimeter KMZs per fire. Extracting `dA/dt` at each sonde's launch minute
required three fixes, all found by inspection (`scripts/isochrone_firepower.py`):

1. **A sixth timestamp format** (`20220717 1630`).
2. **Files with no timestamp in `<name>` at all**, carrying it instead in the ArcGIS attribute
   table rendered into `<description>` HTML, under `DiHo` or `FeHo`.
3. **Perimeter times are local, not UTC.** Nothing in the files says so. The campaign sondes
   settle it: Guissona's launched at 15:59 UTC into an already-convecting plume whose isochrones
   run 17:13-19:29, so under a UTC reading the sonde precedes every mapped perimeter of the fire
   it was sampling. Reading them as local puts it inside the span, and does the same for
   Patagual, Junquillos and Vega Honda. Converted per fire via `zoneinfo` (DST-correct), zones
   assigned from portal coordinates, which fall in four well-separated clusters.

Coverage went from 11/26 to **22/27**; growth rates from 18 fires to **24**, still agreeing with
the portal's independently reported peak burn ratio (guissona 7869 vs 6000 ha/h, lavrio 860 vs
795, martorell 103 vs 90). The 5 unresolved sondes belong to fires whose files hold only a final
perimeter -- no time series exists to extract.

### 21.5 Why the improvement is not real: a timescale mismatch

The gate fires on **3 of 27 sondes**:

| sonde | obs | raw | PFT (GW) | FP peak | FP moment | effect |
|:--|--:|--:|--:|--:|--:|:--|
| pauls | 2 | 3 | 518 | 34 | 4 | demoted, **correct** (both variants) |
| pauls | 2 | 4 | 1278 | 34 | 17 | demoted, **correct** (both variants) |
| guissona | 4 | 4 | 139 | 261 | 20 | per-moment demotes, **wrong** |

That arithmetic is the whole table above: the peak variant fires twice and is right twice
(12 -> 14); the per-moment variant fires three times, right twice and wrong once (12 -> 13).
**The entire apparent improvement is two sondes of a single fire.** At n = 27 that is noise.

The failure is the informative half. Guissona is the sample's one observed PyroCb. Its sonde
was launched into an isochrone interval growing at ~460 ha/h (20 GW) against a fire peak of
6000 ha/h (261 GW) -- **it blew up after the launch**. The per-moment gate correctly reports
that the fire was not yet capable of pyroCb, and is scored wrong for it, because the label
describes the fire's eventual behaviour. Pauls rewards the gate and Guissona punishes it for the
same reason.

So the fire-level peak does not "work better" because it is more physical. It works better
because it matches the *timescale of the label*.

**A moment-level gate cannot be validated against a fire-level label.** This is not a data
volume problem, and no quantity of additional isochrones fixes it. Testing per-moment
conditioning requires per-moment observed classes -- the plume behaviour at each launch -- which
is precisely what GRAF holds operationally and what the published table does not carry. That is
a concrete, answerable request to put to them, and it is a better one than the definitional
questions in §15 because it names the exact variable.

### 21.6 Standing conclusion

**`lcl_ratio_overshoot_max = 1.60` remains the default**, but the reason has changed and is
worth stating precisely. It is *not* that the empirical proxy outperforms the physics -- with
the PFT actually resolving, the explicit gate edges ahead. It is that **the available labels
cannot resolve the comparison**: the one metric that separates the variants is decided by two
sondes of one fire, and the physically preferable variant is penalised by a timescale mismatch
built into the label set.

That is a stronger position than the previous one. The proxy is retained on an argued
measurement limitation with a named remedy, not on an unexplained empirical win.

The corollary for the rest of this document stands, with the bar re-measured on the corrected
27-sonde sample: several sections propose mechanisms better *motivated* than what they would
replace. Motivation is not evidence. Each should clear **12/27 exact, +0.04 bias** -- and clear
it by more than two sondes of one fire -- before becoming a default.

## 22. Provenance

- GRAF map: `Tipus de piroconvecció - ICON-EU 26jul2026 00Z`, Bombers de la Generalitat de
  Catalunya, supplied as a WhatsApp image on 2026-07-26.
- pyflam product: `docs/forecast_2026-07-26_3day/`, generated by
  `scripts/forecast_3day_tuscany.py 2026-07-26 0`, effective source `hybrid`, ladder `shear`
  on all three days.
- Method: Castellnou et al. (2022), *JGR-Atmos*; dry-pyrocloud diagnostic after Castellnou
  Ribau et al. (2024).

## 20. Product reproducibility -- a hard operational limit

Found while regenerating after the physics corrections. **DWD retains ICON-EU open data for
about 24 hours.** On 2026-07-27 the 25 and 26 July 00Z runs both return HTTP 404; only the
current day is fetchable.

The consequence is that a hybrid product can only be regenerated **on the day it was made**.
After that its rasters are frozen, whatever is later found wrong with the physics that produced
them. Specifically:

| product | state |
|:--|:--|
| `docs/forecast_2026-07-25_3day*/` | predates every correction; source data gone; **cannot be regenerated** |
| `docs/forecast_2026-07-26_3day/` | carries the θ′ and nodata corrections, but predates the PFT wiring, so it has no `pft_gw` rasters; source data gone |
| `docs/forecast_2026-07-27_3day/` | first product with the complete corrected physics |

This is not a filing inconvenience. It means **published pyflam products carry the physics of
the day they were run and cannot be brought forward**, so any error found later is permanent in
the archive. Two implications worth acting on:

* Products should record the commit they were generated at, so a reader can tell which
  corrections they predate. Nothing in the current output does this.
* Where a product matters beyond the day, the *inputs* need archiving alongside it -- the ICON
  GRIB cache is ~8 GB per 3-day run, which is the real cost of reproducibility here.

Until then the safe reading is that anything under `docs/forecast_2026-07-2{5,6}_*` is a
historical artefact of superseded physics, not a current forecast.
