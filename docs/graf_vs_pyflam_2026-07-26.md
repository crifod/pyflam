# GRAF vs pyflam hybrid — pyroconvection type, run 2026-07-26 00Z

Comparison of the Catalan fire service (GRAF / Bombers) pyroconvection-type product against
the pyflam hybrid product over central Italy, plus the code change the comparison motivated.

**Status: for evaluation.** Written 2026-07-26; extended the same day with §§7-10 after a
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

## 11. Revised open items

| # | Item | Status |
|:--|:--|:--|
| 1 | Remove the unfounded 600 m ABL gate | **done, §9** |
| 2 | `bulk_richardson_abl_height` Ri_c 0.25 → 0.33; reconcile the 200/400 m start | **done, §10.1** |
| 3 | Implement **Δθ** (entrainment jump) | **done as diagnostic, §10.2** — gate blocked on item 8 |
| 4 | Implement **FireCAPE** (Potter 2005) | **done as diagnostic, §10.2** |
| 5 | Define the ML gradient on the residual layer | **done, §10.4** — evening classifiable 48 → 81 % |
| 6 | Re-run the GRAF comparison once the above are in | now unblocked |
| 7 | Confirm the GRAF time convention (§4) | open, needs the authors |
| 8 | **Calibrate the fire θ-excess** against the campaign's ten measured values; settle first whether the in-plume measurement and `fire_parcel_theta_excess` are the same quantity | open — **precondition for the penetration gate, §10.3** |
| 9 | Regenerate the published products: every figure and raster in `docs/forecast_2026-07-2{5,6}_3day*/` predates items 1-5 and the classes have since changed | open |

Items 2 and 5 were corrections; 3 and 4 added physics taken from the published method rather
than calibration against GRAF output, so model independence is preserved. Item 8 would
calibrate against the campaign's *published measurements*, not against the GRAF product. **No threshold in
this codebase has been altered to match the GRAF product**, and none should be.

## 11. Provenance

- GRAF map: `Tipus de piroconvecció - ICON-EU 26jul2026 00Z`, Bombers de la Generalitat de
  Catalunya, supplied as a WhatsApp image on 2026-07-26.
- pyflam product: `docs/forecast_2026-07-26_3day/`, generated by
  `scripts/forecast_3day_tuscany.py 2026-07-26 0`, effective source `hybrid`, ladder `shear`
  on all three days.
- Method: Castellnou et al. (2022), *JGR-Atmos*; dry-pyrocloud diagnostic after Castellnou
  Ribau et al. (2024).
