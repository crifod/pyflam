"""Italian/English strings for the pyroconvection products.

Both report generators -- ``tests/pyroconv_daily.py`` (one valid day) and
``scripts/forecast_3day_tuscany.py`` (three) -- emit a parallel Italian edition alongside the
English one, from the same rasters and the same run. The two editions must never disagree
about a number, so nothing here computes anything: every entry is prose or a label, and every
figure/table value is formatted by the caller and substituted in.

Conventions for the Italian, chosen to match how the terms are actually used in Italian fire
services and in the Catalan/Spanish literature this method comes from:

* **pyroCu / pyroCb / pyroconvection** stay as they are -- they are the international terms,
  and "cumulo pirogenico" would be a coinage nobody uses operationally.
* **entrainment**, **overshooting**, **spotting** stay in English for the same reason.
* *fuel gate* -> "filtro combustibile"; *fireline intensity* -> "intensita del fronte";
  *head fire* -> "fronte di testa"; *mixed layer* -> "strato rimescolato";
  *residual layer* -> "strato residuo"; *shear* -> "shear" (kept, as in Italian meteorology).
* The formal "si" impersonal is used rather than the imperative, matching the register of an
  operational bulletin.

Accented characters are written as plain ASCII where the surrounding pipeline is ASCII-only
(figure labels go through matplotlib's default font); the report prose is UTF-8 and uses
proper accents.
"""

LANGS = ("en", "it")

# Suffix appended to the .md/.pdf basename. English keeps the historical name so existing
# links and any downstream automation are untouched; Italian gets the tag.
SUFFIX = {"en": "", "it": "_it"}

# --- figure strings -------------------------------------------------------------------
# Kept ASCII-safe: these are drawn by matplotlib, which falls back to a substitute glyph for
# characters the default font lacks, and a missing accent in a printed panel is worse than
# an absent one.
FIG = {
    "en": {
        "class_legend_title": "Pyroconvection class (0 = lowest -> 4 = highest)",
        "nc": "n/c  not classifiable (no usable column)",
        "nc_short": "n/c  not classifiable",
        "potential_title": "POTENTIAL -- atmospheric upper bound (assumes a pyroCu-capable fire everywhere)",
        "gated_title": "FUEL-GATED -- expected (fire power >= 10 MW/m on the Tuscany fuels)",
        "potential_title_daily": "POTENTIAL -- atmosphere only (upper bound: assumes a pyroCu-capable fire in every cell)",
        "gated_title_daily": "FUEL-GATED -- expected (only where fire power >= 10 MW/m on the Tuscany .lcp fuels)",
        "type_prefix": "Pyroconvection type",
        "run": "run",
        "valid": "VALID",
        "tuscany": "Tuscany",
        "decoup_title": "Tuscany dry-pyrocloud decoupling -- 3-day forecast -- fireABL / ABL",
        "decoup_title_daily": "Dry-pyrocloud decoupling  fireABL / ABL",
        "decoup_ref": "reference fire: {k:.0f} K plume excess -- DIAGNOSTIC, no class",
        "decoup_cb": "fireABL / ABL   (1 = no decoupling; higher = deeper dry decoupling)",
        "ptop_title_daily": "Predicted plume-top height  --  the cost ladder inverted",
        "ptop_ref": "declared reference fire {fp:.0f} GW -- a scenario, not a per-cell estimate; validated on 1231 MISR plumes, rho +0.46",
        "ptop_cb": "plume top (m above ground)",
        "margin_title": "Tuscany pyroCb firepower margin -- 3-day forecast -- firepower / PFT",
        "margin_title_daily": "PyroCb firepower margin  firepower / PFT",
        "margin_bridge": "head fire {m:.0f} m, convective fraction {c:g}",
        "margin_cb": "firepower / PFT",
        "margin_hit": "margin >= 1 AND fuel gate passed: pyroCb candidate",
        "margin_miss": "margin >= 1 but gate failed: degenerate low-PFT column",
        "nofuel": "no burnable fuel (firepower = 0)",
        "src_hybrid": "ICON-EU model levels + ICON-2I 2.2 km gate",
        "src_icon2i": "ICON-2I 2.2 km",
    },
    "it": {
        "class_legend_title": "Classe di piroconvezione (0 = minima -> 4 = massima)",
        "nc": "n/c  non classificabile (nessuna colonna utilizzabile)",
        "nc_short": "n/c  non classificabile",
        "potential_title": "POTENZIALE -- limite superiore atmosferico (presuppone ovunque un incendio capace di pyroCu)",
        "gated_title": "FILTRATO PER COMBUSTIBILE -- atteso (potenza >= 10 MW/m sui combustibili toscani)",
        "potential_title_daily": "POTENZIALE -- sola atmosfera (limite superiore: presuppone in ogni cella un incendio capace di pyroCu)",
        "gated_title_daily": "FILTRATO PER COMBUSTIBILE -- atteso (solo dove la potenza raggiunge 10 MW/m sui combustibili .lcp toscani)",
        "type_prefix": "Tipo di piroconvezione",
        "run": "corsa",
        "valid": "VALIDO",
        "tuscany": "Toscana",
        "decoup_title": "Disaccoppiamento pirogeno secco in Toscana -- previsione 3 giorni -- fireABL / ABL",
        "decoup_title_daily": "Disaccoppiamento pirogeno secco  fireABL / ABL",
        "decoup_ref": "incendio di riferimento: eccesso di {k:.0f} K nel pennacchio -- DIAGNOSTICO, nessuna classe",
        "decoup_cb": "fireABL / ABL   (1 = nessun disaccoppiamento; valori maggiori = disaccoppiamento secco piu profondo)",
        "ptop_title_daily": "Altezza di cima del pennacchio prevista  --  la scala dei costi invertita",
        "ptop_ref": "incendio di riferimento dichiarato {fp:.0f} GW -- uno scenario, non una stima per cella; validata su 1231 pennacchi MISR, rho +0,46",
        "ptop_cb": "cima del pennacchio (m dal suolo)",
        "margin_title": "Margine di potenza per pyroCb in Toscana -- previsione 3 giorni -- potenza / PFT",
        "margin_title_daily": "Margine di potenza per pyroCb  potenza / PFT",
        "margin_bridge": "fronte di testa {m:.0f} m, frazione convettiva {c:g}",
        "margin_cb": "potenza / PFT",
        "margin_hit": "margine >= 1 E filtro combustibile superato: candidato pyroCb",
        "margin_miss": "margine >= 1 ma filtro non superato: colonna degenere a PFT bassa",
        "nofuel": "nessun combustibile bruciabile (potenza = 0)",
        "src_hybrid": "livelli modello ICON-EU + filtro ICON-2I 2.2 km",
        "src_icon2i": "ICON-2I 2.2 km",
    },
}

# Class labels for the legend. Taken here rather than from
# ``pyflam.atmosphere.PYROCONVECTION_TYPE_LABEL`` so the Italian edition does not depend on
# changing a library constant that regression tests assert on.
CLASS_LABEL = {
    "en": {"surface_plume": "Surface plume", "convection_plume": "Convection plume",
           "overshooting_pyrocu": "Overshooting pyroCu", "resilient_pyrocu": "Resilient pyroCu",
           "deep_pyrocu_pyrocb": "Deep pyroCu / pyroCb"},
    "it": {"surface_plume": "Pennacchio superficiale", "convection_plume": "Pennacchio convettivo",
           "overshooting_pyrocu": "PyroCu overshooting", "resilient_pyrocu": "PyroCu persistente",
           "deep_pyrocu_pyrocb": "PyroCu profondo / pyroCb"},
}

# --- report prose ---------------------------------------------------------------------
# Every entry is a format string; the caller supplies the numbers. Keys are shared between
# the two generators where the text is shared, and suffixed ``_3day``/``_daily`` where the
# two products genuinely say different things.
MD = {
"en": {
"title_3day": "Tuscany Pyroconvection -- 3-Day Forecast -- run {rundate} {run:02d}Z",
"subtitle_3day": "Valid {d0} .. {d2}. Hybrid: ICON-EU model levels + ICON-2I 2.2 km fuel gate. Method after Castellnou et al. (2022).",
"title_daily": "Tuscany Pyroconvection-Type Forecast -- {src} -- VALID {date}",
"subtitle_daily": "3-hourly, 24 h. Run {rundate} {run:02d}Z. Method after Castellnou et al. (2022), JGR-Atmos.",

"h_questions": "The three questions this report answers",
"questions_intro": """Each map answers one question. Read them in order -- but see the warning below question 3:
they are **not** nested, and a cell can pass question 3 while failing question 2.""",
"questions_table": """| # | Question | Map | The test |
|:--|:--|:--|:--|
| 1 | What is the most intense pyroconvection this **column** could sustain? | POTENTIAL | the ladder, with no fire-side condition at all |
| 2 | Is there enough fire here to make a **plume** at all? | FUEL-GATED | Byram intensity >= 10 MW/m on the .lcp fuels |
| 3 | Is there enough fire here to make a **pyroCb** in this specific column? | PFT MARGIN | total firepower >= that column's own PFT |""",
"questions_note": """Question 1 is a property of the atmosphere alone -- no fuel or land information enters it,
not even in deciding where it is defined. Questions 2 and 3 both condition on fire, and both
descend from the same Byram intensity field, but they are **not** the same test and they do
not agree: (2) compares power *per metre of front* against a fixed 10 MW/m constant, while
(3) converts that to a *total* power through an assumed head-fire length and compares it
against a threshold computed for that column. Cells routinely clear (2) and fail (3) by
orders of magnitude. Question 3 is the operationally decisive one and the hardest to satisfy.""",
"questions_fourth": """A fourth map, **dry-pyrocloud decoupling**, supports all three: it is the dry counterpart to
the moist class ladder and carries no class label. Read class *counts* as indicative.""",

"h1": "1. Potential -- what could this column sustain? (atmospheric upper bound)",
"h2": "2. Fuel-gated -- is there enough fire for a plume at all?",
"h3": "3. PFT margin -- is there enough fire for a pyroCb in *this* column?",
"h4": "{n}. Dry-pyrocloud decoupling -- DIAGNOSTIC (no class label)",
"h_ptop": "{n}. Predicted plume-top height -- the only field validated against observation",
"cap_ptop": "PREDICTED PLUME TOP. The cost ladder solved for height instead of firepower, for a declared {fp:.0f} GW fire. Continuous, no class label.",
"p_ptop": """The three questions above ask *how much fire* a given height costs. Solving the same
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

**Read the firepower as a scenario.** The map is drawn for a declared {fp:.0f} GW *total*
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
include columns the classifier rejected. The other `diag_*.tif` stay raw, as they always have.""",
"h_margin_table": "PyroCb firepower margin -- daytime, % of the burnable classifiable land",

"cap1": "Question 1 -- POTENTIAL. Pyroconvection class from the atmosphere alone, no fire-side condition. Upper bound, not an expectation.",
"cap2": "Question 2 -- FUEL-GATED. The same ladder with every cell below 10 MW/m of Byram intensity on the Tuscany .lcp fuels forced to class 0. The expected product.",
"cap3": "Question 3 -- PFT MARGIN. Firepower over each column's own PyroCb Firepower Threshold, log scale pivoting at the criterion. Rings meet the criterion and pass the fuel gate; crosses meet it only because the PFT has degenerated.",
"cap4": "Supporting diagnostic -- DRY-PYROCLOUD DECOUPLING. fireABL/ABL for a reference intense fire. Continuous, no class label.",

"p1": """The ladder with the fire-side condition switched off entirely, so this is the atmosphere's
answer and nothing else: the most intense pyroconvection the column would support **given a
fire able to exploit it**. That fire is assumed, not computed -- it enters through the ladder's
calibrated thresholds, not through any fuel map. Use this map to see where the atmosphere is
the binding constraint, never as an expectation. Note that the classifier runs on any usable
column, including over water; the sea is masked in the figure but the underlying
`pyroconv_potential_*.tif` rasters are not, so zonal statistics taken straight off them must
be restricted to land.""",
"p2": """The same ladder, on the same columns, with each cell forced to class 0 where the Tuscany .lcp
fuels and the forecast surface weather do not support **10 MW/m** of Byram fireline intensity
(Tedim et al. 2018). The gate is binary and one-directional: above the threshold the fuels
have no further influence, below it the cell reads surface plume whatever the atmosphere says.
This map therefore cannot exceed the potential map anywhere -- every difference between the
two is a cell knocked to 0 -- and it is the expected, operational product.""",
"p3a": """Per-cell **total** firepower divided by that column's own **PyroCb Firepower Threshold**
(Tory & Kepert 2021, eq. 31: `PFT = 0.3 z_fc^2 U_ML dtheta_fc`, in GW). At or above 1 the fire
is powerful enough to force a pyroCb *in that column*; below 1 it is not, however favourable
the atmosphere. The colour pivots exactly at 1 and the scale is logarithmic, because the field
spans four decades.""",
"p3b": """Byram intensity is power per metre of front and the PFT is a total power, so bridging them
needs an active head-fire length. Rather than assume one it is taken from the fire behaviour
the Wildfire Data Portal publishes -- the median over its 20 fires with both burn ratio and
ROS, {m:.0f} m -- with a convective fraction of {c:g} of the released heat entering the plume
(Tory & Kepert app. D). Both are assumptions, and the margin scales linearly with each: a
5 km head fire, which Tory & Kepert give for an extreme case, would multiply every margin here
by about seven. **Read the margin as an order of magnitude, not a number.**""",
"p3c": """This is why the map matters despite that caveat: over most of the burnable domain the margin
is not near 1 but two to three decades below it, so the conclusion "no pyroCb here" is robust
to the head-fire assumption in a way that a marginal cell would not be. Where the margin does
approach 1, the assumption becomes load-bearing and the cell deserves judgement rather than
the map.""",
"p3d": """Three exclusions are drawn differently on purpose: **grey** is a column the classifier
rejected (no atmosphere to test against), **white** is land with no burnable fuel (firepower
exactly zero -- a categorical no, not a small margin), and the sea is masked. The underlying
field is `diag_pft_margin_*.tif`, with its numerator and denominator alongside as
`diag_firepower_gw_*.tif` and `diag_pft_gw_*.tif`.""",

"h_nested": "Questions 2 and 3 are not nested -- read the markers, not the colour",
"p_nested": """A cell can clear the PFT criterion while **failing** the 10 MW/m fuel gate, and on this run
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
the raw margin count instead would overstate the candidate area several-fold.""",
# --- daily-report-only blocks ---------------------------------------------------------
"sup_txt": """Over the classified land the fit is supported in **{pct:.0f}%** of columns at peak of day (i.e. that share holds at least {pts} levels inside the layer; the rest return a nan gradient and leave the class map, which is also what thins the evening panels).""",
"sup_none": "The share of columns supporting the fit was not recorded for this run.",
"ml_para_hybrid": """The **mixed-layer stability** here is a genuine measurement, not a proxy. The ICON-EU
native model levels put ~10 levels inside the mixed layer (this run: {n} median over land at
12-15Z, when the layer is mature -- the night and evening counts are far lower, but nothing is
classified then). {sup} So the
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
weight those hours accordingly.""",
"ml_para_icon2i": """The **mixed-layer stability** diagnostic is the weakest link in this product, and is a
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
measurements.""",
"ml_rows_hybrid": """| ML dtheta/dz (least-squares fit, in mixed layer) | > 1.1e-3 K/m (stable) | Convection plume only -- no pyroCu |
| ML dtheta/dz (measured on model levels) | <= 1.1e-3 K/m | Column is pyroCu-capable (bias ~-0.4e-4 K/m vs radiosondes) |""",
"ml_rows_icon2i": """| ML stability *proxy*: parcel mixing depth (see Method) | depth < 455 m ("stable") | Convection plume only -- no pyroCu |
| ML stability proxy | depth >= 455 m (= dtheta/dz <= 1.1e-3 K/m) | Column is pyroCu-capable (weak filter: spec. 0.29 inland) |""",
"forcing_hybrid": """Forcing: atmosphere from ICON-EU 6.5 km native model levels (DWD open data, CC-BY), lowest ~24 levels; ~{n} inside the mixed layer here (land, 12-15Z), regridded to the ICON-2I 2.2 km grid. Surface fields and fuel gate from ICON-2I 2.2 km (MISTRAL / AgenziaItaliaMeteo). Classifier: pyflam.pyroconvection_type (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta, ABL-top RH; no surface CAPE).""",
"forcing_icon2i": """Forcing: ICON-2I 2.2 km full-Italy GRIB (MISTRAL / AgenziaItaliaMeteo, CC-BY), Tuscany subset: geopotential, temperature, RH and wind on {n} pressure levels plus the 2 m / 10 m state, surface pressure and orography. Classifier: pyflam.pyroconvection_type (bulk-Richardson ABL, Bolton LCL, mixed-layer dtheta/dz proxy, cap gamma-theta, ABL-top RH; no surface CAPE).""",
"shear_partial": """**Ladder actually used for this run: `{lad}`.** The fifth diagnostic -- the distance from the ABL/LCL to the height of maximum wind shear -- resolved over part of the domain only, so the full 5-diagnostic ladder ran there and the reduced one elsewhere. That clause is a *necessary* condition for the top class, so class 4 is somewhat **easier** to reach in the cells that fell back than in the ones that did not. Treat class 4 as an alert to inspect the column, not as a calibrated probability.""",
"shear_full": """**Ladder actually used for this run: `{lad}`.** The full 5-diagnostic method ran: the model levels resolved a shear-maximum height, so class 4 additionally required that maximum to sit within 0.30 ABL of the ABL/LCL -- a *necessary* condition the reduced ladders omit. Class 4 here therefore carries its shear clause; it remains an alert to inspect the column rather than a calibrated probability.""",
"shear_none": """**Ladder actually used for this run: `{lad}`.** The full method has a fifth diagnostic -- the distance from the ABL/LCL to the height of maximum wind shear -- which this source cannot resolve, so the reduced ladder runs. The shear clause is a *necessary* condition for the top class, so omitting it makes class 4 somewhat **easier** to reach here than in the full method. Treat class 4 as an alert to inspect the column, not as a calibrated probability.""",
"shear_row_partial": "| Shear-maximum distance / ABL | <= 0.30 | Required for class 4 where the shear height resolved -- **applied over part of the domain** (see Method) |",
"shear_row_full": "| Shear-maximum distance / ABL | <= 0.30 | Required for class 4 -- **resolved and applied throughout this run** |",
"shear_row_none": "| Shear-maximum distance / ABL | <= 0.30 | Required for class 4 **in the 5-diagnostic ladder only** (not resolvable here -- see Method) |",
"ladder_none": "none (no classifiable cell)",
"h_howto": "How to read this product",
"h_method": "Method and its limits (read this before using the classes)",
"p_method": """The column is classified from a **real vertical profile**, not a standard atmosphere:
{heights}, the mixing depth from the **bulk Richardson number** (first Rib >= 0.33
above 200 m AGL, referenced to the 2 m / 10 m state), the LCL from the exact
**Bolton (1980)** formula, and the humidity at the ABL top from the model's RH. The cap
gamma-theta is taken over ABL+200 m to ABL+1200 m.""",
"heights_hybrid": "per-cell heights from the ICON-EU model-level heights (HHL)",
"heights_icon2i": "per-cell heights from the model geopotential (FI)",
"p_predisposition": """The classes express atmospheric predisposition **given a fire of sufficient power**;
in the gated panel that power is computed, not assumed. They are paper-informed
thresholds, not locally validated ones.""",
"h_classscale": "Class scale (low -> high pyroconvective activity)",
"classscale": """| Level | Colour | Class | Meaning |
|:--:|:--|:--|:--|
| 0 | white | Surface plume | Buoyant smoke plume; no significant cloud development. |
| 1 | green | Convection plume | Plume penetrates a stable mixed layer; condensation possible, no pyroCu. |
| 2 | yellow | Overshooting pyroCu | Brief pyrocumulus; cloud base above the mixing height (LCL/ABL > 1). |
| 3 | orange | Resilient pyroCu | Persistent pyrocumulus in an unstable column (LCL/ABL < 1). |
| 4 | dark red | Deep pyroCu / pyroCb | Deep pyroconvection / pyrocumulonimbus; weak upper cap lets the plume deepen. |""",
"h_thresholds": "Classification thresholds (ladder in use: `{lad}`)",
"thr_head": """| Diagnostic | Threshold | Effect |
|:--|:--|:--|""",
"thr_rest": """| LCL / ABL ratio | 1.0 -- 1.60 | Overshooting pyroCu (brief) |
| LCL / ABL ratio | < 1.0 | Resilient pyroCu (persistent) |
| LCL / ABL ratio | <= 1.10 (+ conditions below) | Admissible for deep pyroCu / pyroCb |
| Cap gamma-theta (ABL+200 m -> ABL+1200 m) | <= 4.2e-3 K/m (weak cap) | Permits deepening to pyroCb |
| Cap gamma-theta | >= 4.8e-3 K/m (strong cap) | Inhibits deepening (resilient at most) |
| RH at the ABL top (mean, ABL +/- 150 m) | >= 60% | Required for classes 3 and 4 (was 80%; relaxed for dry fire weather) |""",
"thr_tail": """| Fireline intensity (fuel gate, gated panel) | >= 10 MW/m | Minimum fire power for any pyroCu (Tedim et al. 2018) |
| ABL depth | < {abl:d} m | Not classifiable (implausible depth). No 600 m gate: it had no basis in Castellnou et al. (2022) and discarded 5 of the 8 campaign fires -- removed 2026-07-26 |
| Usable pressure levels | < 4 | Cell not classified |""",
"h_refcases": "Reference cases (Castellnou et al. 2022, Table 1)",
"p_refcases": """The paper's labelled events, which anchor the published 3-diagnostic ladder
(`pyflam.pyroconvection_type(ladder="castellnou")`, still the library default and
regression-tested against these rows). The profile ladders used for this map are
stricter at the top: they additionally require a moist ABL top and an LCL close to it.""",
"refcases": """| Case | Observed type | LCL/ABL | ML dtheta/dz | gamma-theta (700-500) |
|:--|:--|:--:|:--|:--:|
| T21 | Convection plume | -- | stable | -- |
| SCQ32 | Overshooting pyroCu | > 1 | neutral/unstable | -- |
| M11 | Resilient pyroCu | < 1 | unstable | 4.2e-3 (resilient cap) |
| SCQ41 | pyroCu, not pyroCb | < 1 | unstable | 5.1e-3 (strong cap) |
| SCQ51 | Deep pyroCu / pyroCb | < 1 | unstable | 3.9e-3 (weak cap) |""",
"p_footer_daily": """Gate: Rothermel + Cruz-2005 crown on the .lcp fuels with forecast moisture/wind.
Per-hour diagnostic rasters accompany the classes: ABL, parcel_ml, residual_ml, LCL, LCL/ABL,
ML dtheta/dz, cap gamma-theta, RH-top, fireABL, decoupling, and (hybrid only) delta_theta,
firecape, penetration, pft_gw (PyroCb Firepower Threshold, GW), z_fc, delta_theta_fc, u_ml,
firepower_gw, pft_margin.
Province borders: ISTAT-derived (openpolis geojson-italy).
Generated by tests/pyroconv_daily.py at pyflam `{ver}`, {at} -- frozen with that
commit's physics, since DWD drops the run after ~24 h.""",
"h_intro_3day": "3-day forecast",
"p_intro_3day": """Three valid days from the single {rundate} {run:02d}Z run (day 0, +1, +2), using the day+1/+2
forecast steps in the same ICON-EU + ICON-2I datasets. Each day is a **complete 24 h cycle,
3-hourly, starting at 00Z of the run day** ({hours}), so the three day-rows of every figure
are directly comparable and day 0 opens at the run's own initial time. Panels are the Tuscany
domain with ISTAT province borders; the sea is masked. Atmosphere from ICON-EU native model
levels (bulk-Richardson ABL, Bolton LCL, measured mixed-layer dtheta/dz, cap gamma-theta,
ABL-top RH, shear); the 10 MW/m fuel gate uses ICON-2I's 2.2 km surface fields on the Tuscany
.lcp fuels. See scripts/validation/README.md for the diagnostic validation (radiosondes + ERA5).""",
"p_decoup": """The **decoupling ratio** fireABL / ABL is how high a reference intense fire
(a {k:.0f} K plume temperature excess -- the intense end of the GRAF in-plume measurements,
0.1-13.1 K) would grow its own boundary layer by *sensible heat alone*, divided by the ambient
ABL. It is the **dry** counterpart to the moist class maps above (Castellnou et al. 2022;
Castellnou Ribau et al. 2024): values well above 1 mark deep, hot, dry columns where a fire can
punch through and decouple from the surface *even where the moist ladder scores low*, which is
exactly the situation the gated map under-reports. Like the potential map it assumes a fire
everywhere -- here a fixed reference flux rather than a computed one -- so it is an upper bound,
not an expectation, and no dry/moist LCL split is applied (the +1 km literature offset is not
supported by the GRAF prototype labels). fireABL from
`pyflam.atmosphere.fire_induced_abl_grid`, the parcel intersected with the real theta(z) stack;
cells the classifier rejects (no usable column, or ABL below {abl:d} m) are left blank rather
than painted with a bare quotient.""",
"h_dtbl": "Dry decoupling -- daytime, % of the classifiable land",
"dtbl_head": "| Day | Hour | % grid classifiable | % area ratio >= 2 | % area ratio >= 3 | max ratio |",
"p_dtbl": """*% grid classifiable* is the share of the domain with a usable column at that hour (ABL above
{abl:d} m); the two ratio columns are percentages **of that share**, not of Tuscany.
The denominator collapses towards sunset as the mixed layer decays, so an 18Z row can read
100% off a small residual area -- read it together with the classifiable column.""",
"h_ctbl": "Fuel-gated class distribution -- daytime, % of land cells ({n} land cells)",
"ctbl_head": "| Day | Hour | surface | convect | overshoot | resilient | deep | n/c |",
"p_ctbl": """Classes: 0 surface plume, 1 convection plume, 2 overshooting pyroCu, 3 resilient pyroCu,
4 deep pyroCu/pyroCb. **n/c** is land with no classifiable column (no usable profile, or an
ABL below {abl:d} m) -- distinct from class 0, which is the ladder's finding that the
column supports only a surface plume. Around sunset n/c takes most of the domain; those hours
carry no forecast, and must not be read as quiet ones.""",
"h_evening": "Evening columns and the residual layer",
"p_evening": """After sunset the convective layer collapses into a shallow stable layer, but the air above it
keeps the day's near-neutral profile -- the **residual layer**. A surface-referenced parcel
finds only the stable layer, and a mixed-layer dtheta/dz fitted inside that depth has no levels
to work with, so the whole domain used to fall out as unclassifiable. The gradient is now
fitted over the deeper of the parcel and residual depths (`atmosphere.residual_layer_grid`),
which by day reduces exactly to the parcel depth and only differs in the evening. **This changes
the evening classes**: those hours now carry a forecast where they previously carried none.
Treat them with more caution than the midday ones -- the diagnostic is newer and the residual
depth it returns still looks shallow against what the profile suggests.

Two further diagnostics ride along and are **not** used by the classifier: the entrainment-zone
jump `delta_theta` (the potential-temperature step a plume must cross to escape the mixed layer)
and `firecape` (Potter 2005, the CAPE of a fire-heated parcel). Both come from the source
method's variable set; neither gates a class here, because the penetration test they imply does
not yet discriminate at this product's reference fire strength. See
`docs/graf_vs_pyflam_2026-07-26.md`.""",
"p_footer": """Per-hour class and diagnostic GeoTIFFs accompany this report in `rasters_hybrid_<day>/`:
ABL, parcel_ml, residual_ml, LCL, LCL/ABL, ML dtheta/dz, cap gamma-theta, RH-top, fireABL,
decoupling, delta_theta, firecape, penetration, plus the fire side -- fuel_load,
burnable_fraction, firepower_gw, pft_gw, pft_margin, crit_growth_ha_h, z_fc, u_ml,
delta_theta_fc.
Generated by scripts/forecast_3day_tuscany.py at pyflam `{ver}`, {at}.
This product cannot be regenerated once DWD drops the run (~24 h), so it is frozen with
the physics of that commit.""",
"mtbl_head": "| Day | Hour | % grid with firepower | % margin >= 1 | % margin >= 0.1 | max margin | cells >= 1 | **also gate-passing** |",
"mtbl_head_daily": "| Hour | % grid with firepower | % margin >= 1 | % margin >= 0.1 | max margin | cells >= 1 | **also gate-passing** |",
"p_mtbl": """*% grid with firepower* is the share of the domain that has both a classifiable column and
burnable fuel; the two margin columns are percentages **of that share**, not of Tuscany. A
percentage of the whole domain would be dominated by cells that carry no fuel and can never
contribute, and would move for reasons unrelated to the forecast. *cells >= 1* is the raw
count behind the first percentage, carried because at these rates a percentage rounds to
something that looks like nothing when it is in fact a handful of specific places.
***also gate-passing*** applies the joint criterion of the previous section and is the count
to act on; the gap between it and *cells >= 1* is the degenerate-PFT population.""",
},

"it": {
"title_3day": "Piroconvezione in Toscana -- previsione a 3 giorni -- corsa {rundate} {run:02d}Z",
"subtitle_3day": "Valida dal {d0} al {d2}. Ibrido: livelli modello ICON-EU + filtro combustibile ICON-2I 2.2 km. Metodo secondo Castellnou et al. (2022).",
"title_daily": "Previsione del tipo di piroconvezione in Toscana -- {src} -- VALIDA {date}",
"subtitle_daily": "Trioraria, 24 h. Corsa {rundate} {run:02d}Z. Metodo secondo Castellnou et al. (2022), JGR-Atmos.",

"h_questions": "Le tre domande a cui questo bollettino risponde",
"questions_intro": """Ogni mappa risponde a una domanda. Vanno lette nell'ordine -- ma si veda l'avvertenza sotto la
domanda 3: le domande **non** sono annidate, e una cella puo superare la 3 fallendo la 2.""",
"questions_table": """| # | Domanda | Mappa | Il criterio |
|:--|:--|:--|:--|
| 1 | Qual è la piroconvezione più intensa che questa **colonna** potrebbe sostenere? | POTENZIALE | la scala diagnostica, senza alcuna condizione sul fuoco |
| 2 | C'è abbastanza fuoco per generare un **pennacchio**? | FILTRATO PER COMBUSTIBILE | intensità di Byram >= 10 MW/m sui combustibili .lcp |
| 3 | C'è abbastanza fuoco per generare un **pyroCb** in questa specifica colonna? | MARGINE PFT | potenza totale >= PFT propria di quella colonna |""",
"questions_note": """La domanda 1 è una proprietà della sola atmosfera: nessuna informazione sui combustibili o
sul suolo vi entra, nemmeno per stabilire dove sia definita. Le domande 2 e 3 condizionano
entrambe al fuoco e discendono entrambe dallo stesso campo di intensità di Byram, ma **non**
sono lo stesso criterio e non concordano: la (2) confronta la potenza *per metro di fronte*
con una costante fissa di 10 MW/m, mentre la (3) la converte in potenza *totale* attraverso una
lunghezza di fronte di testa assunta e la confronta con una soglia calcolata per quella colonna.
Le celle superano regolarmente la (2) e falliscono la (3) di ordini di grandezza. La domanda 3
è quella operativamente decisiva, ed è la più difficile da soddisfare.""",
"questions_fourth": """Una quarta mappa, il **disaccoppiamento pirogeno secco**, supporta tutte e tre: è la
controparte secca della scala umida delle classi e non porta alcuna etichetta di classe. I
*conteggi* delle classi vanno letti come indicativi.""",

"h1": "1. Potenziale -- che cosa potrebbe sostenere questa colonna? (limite superiore atmosferico)",
"h2": "2. Filtrato per combustibile -- c'è abbastanza fuoco per un pennacchio?",
"h3": "3. Margine PFT -- c'è abbastanza fuoco per un pyroCb in *questa* colonna?",
"h4": "{n}. Disaccoppiamento pirogeno secco -- DIAGNOSTICO (nessuna etichetta di classe)",
"h_ptop": "{n}. Altezza di cima prevista -- l'unico campo validato contro osservazioni",
"cap_ptop": "CIMA DEL PENNACCHIO PREVISTA. La scala dei costi risolta per l'altezza invece che per la potenza, per un incendio dichiarato di {fp:.0f} GW. Continua, nessuna etichetta di classe.",
"p_ptop": """Le tre domande precedenti chiedono *quanto fuoco* costi una data quota. Risolvendo la
stessa disequazione nell'altro verso -- per la quota massima che un dato incendio puo
permettersi -- si risponde alla domanda che un previsore ha davvero, e si ottiene un campo
continuo invece di una classe.

E anche l'unico prodotto qui **verificato contro osservazioni da un capo all'altro**.
Confrontare la cima di un pennacchio con la quota-bersaglio di un gradino mette quel bersaglio
su entrambi i lati del test e misura per lo piu la propria circolarita; risolvere per l'altezza
mette un numero previsto contro un numero misurato e nient'altro. Valutato cosi su 1231
pennacchi digitalizzati dalle immagini stereo MISR su 7 regioni e 11 biomi, senza nulla di
adattato: Spearman **+0,46** rispetto alla cima osservata, **+125 m** di bias mediano, 576 m di
errore assoluto medio e 12,8 % meglio del predire una costante su un errore adimensionale --
davanti alla sola potenza (+0,37) e alla sola profondita dello strato limite (+0,34), e con le
mediane regionali ordinate a rho +0,61.

**La potenza va letta come scenario.** La mappa e disegnata per una potenza *totale* dichiarata
di {fp:.0f} GW, esplicitata qui perche una potenza totale non si ricava da un'intensita di Byram
senza assumere una lunghezza del fronte di testa, ed e proprio quell'assunzione che questo
filone di lavoro esiste per evitare. Raddoppiare l'incendio di riferimento non raddoppia
l'altezza -- il costo cresce con il cubo della salita -- ma il campo si sposta, e non va letto
come stima per cella di cio che brucera davvero.

**Due limiti da dichiarare.** Il campione di validazione e quello di MISR, quindi e fissato
intorno alle 10:30 solari locali: lo strato limite sta ancora crescendo e gli incendi sono piu
piccoli di quanto saranno nel pomeriggio che questo prodotto prevede. E il campo viene prodotto
solo dove la griglia verticale risolve lo strato su cui si legge la cappa, che in questa corsa
significa la sorgente a livelli modello; un ripiego su livelli di pressione lo lascia vuoto
invece di interpolarlo.

A differenza dei raster di classe, `diag_plume_top_*.tif` e `diag_form_*.tif` sono scritti
**gia mascherati** con il campo `valid` del classificatore: una statistica zonale presa
direttamente su di essi non puo includere colonne che il classificatore ha rifiutato. Gli altri
`diag_*.tif` restano grezzi, come sono sempre stati.""",
"h_margin_table": "Margine di potenza per pyroCb -- ore diurne, % del suolo bruciabile e classificabile",

"cap1": "Domanda 1 -- POTENZIALE. Classe di piroconvezione dalla sola atmosfera, nessuna condizione sul fuoco. Limite superiore, non un'attesa.",
"cap2": "Domanda 2 -- FILTRATO PER COMBUSTIBILE. La stessa scala diagnostica, con ogni cella sotto i 10 MW/m di intensità di Byram sui combustibili .lcp toscani forzata alla classe 0. È il prodotto atteso.",
"cap3": "Domanda 3 -- MARGINE PFT. Potenza dell'incendio rapportata alla PyroCb Firepower Threshold propria di ciascuna colonna, scala logaritmica imperniata sul criterio. I cerchi soddisfano il criterio e superano il filtro combustibile; le croci lo soddisfano solo perché la PFT è degenerata.",
"cap4": "Diagnostico di supporto -- DISACCOPPIAMENTO PIROGENO SECCO. fireABL/ABL per un incendio intenso di riferimento. Continuo, nessuna etichetta di classe.",

"p1": """La scala diagnostica con la condizione sul fuoco completamente disattivata: questa è dunque la
risposta della sola atmosfera. È la piroconvezione più intensa che la colonna sosterrebbe **dato
un incendio in grado di sfruttarla**. Quell'incendio è presupposto, non calcolato: entra
attraverso le soglie calibrate della scala, non attraverso una mappa di combustibili. Questa
mappa serve a vedere dove sia l'atmosfera il vincolo stringente, mai come previsione. Si noti
che il classificatore gira su qualunque colonna utilizzabile, anche sul mare: nella figura il
mare è mascherato, ma i raster `pyroconv_potential_*.tif` sottostanti no, quindi ogni statistica
zonale ricavata direttamente da essi va limitata alla terraferma.""",
"p2": """La stessa scala diagnostica, sulle stesse colonne, con ogni cella forzata alla classe 0 dove i
combustibili .lcp toscani e la meteorologia superficiale prevista non sostengono **10 MW/m** di
intensità del fronte secondo Byram (Tedim et al. 2018). Il filtro è binario e unidirezionale:
sopra la soglia i combustibili non hanno più alcuna influenza, sotto la soglia la cella riporta
un pennacchio superficiale qualunque cosa dica l'atmosfera. Questa mappa non può quindi mai
superare la mappa del potenziale -- ogni differenza fra le due è una cella riportata a 0 -- ed è
il prodotto operativo atteso.""",
"p3a": """Potenza **totale** dell'incendio per cella, divisa per la **PyroCb Firepower Threshold** propria
di quella colonna (Tory & Kepert 2021, eq. 31: `PFT = 0.3 z_fc^2 U_ML dtheta_fc`, in GW). Da 1
in su l'incendio è abbastanza potente da forzare un pyroCb *in quella colonna*; sotto 1 non lo
è, per quanto favorevole sia l'atmosfera. Il colore si imperna esattamente su 1 e la scala è
logaritmica, perché il campo si estende su quattro decadi.""",
"p3b": """L'intensità di Byram è una potenza per metro di fronte e la PFT è una potenza totale: per
collegarle serve una lunghezza attiva del fronte di testa. Anziché assumerne una, la si ricava
dal comportamento del fuoco pubblicato dal Wildfire Data Portal -- la mediana sui suoi 20
incendi che riportano sia il rapporto di superficie bruciata sia la velocità di avanzamento,
{m:.0f} m -- con una frazione convettiva pari a {c:g} del calore rilasciato che entra nel
pennacchio (Tory & Kepert app. D). Sono entrambe assunzioni, e il margine scala linearmente con
ciascuna: un fronte di testa di 5 km, che Tory & Kepert indicano per un caso estremo,
moltiplicherebbe ogni margine qui riportato per circa sette. **Il margine va letto come un
ordine di grandezza, non come un numero.**""",
"p3c": """È questa la ragione per cui la mappa conserva valore nonostante l'avvertenza: sulla maggior
parte del dominio bruciabile il margine non è prossimo a 1, ma due o tre decadi al di sotto,
sicché la conclusione «qui nessun pyroCb» è robusta rispetto all'assunzione sul fronte di testa
in un modo in cui una cella marginale non lo sarebbe. Dove invece il margine si avvicina a 1,
l'assunzione diventa portante e la cella merita un giudizio, non la mappa.""",
"p3d": """Tre esclusioni sono disegnate in modo diverso, deliberatamente: il **grigio** è una colonna che
il classificatore ha rifiutato (nessuna atmosfera su cui applicare il criterio), il **bianco** è
suolo privo di combustibile bruciabile (potenza esattamente nulla: un no categorico, non un
margine piccolo), e il mare è mascherato. Il campo sottostante è `diag_pft_margin_*.tif`, con
numeratore e denominatore a fianco come `diag_firepower_gw_*.tif` e `diag_pft_gw_*.tif`.""",

"h_nested": "Le domande 2 e 3 non sono annidate -- si leggano i simboli, non il colore",
"p_nested": """Una cella può soddisfare il criterio PFT **fallendo** il filtro dei 10 MW/m, e in questa corsa
è quanto accade per la maggior parte di esse. Il meccanismo sta al denominatore: la PFT è
`0.3 z_fc^2 U_ML dtheta_fc`, sicché in una colonna bassa, debolmente tappata e con vento
leggero essa collassa -- a fronte di una mediana di dominio prossima a 90 GW, le celle con
margine >= 1 hanno PFT di pochi GW. Un incendio da 2-8 MW/m «supera» allora una soglia che un
incendio di quelle dimensioni non può fisicamente dirsi aver battuto, perché a quell'intensità
non c'è alcun pyroCu da approfondire.

La mappa distingue i due casi, e solo il primo è una previsione:

* **cerchio nero** -- margine >= 1 *e* filtro combustibile superato. Un autentico candidato pyroCb.
* **croce grigia** -- margine >= 1 ma filtro non superato. La PFT è degenerata, non è stata
  battuta. Vanno letti come un artefatto della forma della soglia a piccoli `z_fc`, non come
  un segnale.

La colonna **anche filtro superato** della tabella che segue è dunque il conteggio operativo.
Prendere invece il conteggio grezzo del margine sovrastimerebbe di parecchie volte l'area
candidata.""",
# --- daily-report-only blocks ---------------------------------------------------------
"sup_txt": """Sul suolo classificato l'interpolazione è supportata nel **{pct:.0f}%** delle colonne al culmine della giornata (cioè quella quota contiene almeno {pts} livelli dentro lo strato; le restanti restituiscono un gradiente nan e abbandonano la mappa delle classi, ed è anche ciò che assottiglia i pannelli serali).""",
"sup_none": "La quota di colonne che supportano l'interpolazione non è stata registrata per questa corsa.",
"ml_para_hybrid": """La **stabilità dello strato rimescolato** è qui una misura autentica, non un proxy. I livelli
modello nativi ICON-EU collocano ~10 livelli dentro lo strato rimescolato (in questa corsa: {n}
in mediana sulla terraferma fra le 12 e le 15Z, quando lo strato è maturo -- i conteggi notturni
e serali sono assai più bassi, ma in quelle ore non si classifica nulla). {sup} Il dtheta/dz
dello strato rimescolato è dunque una vera interpolazione ai minimi quadrati su quei livelli.
Validata sui radiosondaggi IGRA (GLA 12Z, intero periodo di registrazione), l'interpolazione sui
livelli modello riduce il bias del gradiente a ~-0,4e-4 K/m (da ~-4,7e-4 sui 5 livelli di
pressione ICON-2I) e il bias dell'ABL da Rib a ~-70 m (da ~-700 m). È questa la ragione per cui
il prodotto ibrido esiste: i dati aperti ICON-2I a 2,2 km non risolvono lo strato rimescolato,
ICON-EU sì.

Il prezzo è la risoluzione orizzontale: l'atmosfera è a 6,5 km (riportata sulla griglia a
2,2 km), mentre il **filtro combustibile conserva i campi superficiali ICON-2I a 2,2 km**, dove
il dettaglio orografico conta. Il gradiente dello strato rimescolato è ora misurato, ma la
soglia di 1,1e-3 K/m ricade ancora dentro la barra d'errore validata (+/- ~2,6e-4): i *conteggi*
delle classi vanno perciò trattati come indicativi, non esatti.

Dopo il tramonto lo strato su cui il gradiente viene interpolato è lo **strato residuo** --
l'aria quasi neutra che lo strato convettivo in decadimento lascia dietro di sé -- anziché il
sottile strato stabile notturno che una particella riferita alla superficie troverebbe. Di
giorno i due coincidono. È questo che consente alle ore serali di portare una previsione; è una
diagnostica più recente del resto della scala, e quelle ore vanno pesate di conseguenza.""",
"ml_para_icon2i": """La diagnostica di **stabilità dello strato rimescolato** è l'anello debole di questo prodotto,
ed è un *proxy*, non una misura. Validata sui radiosondaggi IGRA (GLA 12Z, intero periodo di
registrazione, n=2835), una colonna a 5 livelli di pressione contiene solo 0-2 livelli modello
dentro lo strato rimescolato -- uno o nessuno nel 92% delle colonne costiere -- e il più basso
cade nello strato superficiale superadiabatico. Il dtheta/dz dello strato rimescolato **non è
quindi misurabile da questo archivio**: interpolarlo sui livelli interni allo strato non ha
alcuna capacità discriminante (J di Youden ~ 0,00), e misurare theta fino alla sommità dell'ABL
da Rib ingloba il salto di entrainment (Delta-theta), che Castellnou et al. (2022, sez. 2.1.1)
trattano come variabile *distinta* dal gradiente su cui la scala condiziona.

Si usa invece come proxy la **profondità di rimescolamento della particella**: alla soglia di
1,1e-3 K/m il criterio si riduce a «strato ben rimescolato profondo almeno 455 m». È l'unico
candidato con capacità discriminante (J = 0,29 nell'interno / 0,55 sulla costa, r = +0,50
rispetto al vero dei radiosondaggi). **Sovrassegnala**: la specificità nell'interno è 0,29,
sicché fra le colonne realmente stabili ne dichiara ancora ~71% capaci di pyroCu (sensibilità
0,99 -- raramente manca una colonna capace).

**Conseguenza: le classi vanno lette come segnale di screening, non come conteggi calibrati.**
I *totali* di classe su queste mappe non sono quantitativamente affidabili. Il prodotto ibrido
(`PYROCONV_SOURCE=hybrid`, livelli modello ICON-EU) misura correttamente questo gradiente; su
questa sorgente a 5 livelli si imposti `PYROCONV_ML_METHOD=surface_to_abl` o `mid_layer` per le
misure superate.""",
"ml_rows_hybrid": """| dtheta/dz strato risc. (minimi quadrati, dentro lo strato) | > 1,1e-3 K/m (stabile) | Solo pennacchio convettivo -- nessun pyroCu |
| dtheta/dz strato risc. (misurato sui livelli modello) | <= 1,1e-3 K/m | Colonna capace di pyroCu (bias ~-0,4e-4 K/m vs radiosondaggi) |""",
"ml_rows_icon2i": """| *Proxy* di stabilità: profondità di rimescolamento della particella (v. Metodo) | profondità < 455 m («stabile») | Solo pennacchio convettivo -- nessun pyroCu |
| Proxy di stabilità | profondità >= 455 m (= dtheta/dz <= 1,1e-3 K/m) | Colonna capace di pyroCu (filtro debole: spec. 0,29 nell'interno) |""",
"forcing_hybrid": """Forzanti: atmosfera dai livelli modello nativi ICON-EU a 6,5 km (dati aperti DWD, CC-BY), i ~24 livelli più bassi; ~{n} qui dentro lo strato rimescolato (terraferma, 12-15Z), riportati sulla griglia ICON-2I a 2,2 km. Campi superficiali e filtro combustibile da ICON-2I 2,2 km (MISTRAL / AgenziaItaliaMeteo). Classificatore: pyflam.pyroconvection_type (ABL da Richardson bulk, LCL di Bolton, dtheta/dz misurato nello strato rimescolato, gamma-theta del tappo, RH alla sommità dell'ABL; nessuna CAPE superficiale).""",
"forcing_icon2i": """Forzanti: GRIB ICON-2I 2,2 km sull'intera Italia (MISTRAL / AgenziaItaliaMeteo, CC-BY), sottoinsieme Toscana: geopotenziale, temperatura, RH e vento su {n} livelli di pressione più lo stato a 2 m / 10 m, pressione al suolo e orografia. Classificatore: pyflam.pyroconvection_type (ABL da Richardson bulk, LCL di Bolton, proxy del dtheta/dz dello strato rimescolato, gamma-theta del tappo, RH alla sommità dell'ABL; nessuna CAPE superficiale).""",
"shear_partial": """**Scala effettivamente usata in questa corsa: `{lad}`.** La quinta diagnostica -- la distanza fra ABL/LCL e la quota di massimo shear del vento -- si è risolta solo su una parte del dominio, sicché lì è girata la scala completa a 5 diagnostiche e altrove quella ridotta. Quella clausola è una condizione *necessaria* per la classe più alta, perciò la classe 4 è in qualche misura **più facile** da raggiungere nelle celle ricadute sulla scala ridotta che nelle altre. La classe 4 va trattata come un allerta a ispezionare la colonna, non come una probabilità calibrata.""",
"shear_full": """**Scala effettivamente usata in questa corsa: `{lad}`.** È girato il metodo completo a 5 diagnostiche: i livelli modello hanno risolto una quota di massimo shear, sicché la classe 4 ha richiesto in aggiunta che quel massimo si collocasse entro 0,30 ABL dall'ABL/LCL -- una condizione *necessaria* che le scale ridotte omettono. La classe 4 porta qui dunque la propria clausola di shear; resta comunque un allerta a ispezionare la colonna, non una probabilità calibrata.""",
"shear_none": """**Scala effettivamente usata in questa corsa: `{lad}`.** Il metodo completo prevede una quinta diagnostica -- la distanza fra ABL/LCL e la quota di massimo shear del vento -- che questa sorgente non può risolvere, sicché gira la scala ridotta. La clausola di shear è una condizione *necessaria* per la classe più alta, e ometterla rende la classe 4 in qualche misura **più facile** da raggiungere qui che nel metodo completo. La classe 4 va trattata come un allerta a ispezionare la colonna, non come una probabilità calibrata.""",
"shear_row_partial": "| Distanza del massimo di shear / ABL | <= 0,30 | Richiesta per la classe 4 dove la quota di shear si è risolta -- **applicata su parte del dominio** (v. Metodo) |",
"shear_row_full": "| Distanza del massimo di shear / ABL | <= 0,30 | Richiesta per la classe 4 -- **risolta e applicata in tutta questa corsa** |",
"shear_row_none": "| Distanza del massimo di shear / ABL | <= 0,30 | Richiesta per la classe 4 **solo nella scala a 5 diagnostiche** (non risolvibile qui -- v. Metodo) |",
"ladder_none": "nessuna (nessuna cella classificabile)",
"h_howto": "Come si legge questo prodotto",
"h_method": "Il metodo e i suoi limiti (da leggere prima di usare le classi)",
"p_method": """La colonna è classificata a partire da un **profilo verticale reale**, non da
un'atmosfera standard: {heights}, la profondità di rimescolamento dal **numero di Richardson
bulk** (primo Rib >= 0,33 sopra i 200 m dal suolo, riferito allo stato a 2 m / 10 m), l'LCL
dalla formula esatta di **Bolton (1980)**, e l'umidità alla sommità dell'ABL dalla RH del
modello. Il gamma-theta del tappo è preso fra ABL+200 m e ABL+1200 m.""",
"heights_hybrid": "quote per cella dalle quote dei livelli modello ICON-EU (HHL)",
"heights_icon2i": "quote per cella dal geopotenziale del modello (FI)",
"p_predisposition": """Le classi esprimono la predisposizione atmosferica **dato un incendio di potenza
sufficiente**; nel pannello filtrato quella potenza è calcolata, non presupposta. Sono soglie
derivate dalla letteratura, non validate localmente.""",
"h_classscale": "Scala delle classi (da minima a massima attività piroconvettiva)",
"classscale": """| Livello | Colore | Classe | Significato |
|:--:|:--|:--|:--|
| 0 | bianco | Pennacchio superficiale | Pennacchio di fumo galleggiante; nessuno sviluppo nuvoloso significativo. |
| 1 | verde | Pennacchio convettivo | Il pennacchio penetra uno strato rimescolato stabile; condensazione possibile, nessun pyroCu. |
| 2 | giallo | PyroCu overshooting | Pirocumulo di breve durata; base della nube sopra l'altezza di rimescolamento (LCL/ABL > 1). |
| 3 | arancione | PyroCu persistente | Pirocumulo persistente in colonna instabile (LCL/ABL < 1). |
| 4 | rosso scuro | PyroCu profondo / pyroCb | Piroconvezione profonda / pirocumulonembo; il tappo superiore debole lascia approfondire il pennacchio. |""",
"h_thresholds": "Soglie di classificazione (scala diagnostica in uso: `{lad}`)",
"thr_head": """| Diagnostica | Soglia | Effetto |
|:--|:--|:--|""",
"thr_rest": """| Rapporto LCL / ABL | 1,0 -- 1,60 | PyroCu overshooting (breve) |
| Rapporto LCL / ABL | < 1,0 | PyroCu persistente |
| Rapporto LCL / ABL | <= 1,10 (+ condizioni sotto) | Ammissibile per pyroCu profondo / pyroCb |
| Gamma-theta del tappo (ABL+200 m -> ABL+1200 m) | <= 4,2e-3 K/m (tappo debole) | Consente l'approfondimento a pyroCb |
| Gamma-theta del tappo | >= 4,8e-3 K/m (tappo forte) | Inibisce l'approfondimento (al più persistente) |
| RH alla sommità dell'ABL (media, ABL +/- 150 m) | >= 60% | Richiesta per le classi 3 e 4 (era 80%; rilassata per condizioni di fire weather secche) |""",
"thr_tail": """| Intensità del fronte (filtro combustibile, pannello filtrato) | >= 10 MW/m | Potenza minima per qualsiasi pyroCu (Tedim et al. 2018) |
| Profondità dell'ABL | < {abl:d} m | Non classificabile (profondità implausibile). Nessun filtro a 600 m: non aveva base in Castellnou et al. (2022) e scartava 5 degli 8 incendi di campagna -- rimosso il 26/07/2026 |
| Livelli di pressione utilizzabili | < 4 | Cella non classificata |""",
"h_refcases": "Casi di riferimento (Castellnou et al. 2022, Tabella 1)",
"p_refcases": """Gli eventi etichettati nell'articolo, che ancorano la scala pubblicata a 3 diagnostiche
(`pyflam.pyroconvection_type(ladder="castellnou")`, tuttora predefinita nella libreria e
sottoposta a test di regressione su queste righe). Le scale su profilo usate per questa mappa
sono più severe in cima: richiedono in aggiunta una sommità dell'ABL umida e un LCL vicino ad
essa.""",
"refcases": """| Caso | Tipo osservato | LCL/ABL | dtheta/dz strato risc. | gamma-theta (700-500) |
|:--|:--|:--:|:--|:--:|
| T21 | Pennacchio convettivo | -- | stabile | -- |
| SCQ32 | PyroCu overshooting | > 1 | neutro/instabile | -- |
| M11 | PyroCu persistente | < 1 | instabile | 4,2e-3 (tappo persistente) |
| SCQ41 | PyroCu, non pyroCb | < 1 | instabile | 5,1e-3 (tappo forte) |
| SCQ51 | PyroCu profondo / pyroCb | < 1 | instabile | 3,9e-3 (tappo debole) |""",
"p_footer_daily": """Filtro: Rothermel + chioma Cruz-2005 sui combustibili .lcp con umidità e vento previsti.
Raster diagnostici orari accompagnano le classi: ABL, parcel_ml, residual_ml, LCL, LCL/ABL,
dtheta/dz dello strato rimescolato, gamma-theta del tappo, RH alla sommità, fireABL, decoupling
e (solo ibrido) delta_theta, firecape, penetration, pft_gw (PyroCb Firepower Threshold, GW),
z_fc, delta_theta_fc, u_ml, firepower_gw, pft_margin.
Confini provinciali: derivati ISTAT (openpolis geojson-italy).
Generato da tests/pyroconv_daily.py con pyflam `{ver}`, {at} -- congelato con la fisica di quel
commit, poiché DWD ritira la corsa dopo ~24 h.""",
"h_intro_3day": "Previsione a 3 giorni",
"p_intro_3day": """Tre giorni validi dalla singola corsa {rundate} {run:02d}Z (giorno 0, +1, +2), usando gli step
previsionali a +1/+2 giorni degli stessi dataset ICON-EU + ICON-2I. Ogni giorno è un **ciclo
completo di 24 h, triorario, che inizia alle 00Z del giorno di corsa** ({hours}), sicché le tre
righe-giorno di ogni figura sono direttamente confrontabili e il giorno 0 si apre sull'istante
iniziale della corsa stessa. I pannelli coprono il dominio toscano con i confini provinciali
ISTAT; il mare è mascherato. Atmosfera dai livelli modello nativi ICON-EU (ABL da numero di
Richardson bulk, LCL di Bolton, dtheta/dz misurato nello strato rimescolato, gamma-theta del
tappo, umidità relativa alla sommità dell'ABL, shear); il filtro combustibile a 10 MW/m usa i
campi superficiali ICON-2I a 2,2 km sui combustibili .lcp toscani. Per la validazione
diagnostica (radiosondaggi + ERA5) si veda scripts/validation/README.md.""",
"p_decoup": """Il **rapporto di disaccoppiamento** fireABL / ABL indica quanto in alto un incendio intenso di
riferimento (un eccesso di temperatura di {k:.0f} K nel pennacchio, l'estremo intenso delle
misure in-plume GRAF, 0,1-13,1 K) accrescerebbe il proprio strato limite per *solo calore
sensibile*, diviso per l'ABL ambientale. È la controparte **secca** delle mappe umide delle
classi qui sopra (Castellnou et al. 2022; Castellnou Ribau et al. 2024): valori ben superiori a
1 segnalano colonne profonde, calde e secche, in cui un incendio può sfondare e disaccoppiarsi
dalla superficie *anche dove la scala umida assegna punteggi bassi* -- che è esattamente la
situazione che la mappa filtrata sottostima. Come la mappa del potenziale, presuppone un incendio
ovunque -- qui con un flusso di riferimento fisso anziché calcolato -- sicché è un limite
superiore, non una previsione; non si applica alcuna distinzione LCL secco/umido (lo scarto di
+1 km presente in letteratura non è supportato dalle etichette prototipo GRAF). fireABL da
`pyflam.atmosphere.fire_induced_abl_grid`, la particella intersecata con la reale colonna
theta(z); le celle che il classificatore rifiuta (nessuna colonna utilizzabile, o ABL sotto i
{abl:d} m) restano bianche anziché essere dipinte con un quoziente nudo.""",
"h_dtbl": "Disaccoppiamento secco -- ore diurne, % del suolo classificabile",
"dtbl_head": "| Giorno | Ora | % griglia classificabile | % area rapporto >= 2 | % area rapporto >= 3 | rapporto max |",
"p_dtbl": """*% griglia classificabile* è la quota di dominio con una colonna utilizzabile a quell'ora (ABL
sopra i {abl:d} m); le due colonne del rapporto sono percentuali **di quella quota**, non della
Toscana. Il denominatore collassa verso il tramonto man mano che lo strato rimescolato decade,
sicché una riga delle 18Z può segnare 100% su una piccola area residua: va letta insieme alla
colonna del classificabile.""",
"h_ctbl": "Distribuzione delle classi filtrate per combustibile -- ore diurne, % delle celle di terra ({n} celle)",
"ctbl_head": "| Giorno | Ora | superficiale | convettivo | overshooting | persistente | profondo | n/c |",
"p_ctbl": """Classi: 0 pennacchio superficiale, 1 pennacchio convettivo, 2 pyroCu overshooting, 3 pyroCu
persistente, 4 pyroCu profondo/pyroCb. **n/c** è suolo privo di colonna classificabile (nessun
profilo utilizzabile, oppure ABL sotto i {abl:d} m) -- cosa distinta dalla classe 0, che è invece
il responso della scala diagnostica: la colonna sostiene solo un pennacchio superficiale. Verso
il tramonto n/c occupa la maggior parte del dominio; quelle ore non portano previsione, e non
vanno lette come ore tranquille.""",
"h_evening": "Colonne serali e strato residuo",
"p_evening": """Dopo il tramonto lo strato convettivo collassa in un sottile strato stabile, ma l'aria
soprastante conserva il profilo quasi neutro della giornata: è lo **strato residuo**. Una
particella riferita alla superficie trova solo lo strato stabile, e un dtheta/dz di strato
rimescolato calcolato dentro quella profondità non ha livelli su cui lavorare, sicché l'intero
dominio finiva per cadere fra i non classificabili. Il gradiente viene ora calcolato sulla
maggiore fra la profondità della particella e quella dello strato residuo
(`atmosphere.residual_layer_grid`), che di giorno si riduce esattamente alla profondità della
particella e differisce solo di sera. **Questo cambia le classi serali**: quelle ore ora portano
una previsione dove prima non ne portavano alcuna. Vanno trattate con più cautela di quelle
centrali -- la diagnostica è più recente e la profondità residua che restituisce appare ancora
sottile rispetto a quanto suggerisce il profilo.

Due ulteriori diagnostiche viaggiano insieme e **non** sono usate dal classificatore: il salto
di zona di entrainment `delta_theta` (il gradino di temperatura potenziale che un pennacchio deve
superare per uscire dallo strato rimescolato) e `firecape` (Potter 2005, la CAPE di una
particella riscaldata dall'incendio). Entrambe provengono dall'insieme di variabili del metodo
originale; nessuna delle due filtra una classe qui, perché il criterio di penetrazione che
implicano non discrimina ancora alla potenza di incendio di riferimento di questo prodotto. Si
veda `docs/graf_vs_pyflam_2026-07-26.md`.""",
"p_footer": """I GeoTIFF orari delle classi e delle diagnostiche accompagnano questo bollettino in
`rasters_hybrid_<giorno>/`: ABL, parcel_ml, residual_ml, LCL, LCL/ABL, dtheta/dz dello strato
rimescolato, gamma-theta del tappo, RH alla sommità, fireABL, decoupling, delta_theta, firecape,
penetration, oltre al lato fuoco -- fuel_load, burnable_fraction, firepower_gw, pft_gw,
pft_margin, crit_growth_ha_h, z_fc, u_ml, delta_theta_fc.
Generato da scripts/forecast_3day_tuscany.py con pyflam `{ver}`, {at}.
Questo prodotto non è rigenerabile una volta che DWD ritira la corsa (~24 h), ed è quindi
congelato con la fisica di quel commit.""",
"mtbl_head": "| Giorno | Ora | % griglia con potenza | % margine >= 1 | % margine >= 0,1 | margine max | celle >= 1 | **anche filtro superato** |",
"mtbl_head_daily": "| Ora | % griglia con potenza | % margine >= 1 | % margine >= 0,1 | margine max | celle >= 1 | **anche filtro superato** |",
"p_mtbl": """*% griglia con potenza* è la quota di dominio che ha insieme una colonna classificabile e
combustibile bruciabile; le due colonne del margine sono percentuali **di quella quota**, non
della Toscana. Una percentuale sull'intero dominio sarebbe dominata dalle celle prive di
combustibile, che non potranno mai contribuire, e si muoverebbe per ragioni estranee alla
previsione. *celle >= 1* è il conteggio grezzo dietro la prima percentuale, riportato perché a
questi tassi una percentuale si arrotonda a un valore che sembra nulla quando in realtà si
tratta di una manciata di luoghi precisi. ***anche filtro superato*** applica il criterio
congiunto della sezione precedente ed è il conteggio su cui agire; lo scarto fra esso e
*celle >= 1* è la popolazione a PFT degenerata.""",
},
}


def t(lang, key, **kw):
    """One localized block, formatted. Falls back to English on a missing key rather than
    raising: a report with one English paragraph is a far smaller failure than no report."""
    s = MD.get(lang, MD["en"]).get(key) or MD["en"][key]
    return s.format(**kw) if kw else s


def f(lang, key, **kw):
    """One localized figure label, formatted. Same fallback rule as :func:`t`."""
    s = FIG.get(lang, FIG["en"]).get(key) or FIG["en"][key]
    return s.format(**kw) if kw else s
