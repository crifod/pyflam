# Simulated Peer Review — *pyflam: A Scientific, Technical and Operational Report*

*Simulated multi-perspective peer review (EIC + 3 reviewers + Devil's Advocate) generated as a pre-submission diagnostic. Venue framing: a computational-model/software paper for a fire-science or geoscientific-model-development journal (GMD / Environmental Modelling & Software / International Journal of Wildland Fire). This is a simulation to surface likely objections — not an actual journal decision.*

**Manuscript reviewed:** `docs/pyflam_scientific_report.md` (732 lines)
**Overall simulated decision:** **Major Revision** (as a *model-description* paper) — or reframe/split (see EIC letter).

---

## Phase 0 — Reviewer Configuration Card

| Reviewer | Simulated identity | Lens |
|---|---|---|
| **EIC** | Editor, geoscientific model-description journal; publishes open models with mandatory evaluation + code/data DOI | Scope, originality, evaluation completeness, reproducibility |
| **R1 — Methodology** | Numerical analyst (Hamilton–Jacobi/eikonal solvers) + model-evaluation statistics | Correctness of methods; verification vs validation; statistics; convergence |
| **R2 — Domain** | Operational fire-behavior scientist; coupled fire–atmosphere modelling (WRF-SFIRE, FIRETEC, ForeFire) | Physical fidelity of the novel fire-science claims; literature currency |
| **R3 — Perspective** | Research-software engineer + fire-management decision-support practitioner | Software quality, usability, the open-source competitive landscape, operational value |
| **Devil's Advocate** | Skeptic of "validated operational tool" claims | Strongest counter-argument; logical gaps; "so what?" |

---

## Phase 1 — Independent Reviews

### Reviewer EIC — Editorial Assessment

**Summary of submission.** A single-author paper describing pyflam, an open-source, cross-platform Python reimplementation of FlamMap-class landscape fire behavior, extended with weather-driven per-cell fuel moisture, native terrain-wind and buoyant-CFD plume coupling, a Finsler-eikonal spread solver, a literature-current crown-fire model, physics-based spotting, and vertical-profile pyroconvection diagnostics.

**Strengths.**
- Genuinely broad, coherent scientific scope, each layer tied to cited primary literature.
- The deterministic surface core is validated cell-by-cell against a real FlamMap run (ROS ~3%, direction ~1°, conditional fireline intensity ~2%) — this is a real, quantitative result.
- Strong software-engineering posture (≈568 tests, CI on 3 Python versions, optional-dependency gating, open license).
- Transparent about what is and isn't done (§11 labels open items; AI-assistance disclosed).

**Editorial concerns.**
1. **Genre mismatch (MAJOR).** This is written as a *technical/operational report*, not a journal article. There is no Abstract / Introduction / Methods / Results / Discussion / Conclusions structure, no research question, and the tone is promotional in places ("robust," "goes well beyond," "most distinctive contribution"). For a model-description journal it must be restructured; promotional adjectives must go.
2. **Evaluation completeness (MAJOR / BLOCKING for GMD-tier).** §11 is titled "Validation *strategy*." For a model-description journal, an evaluation *plan* is not publishable — the model must be evaluated. The validated part (surface core vs FlamMap) is only the least-novel layer; the six headline novelties are not validated against independent data.
3. **Single evaluation case (MAJOR).** All quantitative results come from one landscape (Tuscany) under essentially one weather condition. Generalization claims need ≥2–3 contrasting landscapes/climates.
4. **Reproducibility artifact (MAJOR).** No archived, version-pinned code+data DOI (Zenodo). The FlamMap reference rasters that underpin every validation number are not shown to be shareable — if they cannot be released, the central results are not reproducible by reviewers.
5. **Scope-vs-depth tension.** Seven novelties in one paper means none is evaluated in depth. Consider splitting (see recommendation).

**Recommendation:** **Major Revision** as a model-description paper — *conditional on adding real evaluation of at least the flagship novelties*. If the author's near-term goal is publication of the *software*, EIC would redirect to a software-paper venue (JOSS/SoftwareX), which this manuscript could satisfy after restructuring, and reserve the science claims for a later, evaluation-backed submission.

---

### Reviewer 1 — Methodology

Focus: are the methods correct, and is the evidence *validation* or merely *verification*?

**Positive.**
- The surface-core validation is methodologically sound: cell-by-cell diff, OLS slope/r, log-space and within-X% statistics, burn/no-burn classification. `validate.py` is well-designed and the circularity risk of wind-fitting (`scan_parameter`) is explicitly acknowledged and quarantined — good practice.
- The eikonal formulation is correct: fire arrival time as a static anisotropic Hamilton–Jacobi equation with a Randers–Finsler metric, semi-Lagrangian causal update (Sethian–Vladimirsky / Mirebeau lineage). The theory is right.

**Major methodological issues.**
1. **Verification ≠ validation, and the paper conflates them (CRITICAL).** Every novel claim rests on one of: an *analytic* benchmark on a *uniform synthetic* field (eikonal), a *coefficient* check against a source paper (Cruz, VPD), or a *single plausibility number* (plume: "mean wind 2.43→3.14 m/s"; feedback: "no runaway"). None is a comparison to *observations* or an *independent model*. These establish the code computes the intended equations — not that the equations predict reality. The manuscript's language ("validated," "confirmed") overstates this.
   - *Fix:* systematically relabel verification vs validation; add real validation (below) or explicitly scope claims down to "verified implementation."
2. **Eikonal evaluation is insufficient (MAJOR).** ~4% vs ~12% bias is shown only on a *uniform* field against an analytic solution. Required: (a) a **convergence study** (error vs grid resolution — does it converge at the claimed order?); (b) evaluation on **heterogeneous** landscapes where the metric varies cell-to-cell (the realistic case); (c) a like-for-like cost/accuracy Pareto vs MTT at matched accuracy, not just matched grid.
3. **Coupling and pyroconvection have no quantitative evaluation (MAJOR).** "2.43→3.14 m/s" is a sanity check, n=1. There is no comparison to a coupled reference (WRF-SFIRE, FIRETEC, or observed plume/indraft data). The pyroCb diagnostics (C-Haines, inverted-V, firepower threshold) are presented with zero event-based skill assessment (no hit-rate / false-alarm on known pyroCb vs non-pyroCb days).
4. **Crown fire has no independent validation (MAJOR).** The paper *deliberately* does not validate against FlamMap (defensible) but then validates against *nothing* of its own — it "inherits" Cruz et al.'s validation. Inheriting a source model's validation is not validating pyflam's implementation on independent fires.
5. **Stochastic spotting needs distributional validation (MAJOR).** A Monte-Carlo model must be evaluated distributionally (e.g. predicted vs observed spot-distance distribution, K-S / QQ), not by calibrating "one weakly-constrained length" to anchors.
6. **No uncertainty quantification (MAJOR).** Point agreement statistics only; no CIs, no sensitivity analysis to the many parameters (χ_c=0.6, shade factor, ΔT_sun, the calibrated spotting length).
7. **Apparent internal inconsistency (MINOR but must fix).** Flame length is given as `L = 0.45·I^0.46` in §2 and `L = 0.0775·I^0.46` in Appendix A.1. Presumably different unit conventions (Byram ft/BTU vs metric), but as written it reads as a contradiction — state units and reconcile.

**Recommendation:** Major Revision. The core is sound; the novel claims are not yet supported at the level the wording asserts.

---

### Reviewer 2 — Domain (fire behavior & coupled fire–atmosphere)

**Positive.**
- Crown-fire critique is current and correct: the Cruz & Alexander (2010) under-prediction argument is real, and offering Cruz 2005 + CFIS 2004 as selectable alternatives is a genuine operational contribution.
- The pyroCb framing is scientifically up-to-date: recognizing that deep pyroconvection is a *vertical* (mid-tropospheric moisture / inverted-V / C-Haines) problem rather than a surface-CAPE problem is the correct modern view (Peterson, Tory & Kepert, Castellnou). Wiring a per-cell convective plume factor to this geometry is a nice idea.

**Major domain issues.**
1. **The competitive landscape is mis-set (MAJOR).** The paper positions pyflam almost entirely against FlamMap (closed, Windows). But the relevant scientific comparators are the *open* and *coupled* models: **Cell2Fire, ELMFIRE, SPARK, ForeFire** (open landscape spread) and **WRF-SFIRE, FIRETEC, QUIC-Fire, ForeFire/MesoNH** (coupled fire–atmosphere). Several are open-source and directly overlap pyflam's claims. Almost none are cited. The novelty claims ("absent from the operational tools," "to our knowledge a contribution not present") are only defensible against FlamMap, not against this field — and will be challenged.
2. **"Coupled fire–atmosphere" is a strong claim for a RANS/Boussinesq steady plume (MAJOR).** WRF-SFIRE/FIRETEC solve transient, compressible/large-eddy dynamics. A steady `buoyantBoussinesqSimpleFoam` plume re-solved every dt is a *quasi-steady diagnostic* coupling. That's legitimate and cheaper, but the paper should (a) say so plainly, (b) justify the quasi-steady assumption for fast-evolving plumes, and (c) show it against at least one transient reference or observation.
3. **χ_c ≈ 0.6 convective fraction is asserted, not sourced per-fuel (MODERATE).** Convective/radiative partitioning varies widely; a single constant drives the entire plume forcing. Needs a citation and a sensitivity test.
4. **GEDI-derived CBH/CBD (MODERATE).** Deriving canopy base height and bulk density from height+cover via heuristics is reasonable and honestly flagged, but it injects large, uncharacterized uncertainty into exactly the crown-fire results used to argue the Cruz-vs-Rothermel difference. That difference may partly reflect the CBD heuristic, not the spread model. Disentangle.
5. **Live-fuel and its role in pyroCb/crown (MINOR).** Live fuel moisture treatment is thin relative to its importance in crown fire and plume energy.

**Recommendation:** Major Revision. Strong, current fire science, but the novelty framing needs the open/coupled-model literature, and "coupled" must be calibrated to what the solver actually does.

---

### Reviewer 3 — Perspective (research software + operational practice)

**Positive.**
- As *software*, this is strong: cross-platform, scriptable, community formats (.lcp/.fms/GeoTIFF/GeoJSON), tested, CI, runtime-optional heavy engines. The operational NRT layer (`meteo_report`, `operative`, `nrt.run_realtime`) with GeoJSON export is genuinely useful and deployable.
- Reproducible end-to-end scenario (`final_run_tuscany.py`) on real 5.35M-cell data is a real asset.

**Issues.**
1. **Wrong venue risk (MAJOR, constructive).** The manuscript's greatest, most-defensible contribution *today* is the open, tested, deployable tool — which is exactly what **JOSS/SoftwareX** exist to publish, and which does not require the scientific-novelty validation that a model journal demands. Submitting the current text to GMD invites rejection on evaluation grounds; submitting a restructured version to JOSS is a near-fit now.
2. **Software availability specifics missing (MAJOR).** No DOI, no explicit version tag in the paper, no archived environment (exact dependency versions), no statement of how the FlamMap reference data can be obtained. For any of these venues, add a Code/Data Availability section with a Zenodo DOI for a tagged release.
3. **Usability evidence is asserted, not shown (MINOR).** Claims of operational value would be far stronger with one worked decision-support case and, ideally, a practitioner's account.
4. **Decision-support risk framing (MINOR).** A tool aimed at suppression decisions should state validation limits and intended-use caveats prominently (liability/ethics of operational fire predictions).

**Recommendation:** Accept-ish *as a software paper after restructuring* (JOSS/SoftwareX); Major Revision as a science paper.

---

### Devil's Advocate

**Strongest counter-argument.**
> The paper's central, repeated claim is that pyflam is a *validated* operational fire tool that *goes beyond* the established models. Strip the language and what is actually validated is the one layer that is **not** novel — the Rothermel surface core, checked against FlamMap, i.e. against *the very tool pyflam reimplements*. That is a reproduction check, not a scientific advance. Every genuinely novel component — the eikonal solver, the fire–atmosphere coupling, pyroCb diagnostics, the crown correction, physics-based spotting — is supported only by internal verification (analytic toy problems, coefficient checks, single sanity numbers, "no runaway"). So the paper is strongest exactly where it is least original, and weakest exactly where it claims to contribute. Until a novel component is shown to predict a *real fire* better than the incumbent it critiques, the headline contribution is *potential*, not demonstrated.

**Issue list.**
- **CRITICAL** — Validation is concentrated on the non-novel layer; novel claims unvalidated against observations. (Contribution/Evidence; §§4.2, 6, 7, 8, 11) → per review IRON RULE, this alone precludes an "Accept."
- **CRITICAL** — "Correcting the operational under-prediction" (§6 title) is asserted but never demonstrated on data: pyflam is not shown to match observed crown-fire spread better than FlamMap. It reproduces Cruz's equations; it does not independently show the correction is real *in pyflam*.
- **MAJOR** — Novelty claims ("absent from the operational tools," "to our knowledge not present") are FlamMap-relative and ignore open/coupled competitors (Cell2Fire, ELMFIRE, ForeFire, WRF-SFIRE). Some claims may not survive contact with that literature.
- **MAJOR** — "So what?" for the eikonal ~4% vs ~12%: on a synthetic uniform field. Does the bias reduction change any operational decision on a real landscape? Not shown.
- **MAJOR** — Single landscape + single weather = one data point; robustness/overfitting to Tuscany unaddressed.
- **MODERATE** — Confirmation-friendly presentation: results are reported where favorable (ROS 3%, direction 1°) and downgraded to "partially recoverable" / "open" where not (burn probability, perimeter/arrival, crown). Fine, but the abstract/executive summary should carry the negatives with equal prominence.

**Missing stakeholder perspective.** Incident commanders relying on spread/spotting predictions — the decision-risk of an unvalidated coupled model used operationally is not discussed.

**Observation (non-defect).** The honesty of §11 and the AI disclosure is commendable and works in the author's favor.

---

## Phase 2 — Editorial Synthesis & Decision

### Consensus (all reviewers agree)
1. The **surface core is genuinely validated**; the software engineering is strong.
2. The **six novel contributions are verified, not validated** — this is the decisive, unanimous finding (R1 CRITICAL, DA CRITICAL, EIC MAJOR, R2 MAJOR).
3. **One landscape / one weather** is insufficient for the generalization claimed.
4. **Reproducibility needs a DOI-archived code+data artifact** and clarity on the FlamMap reference data.
5. **Genre/venue**: as written it is a report, not an article.

### Disagreement (arbitrated)
- **R3** would accept it *now* as a software paper (JOSS); **EIC/R1/R2/DA** would not accept it as a *science* paper. → **Arbitration:** both are right about different venues. The manuscript is ~publishable as *software* after restructuring, and *not yet* publishable as a *model-evaluation* paper. The author should choose the target and the claims must match it.

### Decision
Because the Devil's Advocate raised **CRITICAL** issues (IRON RULE #4), the decision cannot be Accept.

**Simulated decision: MAJOR REVISION**, with an explicit fork:

- **Path A — Software paper now (fastest):** restructure to JOSS/SoftwareX format, scope every claim to "verified implementation," add Code/Data Availability + Zenodo DOI, strip promotional tone. *Likely acceptable after minor work.*
- **Path B — Model-description/science paper (higher impact, more work):** keep GMD/EMS/IJWF as target but add real validation of at least 2 flagship novelties, on ≥2 landscapes, with UQ. *This is the Major Revision proper.*
- **Path C — Split:** JOSS software paper now (cite it) + a focused science paper later on **one** validated novelty (recommended: the eikonal solver, which is closest to demonstrable — or the crown correction validated on observed fires).

### Decisive blocking issues (must resolve before a model-journal Accept)
1. **Validate ≥2 headline novelties against observations or an independent model** (not analytic/coefficient/sanity checks). Highest-value, most-tractable targets:
   - *Crown correction:* pyflam vs **observed** crown-fire ROS (a small observed dataset) — directly substantiates the §6 claim.
   - *Eikonal:* convergence study + heterogeneous-landscape accuracy vs a fine-grid reference; decision-relevant impact.
   - *pyroCb diagnostics:* hindcast skill (hit/false-alarm) on a set of known pyroCb vs non-pyroCb fire days.
2. **Second and third evaluation landscapes** in contrasting fuel/terrain/climate.
3. **Reproducibility:** tagged release + Zenodo DOI; state FlamMap-reference-data availability; pin the software environment.
4. **Relabel verification vs validation throughout**; remove promotional language; add UQ/sensitivity for key parameters (χ_c, spotting length, shade/ΔT_sun).
5. **Reposition novelty against the open/coupled-model literature** (Cell2Fire, ELMFIRE, ForeFire, WRF-SFIRE, FIRETEC, QUIC-Fire); qualify "coupled" to the quasi-steady RANS reality.
6. **Fix the flame-length coefficient inconsistency** (§2 vs A.1) with explicit units.

### Prioritized Revision Roadmap
| # | Priority | Action | Addresses |
|---|---|---|---|
| 1 | P0 | Decide venue/path (A/B/C) and align all claims to it | EIC, R3, DA |
| 2 | P0 | Relabel verified vs validated; remove promotional tone | R1, DA, EIC |
| 3 | P0 | Code/Data Availability + Zenodo DOI + env pin + FlamMap data statement | EIC, R3 |
| 4 | P1 | Validate crown correction vs observed fire ROS | DA-CRITICAL, R2 |
| 5 | P1 | Eikonal convergence + heterogeneous-landscape evaluation | R1 |
| 6 | P1 | Add 2nd/3rd landscape | EIC, DA |
| 7 | P1 | pyroCb hindcast skill (if keeping the claim) | R1, R2 |
| 8 | P2 | Reposition vs open/coupled models; qualify "coupled" | R2, DA |
| 9 | P2 | UQ / sensitivity for key constants | R1 |
| 10 | P2 | Distributional validation of spotting; fix flame-length units | R1 |

---

*End of simulated review. Reviewers acted independently; the synthesis cites only points raised in the reports above. This is a pre-submission diagnostic, not an actual editorial decision.*
