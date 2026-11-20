# Gauge Freedom in the Thermodynamics of pyflam's Pyroconvection Diagnostics

**A literature review with a computational audit of reference-state covariance**

Deep-research report (ARS `full` mode) · 2026-11-19 · pyflam at commit `f93b693` (branch `feat/fire-weather-fwi-pipeline`)

---

## Abstract

"Gauge theory of thermodynamics" refers to three distinct bodies of work: (A) the contact-geometric formulation of equilibrium thermodynamics, in which the choice of potential is a gauge freedom; (B) network non-equilibrium thermodynamics, in which entropy production is a gauge-invariant observable of Markov dynamics; and (C) the dependence of atmospheric thermodynamic variables on arbitrary reference constants. This report asks which of the three can change a number that pyflam produces. Readings A and B are reviewed briefly. Both are mathematically sound, and neither connects to an operational quantity in pyflam's pipeline. Reading C does. Treating the reference pressure p₀ of potential temperature as a scale gauge (θ → cθ, c = (p₀′/p₀)^κ) gives a simple weight-counting rule: a formula is reference-independent only if every term carries the same gauge weight. Applying this rule to `src/pyflam/atmosphere.py` identifies three weight mismatches. The PyroCb Firepower Threshold (PFT, Tory & Kepert, 2021) adds a fire moisture increment proportional to θ, outputs a θ difference rather than a temperature difference, and inherits both mismatches into the cost ladder. `fire_cape_grid` adds a temperature-like excess to θ_v. A single-column sensitivity audit puts the effect of the global choice of p₀ at 0.3%, which is negligible. The terrain-correlated part reaches −4% to +9% of PFT over 750–925 hPa surfaces. That is about half the spread caused by the assumed fire heat-to-moisture ratio, and larger than the β-discretisation noise. The recommended action is narrow: declare the gauge once, fix the duplicated Poisson exponent, add a covariance property test, and resolve which reference Tory & Kepert intended before trusting PFT maps over the Pyrenees and Apennines.

**Keywords:** gauge invariance, potential temperature, reference state, pyrocumulonimbus firepower threshold, contact geometry, stochastic thermodynamics

---

## 1. Introduction

pyflam's atmospheric layer is quasi-equilibrium parcel thermodynamics. The core is potential temperature with a fixed 1000 hPa reference and a dry Poisson exponent (`theta_kelvin`, `atmosphere.py:891`). On top of it sit θ_v, a pseudo-adiabatic lapse rate, Bolton's (1980) LCL, Continuous Haines, a dry fire-CAPE, and the PFT and cost ladder built on Tory et al.'s (2018) saturation-point construction. The request was to evaluate gauge theory as a new approach to this thermodynamics.

The word "gauge" goes back to Weyl's 1918 scale-gauge geometry (*Eichung*, calibration of length). In that geometry a change of local length scale shifts a connection by the gradient of the log scale factor (Straumann, 1996). The thermodynamic literature uses the word in looser senses. So the useful question is not whether thermodynamics can be written as a gauge theory, because it can, in several ways. The useful question is whether any of those formulations produces a falsifiable change in a pyflam output.

**Primary research question.** Are pyflam's pyroconvection diagnostics (PFT, cost ladder, fire-CAPE, C-Haines) invariant under the gauge choices of atmospheric thermodynamics? These are the reference pressure, entropy reference values, dry versus moist exponent, and the conserved variable. If they are not invariant, how large is the dependence relative to the uncertainties already in the pipeline?

**Sub-questions.**
- SQ1: Which diagnostics depend on reference constants, and by how much?
- SQ2: Does a reference-free moist entropy formulation (θ_s; Marquet, 2011) matter for these diagnostics?
- SQ3: Do readings A or B offer a testable advantage over C?

## 2. Method

### 2.1 Literature search

Sources were searched with web queries targeting arXiv (cond-mat.stat-mech, physics.ao-ph), *Reviews of Modern Physics*, *EPL*, *Eur. Phys. J. B*, *Entropy*, *QJRMS*, *JAS* and *MWR*. The seed terms were "thermodynamics gauge theory", "contact geometry thermodynamics", "Legendre submanifold", "nonequilibrium gauge", "moist entropy reference state", "equivalent potential temperature computation", and "pyrocumulonimbus firepower threshold". Inclusion required (i) a peer-reviewed or archived-preprint primary source and (ii) bibliographic confirmation (journal, volume, pages or DOI) from at least one independent record. Sources that could not be confirmed were dropped. Abstracts and bibliographic records were checked for every source. Full texts were **not** read for most of them, and Section 6 states where this limits a claim.

### 2.2 Computational audit

The audit has two parts.

1. **Static weight analysis.** Every expression in the PFT chain and in `fire_cape_grid` was assigned a gauge weight under p₀ rescaling (Section 3.3).
2. **Sensitivity runs.** The tests' `_synthetic_column` fixture was used: a well-mixed ABL under a constant-lapse free troposphere. Two cases were run, hot-dry (RH 25%, ABL 1800 m) and moist (RH 45%, ABL 1200 m). Each was shifted to surface pressures of 1000, 925, 850 and 750 hPa, keeping the θ profile and approximate RH. The PFT was then recomputed under five re-gaugings, with the β grid at 60, 600 and 1200 steps. Two physical comparators set the scale: β-grid discretisation, and the fire heat-to-moisture ratio varied from 12 to 18 K per g/kg around pyflam's default of 15.

The scripts and raw output are in `docs/gauge_audit/`.

**Validity criterion (set in Phase 1).** A reformulation is worth adopting only if its effect exceeds the existing noise floor (discretisation) and is a material fraction of the physical uncertainty (heat-to-moisture), or if it improves skill on the reference fires. The last test was not run here (Section 6).

## 3. Theoretical framework

### 3.1 Reading A: contact geometry and Legendre gauge

Equilibrium states form a Legendre submanifold of a contact manifold (thermodynamic phase space). Moving between the energy, enthalpy, free-energy or entropy representations is a contact transformation, so the choice of potential is a gauge freedom of description (Mrugała et al., 1991). Balian and Valentin (2001) make this explicit: a thermodynamic system is a gauge-invariant Lagrangian submanifold of a contact-type space, processes, including dissipative ones, are Hamiltonian trajectories on it, and changes of equation of state are canonical transformations. Bravetti (2017) reviews contact Hamiltonian dynamics in reversible and irreversible thermodynamics. Bravetti et al. (2015) define a conformal gauge transformation of the phase-space structure that maps Weinhold's energy metric onto Ruppeiner's entropy metric. They show that the curvature of the associated thermodynamic connection is *not* conformally gauge invariant. The metric strand itself goes back to Ruppeiner (1995), whose thermodynamic curvature bounds the length scale below which classical fluctuation theory fails.

*Relevance to pyflam.* Low. Parcel thermodynamics in pyflam is already well-posed in one representation (T, p, q). Contact geometry would reformulate the same equations of state without changing a parcel's buoyancy. The one place the strand touches something measurable is the Bravetti et al. (2015) result: curvature-based quantities depend on representation. That is a caution against building diagnostics on thermodynamic-metric curvature, not an invitation to do so.

### 3.2 Reading B: non-equilibrium thermodynamics as a gauge theory

Schnakenberg (1976) built the macroscopic thermodynamics of master-equation systems from graph theory. Steady states come from Kirchhoff's theorem, and cycle affinities serve as generalised forces. Polettini (2012) showed that Markov dynamics on a finite graph has a gauge symmetry, and that the thermodynamic sector, entropy production in particular, follows from requiring observables to be gauge invariant. The driving force acts as a gauge potential, and Clausius's irreversibility along a cycle becomes a geometric phase. Polettini and Vedovato (2013) extend this to diffusions and frame the gauge as the freedom to choose a reference ("subjective") prior. GENERIC (Grmela & Öttinger, 1997) is the continuum counterpart that pairs a reversible Poisson structure with an irreversible dissipative one, though it is not usually called a gauge theory.

*Relevance to pyflam.* Low to moderate in principle, nil in practice. Fire–atmosphere coupling is genuinely non-equilibrium. The fire's entropy production feeds the plume, and the plume feeds back on the fire through indraught wind (`pyroconvection.fire_atmosphere_march`). But pyflam does not model this coupling as a Markov network or a GENERIC system, and none of its outputs is an entropy production. Applying reading B would mean building a new model, not re-expressing an existing one. Polettini and Vedovato's lesson does carry over directly to reading C: a gauge is a choice of reference, and physical outputs must not depend on it.

### 3.3 Reading C: the reference pressure as a scale gauge

With θ = T (p₀/p)^κ and Exner function Π = (p/p₀)^κ, temperature factorises as T = θΠ. Changing the arbitrary reference p₀ → p₀′ acts as a one-parameter dilation:

| Quantity | Transformation | Gauge weight w |
|---|---|---|
| θ, θ_v, Δθ | × c, with c = (p₀′/p₀)^κ | +1 |
| Π | × c⁻¹ | −1 |
| T, ΔT, p, q, z | unchanged | 0 |
| Δθ/θ, β (fractional excess), buoyancy (θ_p − θ_e)/θ_e | unchanged | 0 |

**Covariance rule.** A sum is gauge-covariant only if all its terms have equal weight. A diagnostic is reference-independent only if its total weight is zero. Adding a weight-0 quantity (a measured ΔT, or a constant expressed "per K of temperature") to a weight-1 quantity breaks covariance. The result then depends on p₀, and the error scales with (p₀/p_s)^κ.

The symmetry is *global* as long as p₀ is the same constant everywhere, and pyflam uses 1000 hPa throughout. It becomes a *local* gauge if p₀ varies in space, for example a terrain-relative p₀ = p_s(x). Then ∇θ′ = c(∇θ + θ ∇ln c), and horizontal gradients need a compensating connection A = ∇ln c, with covariant derivative Dθ = ∇θ − θA. This is the same structure as Weyl's length connection (Straumann, 1996). It is the only sense in which atmospheric reference-state choice is *literally* a gauge field. It matters for pyflam only if surface-referenced θ maps are ever differenced horizontally.

**Moist and entropy references.** Moist entropy contains additive reference constants for each component. These cancel in a closed parcel but not when the total-water content changes. Fixing them requires third-law reference values, which is the basis of Marquet's (2011) θ_s. Hauf and Höller (1987) instead adopted a reference state that depends on composition (Marquet & Stevens, 2022). Marquet and Stevens (2022) show that, of θ_e, θ_l and θ_s, only θ_s measures parcel entropy. Separately, the various θ_e formulae (Bolton, 1980; Bryan, 2008) are approximations with comparable errors. They are approximation choices, not gauges. Finally, the dry κ = R_d/c_pd differs from the moist exponent by a fraction of a percent at typical mixing ratios.

## 4. Findings

### 4.1 Static weight analysis (SQ1)

| Location | Expression | Weight check | Verdict |
|---|---|---|---|
| `theta_kelvin` (`atmosphere.py:891`) | T(1000/p)^0.286 | Defines the gauge | OK, but the literal 0.286 differs from `_KAPPA_DRY` = 287.05/1004 = 0.28591 used to invert θ → T in `_mixed_layer_state` and `_saturation_point`. The forward and inverse maps use different exponents. **Fixed 2026-11-19** (Recommendation 1): one gauge, `_P0_HPA` and `_KAPPA_DRY`. |
| `_saturation_point` (`:1808`) | θ_p = (1+β)θ_ML | w = 1 on both sides | Covariant |
| `_saturation_point` | q_p = q_ML + βθ_ML/h, h = 15 K per g/kg | w = 0 (q) + w = 1 (βθ_ML/h, if h is a temperature ratio) | **Breaks.** The fire moisture depends on p₀. Over elevated terrain it is (p₀/p_s)^κ larger than the same ratio expressed per K of surface temperature. *Update 2026-11-19: this is the published convention (Tory & Kepert, 2021, eq 30 and Table 2), which pyflam reproduces. See Rec. 3.* |
| `pyrocb_firepower_threshold_grid` (`:1828`) | Δθ_fc = βθ_ML, PFT ∝ Δθ_fc (eq 31) | Output weight +1 | **Breaks.** PFT scales with c. The global part is uniform. *Update 2026-11-19: the paper fixes P₀ = 10⁵ Pa (Table 2) and derives ΔH = ρ₀c_pΔθ from dry static energy referenced to z(P₀) (eq 10, appendix B). The gauge is part of the published PFT.* |
| `pyroconvection_cost_grid` (`:2344`) | Reuses `_saturation_point` and eq 25 | Inherits both | **Breaks** (same terms) |
| `fire_cape_grid` (`:1416`) | b = (θ_v,sfc + exc − θ_v)/θ_v | Ratio has w = 0, and exc is a θ′ (w = 1) | **Covariant.** *Corrected 2026-11-19:* the GRAF ingest measures θ′ directly (in-plume minus environment θ), and θ′ is conserved along dry adiabats, so it is the right quantity to carry up the column. The earlier "breaks" verdict assumed a ΔT. |
| `fire_parcel_theta_excess` (`:1180`) | w* = (3gQ₀H/2θ_v)^{1/3} with Q₀ in K m/s from c_p | w* ∝ c^{−1/3} | Weak break (≈1.6% in w* at 850 hPa). The docstring already says not to use this as a plume excess. |
| `critical_rh_for_cloud_persistence` (`:1503`) | RH_crit from θ_excess·c_p / (L_v q_s) | Equated a w = 1 excess with a w = 0 latent cooling | **Broke; fixed 2026-11-19.** The excess is now converted, ΔT = Δθ(p/p₀)^κ. RH_crit was 0.5–4.4 points too low (54.7 vs 59.1% for 10 K at 700 hPa). |
| `pyflam_gui/core/pyroconv.classify` | θ_sfc = θ(T2m, **1000 hPa**) | Wrong pressure for the 2 m state over terrain | **Broke; fixed 2026-11-19.** θ_sfc was 4.4 / 9.1 / 14.2 K too cold at 950 / 900 / 850 hPa, 3–10× the ML-stable gate, biasing the Fire Weather Preview quiet over hills. Now uses surface pressure. |
| `continuous_haines` | T, T_d at 850/700 hPa | w = 0 | Invariant. Its fixed isobaric levels are a coordinate choice, not a gauge (it returns nan below ground). |
| `briggs_buoyancy_flux` | gQ/(πρc_pT_ref) | w = 0 | Invariant |

Two minor consistency issues also turned up. The air heat capacity appears as both `_CP_AIR = 1005` and `_CP_DRY = 1004`. And `_KAPPA = 0.4` is the von Kármán constant, as the comment at `:752` notes, which makes the 0.286 literal easy to confuse with it.

### 4.2 Sensitivity runs

All PFT values below are in GW with a 600-step β grid unless stated. Raw output: `docs/gauge_audit/audit_output.txt` and `terrain_output.txt`.

**Global gauge and noise floor (surface at 1000 hPa).**

| Perturbation | Hot-dry ΔPFT | Moist ΔPFT |
|---|---|---|
| p₀ = 1013.25 hPa instead of 1000 | +0.31% | −0.19% to +0.34% |
| Single κ = R_d/c_pd, forward and inverse | −0.04% to +0.26% (now 0 by construction, see Rec. 1) | −0.05% |
| β grid 60 vs 600 steps (noise floor) | 2.3% | 0.4% |
| h = 12 / 18 K per g/kg (physical comparator) | −9.8% / +7.8% | −7.6% / +5.6% |

The global choice of p₀ and the κ mismatch are both below the discretisation noise. On flat terrain they do not matter.

**Terrain-correlated part (1200-step β grid).** Changes are relative to the PFT as currently coded. Re-run after Recommendation 1; values moved by ≤ 0.3 percentage points.

| Surface | (p₀/p_s)^κ | Hot-dry: p₀ = p_s | Hot-dry: moisture per K of T | Moist: p₀ = p_s | Moist: moisture per K of T |
|---|---|---|---|---|---|
| 925 hPa | 1.023 | −1.1% | +1.2% | −1.4% | +0.8% |
| 850 hPa | 1.048 | +4.5% | +9.4% | −3.0% | +1.6% |
| 750 hPa | 1.086 | −3.7% | +4.5% | −4.6% | +3.6% |

Three points follow from this table.

1. The two broken terms push in opposite directions. Referencing the moisture increment to surface T alone raises PFT, because the parcel gets less moisture and condenses higher. Re-gauging p₀ to the surface also shrinks eq 31's Δθ_fc by the same factor, which partly cancels the first effect.
2. The response is not monotonic in p_s. z_fc is quantised by the model levels at which the buoyancy check runs, so a small shift in the saturation point can move the critical β across a level.
3. The magnitude, |1–9|% at 750–925 hPa, exceeds the noise floor and is roughly half of the heat-to-moisture uncertainty. It is also spatially structured: it follows orography. A uniform offset leaves map rankings unchanged, but an orography-correlated one does not.

### 4.3 Moist entropy reference (SQ2)

The PFT chain uses only temperature and pressure along a dry ascent to the saturation point, then a pseudo-adiabat. No diagnostic integrates entropy or uses θ_e or θ_s as a conserved variable. The third-law reference problem (Marquet, 2011; Marquet & Stevens, 2022) therefore does not reach any current pyflam output. It would become relevant if pyflam adopted entropy-based plume budgets in which the fire adds water, because the reference constants then stop cancelling.

### 4.4 Readings A and B (SQ3)

No pyflam quantity changes under either formalism (Sections 3.1–3.2). Their value is conceptual. Reading A shows that representation-dependent curvature quantities make poor diagnostics (Bravetti et al., 2015). Reading B supplies the principle "observables must be invariant under the choice of reference" (Polettini, 2012; Polettini & Vedovato, 2013), and that principle drives the reading C audit.

## 5. Discussion

**What "gauge theory" buys pyflam.** As a modelling framework, nothing. As an *audit discipline*, a cheap and general check: give every quantity a weight, require weight-homogeneous sums and weight-zero outputs. The check found three real covariance breaks in code that has been reviewed repeatedly. None of them shows up on a flat-terrain test column, which is why existing tests would not catch them.

**Is the bias physical?** *Resolved 2026-11-19 from the full text of Tory and Kepert (2021): the gauge is theirs. Table 2 sets P₀ = 10⁵ Pa and κ = 0.2854. Eq 10 takes ΔH = ρ₀c_pdΔθ from dry static energy S_d = c_pdT + gΔz = c_pdθ, referenced to the height of P₀ (appendix B). Eq 30 applies the fire ratio to β_T18 θ_ML, and their manual recipe intersects the θ_ML + 15 K adiabat with q_ML + 1 g/kg. pyflam reproduces this exactly. The terrain dependence is a property of the published method, not of the port.* Original text: The audit establishes that PFT depends on an arbitrary constant. It does not establish which gauge is correct. That depends on how Tory et al. (2018) and Tory and Kepert (2021) defined the heat-to-moisture ratio and the Δθ in eq 31. If their derivation is written in 1000 hPa θ with Australian near-sea-level soundings, the ambiguity hardly mattered for their cases: (1000/980)^0.286 ≈ 1.006. It does matter for ICON-EU columns over the Pyrenees, the pre-Pyrenees and the Apennine ridge. A plausible reading of the physics is that the fire heats near-surface air by ΔT and adds Δq, which favours the "moisture per K of T" form. But this should be checked against the papers before the code changes.

**Recommendations, in order of cost.**

1. **Declare the gauge once.** Add a module constant `P0_PA = 100000.0` and use `_KAPPA_DRY` in `theta_kelvin` instead of the literal 0.286. Pick one dry-air c_p. This is hygiene, with an effect below 0.1%. **Done 2026-11-19:** `_P0_HPA`/`_P0_PA`/`_KAPPA_DRY`, c_p merged at 1004, and the test `test_reference_gauge_is_single_and_consistent`. Measured drift: θ ≤ 0.023 K, PFT +0.26% / −0.05% on the synthetic columns. The §4.2 tables were re-run on the new gauge.
2. **Covariance property test.** *Done 2026-11-19 (`tests/test_atmosphere.py`, `tests/test_gui_core.py`); it would have caught both new bugs above.* Add a pytest that re-gauges p₀ (or shifts a column to p_s = 850 hPa while holding T and q fixed in physical terms) and asserts that weight-0 outputs (z_fc, buoyancy, fire-CAPE without excess) are unchanged to within the β-grid tolerance. It should mark the known breaks as `xfail` until they are resolved, so the next additive-K-into-θ term fails loudly.
3. **Resolve the Tory–Kepert convention.** *Resolved: pyflam matches the paper's P₀ = 10⁵ Pa gauge, and a test pins it. A terrain-consistent variant now exists as `moisture_reference="surface_t"` (opt-in). It was scored on the GRAF capability-margin set (27 columns; `docs/gauge_audit/score_moisture_reference.py`): −0.04 to +0.5% on 24 columns, +5.5 to +15.5% on three, Guissona's margin +0.39 → +0.34, and no change in detections, Fisher p or rank correlations. The validation data cannot tell the two forms apart, so the paper form stays the default.* Read eqs 25–31 of Tory and Kepert (2021) and eqs 29–30 of Tory et al. (2018) for the reference of Δθ and of the 15 K per g/kg ratio. Implement the gauge-fixed variant behind a keyword (default unchanged), and score both on the 24-fire reference set before changing the default.
4. **`fire_cape_grid` excess.** *Resolved: the excess is θ′, so no conversion is needed.* Document whether `theta_excess` is a ΔT or a Δθ. If it is a sonde ΔT, as in `scripts/ingest_inplume_sondes.py`, multiply it by (p₀/p_s)^κ.
5. **Do not adopt readings A or B in the pipeline.** Revisit B only if pyflam builds an explicit entropy-production model of fire–plume coupling.
6. **Naming.** In any publication, call this "reference-state (scale-gauge) covariance" and cite Weyl's origin of the term through Straumann (1996). "Gauge theory of thermodynamics" in a title would promise more than the work delivers.

## 6. Limitations

- **Full texts not read.** For most sources only abstracts and bibliographic records were verified. Tory and Kepert (2021) was subsequently read in full (2026-11-19), which resolved Recommendation 3.
- **Synthetic columns only.** Two idealised soundings were run, shifted in pressure by an approximate RH-preserving transform. No ICON-EU, GFS or GRAF radiosonde columns were used, and no skill score on the 24-fire set was computed, so the Phase 1 validity criterion is met only in its noise-floor form.
- **Level quantisation.** The non-monotonic terrain response comes partly from the level-wise buoyancy check. Its interaction with the gauge effect was not separated out.
- **Inter-model comparator missing.** The planned comparison against the ICON-EU versus GFS spread was not run.
- **Scope.** C-Haines' isobaric-level choice and hybrid-height interpolation are coordinate questions, not gauge questions, and were only noted.

## 7. Conclusion

Of the three senses in which thermodynamics is a gauge theory, only reference-state covariance (reading C) touches pyflam's numbers. It finds three weight mismatches in the PFT, cost-ladder and fire-CAPE chain. On flat terrain they are below the noise. Over 750–925 hPa terrain they produce an orography-correlated PFT bias of about −4% to +9%, roughly half the heat-to-moisture uncertainty. The fix is small: declare the gauge, add a covariance test, and settle the Tory–Kepert reference convention against the source papers. Contact geometry and network gauge thermodynamics are worth knowing as background and are not worth building into the pipeline.

---

## References

Balian, R., & Valentin, P. (2001). Hamiltonian structure of thermodynamics with gauge. *The European Physical Journal B, 21*, 269–282. https://doi.org/10.1007/s100510170202

Bolton, D. (1980). The computation of equivalent potential temperature. *Monthly Weather Review, 108*(7), 1046–1053.

Bravetti, A. (2017). Contact Hamiltonian dynamics: The concept and its use. *Entropy, 19*(10), 535. https://doi.org/10.3390/e19100535

Bravetti, A., López-Monsalvo, C. S., & Nettel, F. (2015). Conformal gauge transformations in thermodynamics. *Entropy, 17*(9), 6150–6168. https://doi.org/10.3390/e17096150

Bryan, G. H. (2008). On the computation of pseudoadiabatic entropy and equivalent potential temperature. *Monthly Weather Review, 136*(12), 5239–5245.

Grmela, M., & Öttinger, H. C. (1997). Dynamics and thermodynamics of complex fluids. I. Development of a general formalism. *Physical Review E, 56*(6), 6620–6632. https://doi.org/10.1103/PhysRevE.56.6620

Hauf, T., & Höller, H. (1987). Entropy and potential temperature. *Journal of the Atmospheric Sciences, 44*(20), 2887–2901.

Marquet, P. (2011). Definition of a moist entropy potential temperature: Application to FIRE-I data flights. *Quarterly Journal of the Royal Meteorological Society, 137*(656), 768–791.

Marquet, P., & Stevens, B. (2022). On moist potential temperatures and their ability to characterize differences in the properties of air parcels. *Journal of the Atmospheric Sciences, 79*(4), 1089–1103. https://doi.org/10.1175/JAS-D-21-0095.1

Mrugała, R., Nulton, J. D., Schön, J. C., & Salamon, P. (1991). Contact structure in thermodynamic theory. *Reports on Mathematical Physics, 29*(1), 109–121. *(Confirmed through a citing reference list only.)*

Polettini, M. (2012). Nonequilibrium thermodynamics as a gauge theory. *EPL, 97*(3), 30003. https://doi.org/10.1209/0295-5075/97/30003

Polettini, M., & Vedovato, F. (2013). *Of dice and men: Subjective priors, gauge invariance, and nonequilibrium thermodynamics* (arXiv:1307.2057). arXiv.

Ruppeiner, G. (1995). Riemannian geometry in thermodynamic fluctuation theory. *Reviews of Modern Physics, 67*(3), 605–659.

Schnakenberg, J. (1976). Network theory of microscopic and macroscopic behavior of master equation systems. *Reviews of Modern Physics, 48*(4), 571–585. https://doi.org/10.1103/RevModPhys.48.571

Straumann, N. (1996). *Early history of gauge theories and weak interactions* (arXiv:hep-ph/9609230). arXiv.

Tory, K. J., & Kepert, J. D. (2021). Pyrocumulonimbus firepower threshold: Assessing the atmospheric potential for pyroCb. *Weather and Forecasting.* https://doi.org/10.1175/WAF-D-20-0027.1 *(Published online December 2020. Volume and pages not verified.)*

Tory, K. J., Thurston, W., & Kepert, J. D. (2018). Thermodynamics of pyrocumulus: A conceptual study. *Monthly Weather Review, 146*(8), 2579–2598. https://doi.org/10.1175/MWR-D-17-0377.1

*Dropped by the verification gate:* Weyl (1918) is cited only through Straumann (1996), because its Sitzungsberichte pagination could not be confirmed. Öttinger and Grmela (1997, Part II) and Quevedo's geometrothermodynamics were left out because they are not needed for any claim.

---

## Appendix A. Review panel (Phases 1, 3, 5) and revisions (Phase 6)

**Devil's Advocate, Checkpoint 1 (scoping): REVISE → resolved.** Readings A and B risked being metaphor with no output. They were resolved by restricting them to a short theory section and making reading C carry a falsifiable audit.

**Devil's Advocate, Checkpoint 2 (analysis): PASS with conditions.**
- *Cherry-picking:* only two idealised columns. → Stated in §6. Every perturbation run is reported, including the ones that cancel or flip sign.
- *Alternative explanation:* the non-monotonic terrain response could be level quantisation rather than the gauge. → Acknowledged in §4.2 and §6. The weight analysis (§4.1) shows the dependence exists independently of the numbers.
- *Confirmation bias:* "gauge theory" could be oversold. → The report concludes that A and B should *not* enter the pipeline.

**Editor-in-Chief: MINOR REVISION.** (i) The original title used "gauge theory" unqualified; it was retitled to "gauge freedom" and "reference-state covariance". (ii) The claim that the Tory–Kepert convention is ambiguous needed hedging because the full text was not read; this is now marked unverified in §4.1, §5 and §6. (iii) The table header in §4.2 needed a stated baseline; added.

**Ethics review: CLEARED.** No human subjects, no dual-use concern (fire-danger forecasting is protective). Attribution: every claim about the literature is tied to a verified record, and secondhand confirmations are labelled. AI disclosure below.

**Devil's Advocate, Checkpoint 3 (final): PASS.** *So what?* The uniform part changes nothing. The orography-correlated part is the only result that matters, and it is reported with its size relative to an uncertainty pyflam already carries. *Strongest counter-argument:* NWP temperature and moisture errors over complex terrain probably exceed 5% in PFT anyway. That is accepted, but those errors are random across runs while the gauge bias is systematic, and a systematic error costs nothing to remove.

## Appendix B. Reproducibility

```bash
# from the repository root, commit f93b693
python docs/gauge_audit/audit.py     # global-gauge and comparator runs -> audit_output.txt
python docs/gauge_audit/terrain.py   # terrain sweep 1000/925/850/750 hPa -> terrain_output.txt
```

The scripts monkeypatch `theta_kelvin`, `_mixed_layer_state` and `_saturation_point` with gauge-parameterised copies. With the default arguments they reproduce the shipped functions (baseline rows equal the as-coded PFT). No library code is modified.

## AI disclosure

This report was produced with an AI research assistant (Claude, Anthropic) running the academic-research-skills `deep-research` pipeline: literature search and bibliographic verification by web search, static code analysis, sensitivity runs, drafting, and simulated review checkpoints. All references were checked against at least one independent bibliographic record. The author should read the Tory and Kepert (2021) and Tory et al. (2018) full texts before acting on Recommendation 3.
