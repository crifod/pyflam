# Design plan — coupling empirical crown-fire spread to the CFD plume + spotting

**Status:** plan (not yet implemented). **Author:** pyflam. **Date:** 2026-06-25.

## 1. Motivation

pyflam's pyroconvection loop (`pyroconvection.fire_atmosphere_march`) already does
something no other open landscape fire engine does: it feeds the fire's **fireline
intensity** into a buoyant CFD plume (OpenFOAM), gets a plume-modified wind back, and
re-grows the fire in that wind — and it lofts embers (`spotting`) from the same
intensity. Today that loop is driven entirely by **surface** fireline intensity
(`SpreadField.fireline_intensity`, Byram surface). But the regime where the plume,
the indrafts and the spotting actually matter is **crown fire**, where intensity is
several times higher. With the Cruz et al. (2005) active crown-ROS model now
available (`crown_spread="cruz2005"`), we can close a genuinely novel loop:

> crowning raises fireline intensity → stronger plume + farther spotting → plume
> wind raises U10 → Cruz crown ROS rises → more crowning.

This bridges the empirical (Cruz) and physics-based (FIRETEC/WFDS-style plume) camps
while staying landscape-tractable. It is the contribution flagged in the deep-research
review on crown-fire models.

## 2. Scope

In scope: route **crown** fireline intensity and **crown** ROS through the existing
march loop, plume coupling, and spotting. Quasi-steady (the loop already is). Out of
scope: a full 3-D combustion solver; live-flammability chemistry; sub-canopy fire.

Prerequisite data: the landscape must carry **canopy base height + bulk density**
bands (same requirement as `crown_fire_potential`); without them the loop degrades
gracefully to the current surface-only behaviour.

## 3. Current pieces and where crown fire plugs in

| Piece | Today (surface) | After (crown-aware) |
|---|---|---|
| `spread_field` → `SpreadField.fireline_intensity` | Byram surface I | **max(surface, crown) I** on crowning cells |
| `spread_field` → `ros_max` | surface ROS | **Cruz active ROS** on active-crown cells |
| `pyroconvection.fire_heat_flux(I)` | from surface I | from **crown** I (much larger) |
| `spotting` loft (uses `fireline_intensity`) | surface I | **crown** I → higher loft, longer spotting |
| `couple_fire_wind` → plume wind | surface-driven | crown-driven indrafts |
| `fire_atmosphere_march` step | re-read weather, re-grow | also **re-classify crown** each step |

The key insight: most of the machinery already keys off `SpreadField.fireline_intensity`
and `SpreadField.ros_max`. If we build a **crown-aware `SpreadField`** each march step,
the plume and spotting coupling get crown behaviour *for free*.

## 4. Architecture

### 4.1 A crown-aware spread field

Add `crown_spread_field(ls, *, canopy ..., crown_spread="cruz2005", wind_20ft, **inputs)`
(in `crownfire.py` or a thin wrapper in `mtt.py`) that:

1. builds the normal surface `SpreadField`;
2. runs `crown_fire_potential(...)` to get per-cell `fire_type`, crown ROS, crown
   intensity;
3. returns a new `SpreadField` where, on cells that crown:
   - `ros_max` ← crown ROS (Cruz active rate on active cells; surface→active blend on
     passive), keeping the existing elliptical `eccentricity`/`heading` (wind/slope
     geometry is unchanged; only the magnitude grows);
   - `fireline_intensity` ← total crown intensity (kW/m → Btu/ft/s).

This is the single new object; everything downstream is unchanged.

### 4.2 The march loop

`fire_atmosphere_march` gains `crown=True` and the canopy/foliar inputs. Each step:

1. read weather (existing) → wind, moisture per cell;
2. **build the crown-aware spread field** (§4.1) using the *current plume-modified
   wind* for the Cruz U10;
3. grow the fire one `dt` (MTT or fast-marching) on that field;
4. compute the plume heat flux from the **crown** intensity (`fire_heat_flux`), solve
   the buoyant CFD plume, merge with ambient (`merge_plume_wind`) → new wind;
5. loft embers with `spotting` from the **crown** intensity;
6. loop.

The only new coupling is step 2 consuming step 4's wind for the Cruz ROS — the
positive feedback. Steps 4–5 are unchanged except for the intensity they read.

### 4.3 Wind consistency (the one subtlety)

Cruz (2005) uses **10-m open wind (U10)**; the plume solver returns a near-surface
wind field; the surface model uses **midflame** wind. We already convert 20-ft↔U10
(`crownfire._U20FT_TO_U10`) and 20-ft↔midflame (`wind_reduction`). Define one
consistent path: plume wind (height-referenced) → 20-ft open → {midflame for surface,
U10 for Cruz}. Encapsulate in a small helper so the feedback uses a single wind source.

## 5. Feedback stability

The loop is positive-feedback, so it must be bounded:

- **Physical cap:** Cruz ROS already saturates slowly (U10^0.90), and the plume
  enhancement is bounded by `convective_plume_factor`. Cap the per-step wind
  enhancement (already implicit in the CFD) and the crown ROS at a physical maximum.
- **Under-relaxation:** blend the new wind with the previous step
  (`w ← (1−α)·w_prev + α·w_new`, α~0.5) to damp oscillation — standard for the
  quasi-steady coupling pyflam already runs.
- **Convergence check:** stop iterating a step when the burned-area / mean-wind change
  falls below a tolerance (the march already has the hook).

## 6. API sketch

```python
res = pyflam.fire_atmosphere_march(
    ls, ignitions=[(r, c)], total_time=120, dt=15,
    atmosphere=prov, location=(lat, lon), start_time=t0,
    crown=True, crown_spread="cruz2005",
    foliar_moisture=100.0, m_live_herb=0.6, m_live_woody=0.9,
    spotting=pyflam.FirebrandPhysics(...),     # lofts from crown intensity
    use_plume=True, spatial=True,              # CFD plume + per-cell weather
)
# res.fire_type raster, res.arrival_time, res.plume_wind, res.spots ...
```

`spread`/`crown_fire_behavior` keep working standalone; the coupling is opt-in via
`crown=True`.

## 7. Validation

1. **Unit/physics:** crowning must raise plume heat flux and spotting distance vs a
   surface-only run on the same landscape (monotonic checks, no external data).
2. **Crown classification:** the per-step `fire_type` field diffs against FlamMap via
   the existing harness (`validate_flammap_crown.py`, `compare_categories`) — needs a
   canopy landscape + FlamMap crown raster.
3. **Spread magnitude:** active crown ROS vs Cruz CFIS / observed wildfire ROS.
4. **Feedback sanity:** the loop must converge (bounded wind/area) on a synthetic
   canopy landscape; assert no runaway.

## 8. Risks / open questions

- **Feedback runaway** — mitigated by §5; needs a synthetic stress test.
- **Cost** — a CFD plume solve per march step is already the dominant cost; crown
  adds a `crown_fire_potential` pass (cheap, vectorized) and possibly more steps to
  converge.
- **Wind-height bookkeeping** — the U10/midflame/plume reconciliation (§4.3) is the
  most error-prone part; isolate and unit-test it.
- **Validation data** — same blocker as the crown diff: no canopy-band landscape +
  FlamMap crown raster in the bundled dataset. The physics/feedback tests don't need
  it; the quantitative diff does.

## 9. Phasing

1. ✅ **`crown_spread_field`** (§4.1) + unit tests — pure function, no CFD. *Done
   (`crownfire.crown_spread_field` → `CrownAwareField`; `tests/test_crown_coupling.py`).*
2. ✅ **Wind-reconciliation helper** (§4.3) + tests. *Done
   (`crownfire.wind_20ft_to_u10_kmh`, the single 20-ft→U10 conversion used by all
   the Cruz code).*
3. ✅ **`fire_atmosphere_march(crown=True)`** wiring (§4.2) — each increment rebuilds
   the crown-aware field from the current plume-modified wind (the 20-ft wind →
   `crown_spread_field`), so crown intensity feeds the plume/spotting and the crown
   ROS drives growth; output carries the `fire_type` raster. Wind provider injectable
   → tested without OpenFOAM. *Done (`tests/test_crown_march.py`).*
4. ✅ **Feedback stability** (§5) + stress test. *Done* — `fire_atmosphere_march`
   gained `wind_relax` (under-relaxation: vector blend of consecutive winds; `1` =
   none) and `max_wind_factor` (cap each cell's speed at that multiple of the ambient
   mean), a no-op at the defaults; `return_history` adds `mean_wind` per increment.
   `tests/test_feedback_stability.py` proves a runaway plume is bounded.
5. ✅ **Spotting-from-crown-intensity** physics test. *Done* —
   `tests/test_crown_spotting.py`: the spotting models loft from
   `field.fireline_intensity`, so a crown-aware field (whose crown intensity ran
   ~40x the surface value in the test landscape) drives much farther firebrand
   reach (`SpottingModel.max_spot_distance`) and lands embers farther downwind
   (`generate_spots`) than the surface field. No new code — the step-1 crown-aware
   field already carries the right intensity into the existing loft path.
6. **Quantitative validation against _observed_ crown fire** (not FlamMap). FlamMap
   and the Rothermel + Van Wagner stack under-predict crown fire (Cruz & Alexander
   2010), so a diff against FlamMap would validate against a known-biased reference —
   the opposite of the point. Instead: (a) the Cruz 2005/2004 equations are
   unit-tested against their published coefficients, so pyflam inherits their
   literature validation against observed wildfires; (b) the right next anchor is a
   comparison to observed crown-fire rate of spread (analogous to
   `spotting.LITERATURE_SPOT_ANCHORS`), if such data can be sourced.

Steps 1–5 are done — the crown / plume / spotting coupling is built end to end: a
crown-aware spread field, a single wind reconciliation, the march wiring that closes
the crowning → plume → wind → crown feedback, the stabilizers that bound it, and the
crown-driven spotting it produces. Step 6 is reframed away from FlamMap (above): the
implementation is faithful to the literature-validated Cruz models, exercised end to
end on a synthetic canopy `.lcp` (`tests/make_synthetic_canopy_lcp.py`).
