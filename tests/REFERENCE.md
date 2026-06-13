# Validating pyflam against real FlamMap / BehavePlus

The golden-master numbers in `test_rothermel.py` only lock *our* current output.
They are **not** independently verified. The real acceptance test is to diff
pyflam against the tools FlamMap shares its surface-fire engine with. Do this on
a Windows machine (or VM) where FlamMap/BehavePlus run.

## Option A — BehavePlus (fastest to set up)

BehavePlus uses the same Rothermel surface-fire core as FlamMap and exposes every
input directly, so it's the cleanest reference for this step.

1. Install BehavePlus (USFS, Windows).
2. New run, module **SURFACE**.
3. For each standard fuel model, enter the scenario used in the golden table:
   - 1-h / 10-h / 100-h dead moisture = 6 / 7 / 8 %
   - live herb / live woody moisture = 60 / 90 %
   - midflame wind speed = 5 mi/h
   - slope = 0 %
   - wind/slope directed upslope (head fire)
4. Record **rate of spread**, **flame length**, **fireline intensity**,
   **reaction intensity**.
5. Paste the values into the table below and tighten the tolerance in
   `test_golden_master_ros` (and add tests for FL / I_B / I_R).

## Option C — landscape ROS raster vs a real FlamMap run (automated)

`pyflam.validate` + `tests/validate_flammap_ros.py` diff pyflam's max-spread ROS
against a FlamMap `ROS.tif` cell by cell. FlamMap's "Rate of Spread" is the
head-fire (maximum) rate combining wind and slope, so the matching quantity is
`pyflam.spread_field(...).ros_max` (the wind+slope **vector** combination), not
the scalar `basic_fire_behavior`.

```bash
# sweep wind speeds to identify the run's wind, then compare at the best fit:
PYTHONPATH=src python tests/validate_flammap_ros.py --scan 8,10,12,14,16,20
# compare at a known wind and write a (pyflam - FlamMap) diff raster:
PYTHONPATH=src python tests/validate_flammap_ros.py --wind-mph 14 --wind-dir 45 \
    --diff-out diff.tif
```

The wind scan is a **diagnostic** to recover the (unparsed) project wind, kept
distinct from validation. `compare_fields` reports bias / MAE / RMSE, median
ratio, Pearson r, an OLS fit, log-space RMSE (ROS spans decades), "within X%"
fractions, and — importantly — the burn/no-burn **classification** agreement,
which catches nonburnable / moisture-of-extinction mismatches a numbers-only diff
would miss.

### Result — Tuscany 100 m landscape (2470×2166, 5.35M cells)

FlamMap run: `Toscana.tif` (FM40), moisture 6/7/8/60/90, `ROS.tif` in chains/hr.
The crown bands are zero, so this validates **surface ROS only**, and (height 0)
the wind adjustment factor is the unsheltered, fuel-depth-based one. Best-fit
wind ~16 mi/h, dir 45° (consistent with the `15.000000 ... 45` in the `.fmp`).

The first pass used a flat 0.4 WAF; switching to the **per-fuel-depth WAF**
(`pyflam.wind_reduction`, Albini & Baughman 1979) closes nearly the whole gap:

| metric | flat 0.4 WAF | per-fuel WAF (`--waf auto`) |
|---|---:|---:|
| burn/no-burn agreement | 100.00% | **100.00%** |
| Pearson r (1.64M burning cells) | 0.991 | **0.9998** |
| median ratio pyflam/FlamMap | 1.06 | **0.982** |
| RMSE | 12.4 ft/min | **1.21 ft/min** |
| within 10% / 25% | 27% / 90% | **98.4% / 100%** |
| log10-ratio RMSE | 0.088 (×1.22) | **0.0115 (×1.03)** |
| OLS slope / intercept | 0.74 / 2.2 | **0.981 / 0.04** |

pyflam reproduces FlamMap's surface ROS to **~3% typical** across the whole
1.64M-cell burning landscape — slope 0.98, r 0.9998. This both validates the
Rothermel implementation and confirms the per-fuel wind adjustment factor was the
missing piece. The small residual (~2% low) is within FlamMap's own ROS output
rounding and minor directional / length-to-breadth differences.

## Option D — spread/perimeter (MTT) engine vs a FlamMap MTT run

`tests/validate_flammap_mtt.py` checks the spread engine against a FlamMap "MTT
Random Ignitions" burn-probability run (`burn_prob_test_*`). Its `run_log.txt`
pins everything: 2000 random-ignition fires, 60-min sim time, wind 30 km/h @ 10 m
from 45°, spot probability 0.1.

```bash
PYTHONPATH=src python tests/validate_flammap_mtt.py --direction   # deterministic
PYTHONPATH=src python tests/validate_flammap_mtt.py --burnprob --nfires 2000
```

**Max spread direction (`MAX_SPRE_DIR.tif`, radians) — clean, deterministic.**
pyflam's `spread_field(...).heading` (the wind+slope vector combination) vs
FlamMap, over 1.64M cells:

| metric | value |
|---|---|
| mean \|angular error\| | **0.96°** |
| median \|angular error\| | 0.45° |
| RMSE | 1.68° |
| within 5° / 10° | **97.5% / 99.9%** |

Essentially exact. Together with the ROS match (Option C), both deterministic
inputs to fire growth — per-cell spread *rate* and *direction* — are validated, so
the arrival-time engine rests on validated components.

**Burn probability (`BURN_PROB.tif`) — partly recoverable, with spotting.** That
FlamMap run had **ember spotting on** (spot prob 0.1), and the landscape is mostly
slow fuel (median ROS ~10 ft/min) threaded with fast-grass corridors. Without
spotting, pyflam fires choke at the slow-fuel barriers (a single 60-min fire from
the fastest cell burns only ~99 cells — correct no-spotting behaviour) and the
mean burn probability is ~100× below FlamMap.

`pyflam.spotting` + `spread_with_spotting` model firebrand loft-and-drift so fires
jump barriers (`--spotting`, plus `--physics` for the stochastic model). Two
spotting models:

* `SpottingModel` (parameterized): with tuned coefficients and a +30% load factor
  (`--load-factor 1.3`) mean burn probability rises ~36×, to within ~3× of FlamMap.
* `FirebrandPhysics` (stochastic, physics-based): spotting emerges from plume
  buoyancy and firebrand mechanics, so it is physically *bounded* — with sensible
  physical parameters it raises mean burn probability ~8× (it can't be cranked
  arbitrarily the way the parametric model can). This is the more faithful model
  and the recommended one.

Either way the cell-level *pattern* stays uncorrelated, for reasons that are not
engine error: (1) at a few hundred fires with a different random seed the
burn-probability field is Monte-Carlo noise (FlamMap ran 2000), and (2) neither
pyflam spotting model is FlamMap's exact one. So burn probability is a mechanism
demonstration, not a clean validation. The deterministic checks above (ROS,
direction) remain the rigorous result; a FlamMap single-fire arrival-time /
perimeter export with spotting **off** would give a clean growth diff.

## Option E — firebrand spotting vs measured spot distances

The physics-based `FirebrandPhysics` model is calibrated and checked against
spot-distance data independently of FlamMap. `spot_distance_distribution(I, U)`
returns the per-source landing distribution; `calibrate_front_length(anchors)`
fits the one weakly-constrained physical length to measured `(intensity, wind,
max-distance)` rows (iterated, because longer loft also burns more brands out —
a self-limiting nonlinearity); `spot_distance_report` tabulates modelled vs
measured.

Calibrated (`front_length = 115 m`) against literature anchors (Albini 1979/1983;
Koo et al. 2010, IJWF; Sardoy et al. 2008), p98 of the modelled distribution:

| fire (intensity, wind) | modelled max | literature | ratio |
|---|---:|---:|---:|
| surface (100 kW/m, 5 m/s) | ~95 m | 100 m | ~0.9 |
| shrub (1000 kW/m, 8 m/s) | ~460 m | 500 m | ~0.9 |
| timber (5000 kW/m, 10 m/s) | ~1300 m | 1500 m | ~0.85 |
| crown/extreme (30000 kW/m, 15 m/s) | ~4600 m | 3000 m | ~1.5 |

The three non-extreme classes land within ~15%; crown spotting overshoots, which
is acceptable (extreme spotting is highly variable and routinely exceeds a few km).
The residual is an intensity-scaling shape difference a single length can't remove
— refine the firebrand size distribution / burnout against real data to close it.
**To recalibrate on your own measurements**, pass a list of `(intensity_kW_m,
wind_m_s, distance_m)` rows to `calibrate_front_length` and set
`FirebrandPhysics(front_length=...)`.

## Option F — MTT growth vs a FlamMap single-fire perimeter (spotting off)

The clean growth-engine acceptance test (Option D's burn-probability run was
spotting-confounded). Produce a **FlamMap MTT single fire with spotting OFF** and
export its time-of-arrival raster, then:

```bash
PYTHONPATH=src python tests/validate_flammap_perimeter.py \
    --toa flammap_arrival.tif --toa-units minutes \
    --ignition-xy <x> <y> --wind-mph 16 --wind-dir 45 --perimeter-time 120
```

`pyflam.validate.compare_perimeters` reports burned-area **Jaccard / Dice**, the
area ratio and the symmetric **mean / Hausdorff perimeter distance** (map units);
`compare_arrival_times` reports timing agreement on the cells both fires reach.
The current Tuscany dataset has only spotting-*on* burn-probability runs, so this
awaits a spotting-off export; the metrics are unit-tested on synthetic perimeters.

## Option B — FlamMap "Basic Fire Behavior"

1. Build (or load) a small landscape (`.lcp`) of uniform fuel model and slope.
2. Set constant fuel moistures and a constant wind in the run settings; note
   FlamMap applies a **wind reduction factor** to get midflame wind from the
   20-ft wind — either set the reduction factor to 1.0 or back it out so the
   midflame wind matches what you pass to `pyflam.spread`.
3. Run Basic Fire Behavior; export the ROS / FL / intensity rasters.
4. Compare a uniform-cell value against `pyflam` for the same fuel + inputs.

Watch out for: midflame vs. 20-ft wind, slope as percent vs. rise/run vs.
degrees, and English vs. metric output units.

## Reference table (fill in)

Scenario: dead 6/7/8 %, live 60/90 %, midflame wind 5 mi/h, slope 0.

| FM | pyflam ROS (ft/min) | BehavePlus ROS | pyflam FL (ft) | BehavePlus FL | within 5%? |
|----|--------------------:|---------------:|---------------:|--------------:|:----------:|
| 1  | 103.28 |  |  |  |  |
| 2  | 40.22  |  |  |  |  |
| 3  | 129.56 |  |  |  |  |
| 4  | 89.72  |  |  |  |  |
| 5  | 29.05  |  |  |  |  |
| 6  | 37.23  |  |  |  |  |
| 7  | 32.76  |  |  |  |  |
| 8  | 2.23   |  |  |  |  |
| 9  | 9.65   |  |  |  |  |
| 10 | 10.03  |  |  |  |  |
| 11 | 6.73   |  |  |  |  |
| 12 | 14.63  |  |  |  |  |
| 13 | 17.66  |  |  |  |  |

If any row is off by more than a few percent, the likely culprits (in order):
size-class net-load weighting, live moisture of extinction, or the wind factor
exponents — all isolated in `pyflam/rothermel.py` with equation references.
