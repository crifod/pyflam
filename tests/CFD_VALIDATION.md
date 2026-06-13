# Validating the RANS wind solver (`pyflam.cfd`)

The pure-Python parts (mesh topology, case generation, reader, discovery) are
covered by `test_cfd_*.py` without OpenFOAM. The CFD *physics* needs OpenFOAM and
is exercised by `test_cfd_integration.py` (auto-skipped when OpenFOAM is absent).
This file describes the deeper, run-on-demand validation against canonical
benchmarks — the real evidence the solver is trustworthy, beyond "it runs."

## Prerequisites

ESI OpenFOAM (openfoam.com), invoked via the `openfoam` wrapper. On Apple
Silicon: `brew install gerlero/openfoam/openfoam`. Confirm:

```bash
openfoam -c 'command -v simpleFoam buoyantBoussinesqSimpleFoam checkMesh'
```

pyflam discovers it through `PYFLAM_OPENFOAM` or `openfoam` on `PATH`.

## 1. Mesh quality — `checkMesh`

Generate a case and check the terrain mesh has no errors:

```python
from pyflam.cfd import mesh, case
m = mesh.build_terrain_mesh(dem, cellsize, nz=24, domain_height=2000)
case.write_case("/tmp/case", m, case.CaseConfig(speed=10, direction=270))
```
```bash
openfoam -c "cd /tmp/case && checkMesh -constant"   # expect: "Mesh OK."
```
Watch max non-orthogonality (< ~65) and skewness (< ~4) on steep terrain; raise
`nz` / `domain_height` or lower `expansion_ratio` if they climb.

## 2. Richards & Hoxey horizontal homogeneity (neutral)

Over **flat** ground the inlet ABL profile must be preserved across the domain
(no streamwise decay). `solve_rans` on a flat DEM should return a near-uniform
field at the chosen height with direction unchanged (the automated
`test_neutral_flat_preserves_wind` is a coarse version). For a full check, sample
vertical U/k/epsilon profiles at several x-stations and confirm they overlay the
inlet profile.

## 3. Askervein Hill (neutral terrain benchmark)

The standard real-terrain CFD benchmark. Run `solve_rans` on the Askervein DEM
with the measured upwind reference wind and compare the fractional speed-up
(ΔS) along lines A/AA and at hilltop HT against the field campaign. Expect
hilltop speed-up in the right ballpark (~+80% at 10 m); RANS k-ε typically
under-predicts the lee wake.

## 4. Diurnal slope flows (non-neutral buoyancy)

With `buoyant=True` and a `surface_heat_flux` (W/m^2):
- **Daytime / convective** (`q > 0`): upslope (anabatic) flow strengthens; on a
  ridge the near-surface flow tilts upslope and the crest speed-up grows.
- **Nighttime / stable** (`q < 0`): drainage (katabatic) downslope flow forms;
  with weak ambient wind, near-surface vectors point down the slope.

Compare against ESI tutorials `verificationAndValidation/atmosphericModels/
atmForestStability` (stable / unstable setups) for profile shapes. Use weak
ambient wind (1-2 m/s) so buoyancy dominates; expect to lower relaxation /
raise iterations for convergence (buoyant RANS over terrain is fragile — see the
plan's risk notes).

## 5. Cross-checks

- **vs. the diagnostic solver**: for a neutral case, `pyflam.cfd.solve_rans` and
  `pyflam.windsolver.solve_mass_consistent` should agree qualitatively (ridge
  speed-up, channeling). RANS adds lee separation/wake the mass-consistent model
  cannot produce — that difference is expected and is the reason to use RANS.
- **End-to-end**: `solve_rans` -> `WindField.to_midflame` ->
  `basic_fire_behavior` should yield faster spread on accelerated (crest/gap)
  cells.

## Notes / known limitations

- Staircase terrain (structured hex with terrain-following bottom); very steep
  cells raise non-orthogonality.
- One representative `z0`; per-cell roughness from fuels is a future refinement.
- Diurnal is modelled as a *quasi-steady* heated/cooled snapshot (fixed surface
  flux), not a transient diurnal cycle.
- Single `z0`, neutral inlet ABL; the inlet does not yet carry a stability-scaled
  (Monin-Obukhov) profile.
