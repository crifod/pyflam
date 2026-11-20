# Pyroconvection product: ICON-2I vs hybrid (ICON-EU) — 2026-07-15

A cell-by-cell comparison of the two atmospheric sources in `tests/pyroconv_daily.py`, run on
the same day, same grid, same fuel gate, with the **current shipped defaults** (model-level
hybrid; the shear diagnostic on the ICON-EU levels; `rh_top_moist = 60`):

- **`icon2i`** — atmosphere from ICON-2I 2.2 km, 5 pressure levels. Mixed-layer stability is a
  **mixing-depth proxy** (unmeasurable on 5 levels); no shear (4-diagnostic ladder).
- **`hybrid`** (default) — atmosphere from ICON-EU **native model levels** (regridded to the
  2.2 km grid), with a **measured** mixed-layer dθ/dz and the **shear** diagnostic (5-diagnostic
  ladder). Surface fields and the 10 MW/m fuel gate stay ICON-2I 2.2 km.

Run: 2026-07-15 00Z, valid same day, 3-hourly, Tuscany, 9,907 land cells. Both products and
their diagnostic rasters are committed; the numbers are reproducible via `compare_sources.py`.

> **Summary.** Operationally the two are almost identical once the fuel gate is applied (a
> low-danger day). Atmospherically they now *both* reach the pyroCb classes — but by different
> routes: `icon2i`'s **deeper, proxy-driven ABL** pushes LCL/ABL below 1 and sends cells into
> the **deep** branch, while `hybrid`'s **measured, shallower ABL** (LCL/ABL ≈ 1.2, matching
> ERA5 and inland soundings) sends them into **overshooting**. The independent evidence
> (`scripts/validation/`) places the truth with the hybrid's LCL/ABL.

---

## 1. Potential (atmospheric upper bound) — class distribution

Percent of land cells, per valid hour. 0 surface · 1 convection · 2 overshoot · 3 resilient ·
4 deep pyroCu/pyroCb.

**ICON-2I (proxy, 4-diagnostic):**

| Hour | surface | convect | overshoot | resilient | deep |
|:--|--:|--:|--:|--:|--:|
| 09Z | 9.0 | 7.1 | 53.8 | 10.8 | **19.2** |
| 12Z | 5.5 | 1.4 | 79.8 | 0.5 | **12.8** |
| 15Z | 21.1 | 8.6 | 59.8 | 0.3 | **10.2** |
| 18Z | 98.3 | 1.6 | 0.2 | 0.0 | 0.0 |

**Hybrid (measured, 5-diagnostic):**

| Hour | surface | convect | overshoot | resilient | deep |
|:--|--:|--:|--:|--:|--:|
| 09Z | 4.0 | 44.3 | 50.5 | 1.2 | 0.0 |
| 12Z | 18.5 | 22.0 | 58.5 | 0.4 | **0.7** |
| 15Z | 15.5 | 22.4 | 53.3 | 1.6 | **7.1** |
| 18Z | 100.0 | 0.0 | 0.0 | 0.0 | 0.0 |

Both share the diurnal envelope (nothing overnight, activity 09–15Z, collapse by 18Z) and
both now produce the resilient/deep tail (the `rh_top_moist = 60` relaxation admits pyroCb in
the dry regime where it occurs). The differences:

- **Convection plume** is large in the hybrid (22–44%) and near-absent in icon2i (1–9%): the
  hybrid finds many cells marginally stable or with the LCL above the ABL top (convection
  plume), where the proxy calls them pyroCu.
- **Deep pyroCu/pyroCb** is *higher in icon2i* at midday (12.8% vs 0.7% at 12Z), the reverse of
  overshooting. This is the key point, explained in §4.

## 2. Agreement & 12Z confusion (potential)

| Hour | 00 | 03 | 06 | 09 | 12 | 15 | 18 | 21 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| agreement % | 100 | 100 | 99.4 | **31.8** | 53.3 | 66.8 | 98.3 | 100 |

12Z confusion, rows = ICON-2I class, cols = hybrid class, % of land:

| ICON-2I ↓ / hybrid → | surface | convect | overshoot | resilient | deep |
|:--|--:|--:|--:|--:|--:|
| overshoot | 12.0 | 19.0 | **48.7** | 0.0 | 0.1 |
| **deep** | 3.5 | 0.5 | **8.0** | 0.3 | 0.5 |

Most of icon2i's overshooting is shared (48.7 pts) or drops a rung in the hybrid (12 → surface,
19 → convection). Almost all of icon2i's **deep** (12.8%) becomes **overshooting** in the hybrid
(8.0 pts) — the two disagree on whether these inland afternoon columns are deep-pyroCb-capable or
overshooting.

## 3. Fuel-gated (expected) product — nearly identical

Percent of land cells reaching classes 3–4 (resilient + deep) in the operational, fuel-gated map:

| Hour | ICON-2I | hybrid |
|:--|--:|--:|
| 12Z | 0.08 | 0.02 |
| 15Z | 0.14 | 0.07 |

On this low-danger day the 10 MW/m fuel gate leaves ≥97% of the domain at surface plume in both,
so the large §1 disagreement is almost entirely in cells where no fire could realise it. The two
products are operationally interchangeable today; the divergence would matter on a high-danger day.

## 4. Why they differ — the diagnostics

Land-median at the active hours, ICON-2I (2I) vs hybrid (HY):

| Diagnostic | 12Z 2I / HY | 15Z 2I / HY |
|:--|:--|:--|
| Rib ABL (m) | **2311 / 1510** | 1734 / 1429 |
| **LCL / ABL** | **0.79 / 1.23** | 0.94 / 1.29 |
| ML dθ/dz (K/m) | 3e-4 / 1.4e-3 | 4e-4 / 1.0e-3 |
| RH at ABL top (%) | 52.8 / 45.9 | 54.1 / 46.9 |

The whole story is **LCL/ABL**. icon2i's ABL is ~800 m deeper (its proxy over-mixes), so LCL/ABL
falls **below 1** → cells enter the resilient/deep branch → icon2i's high deep fraction. The
hybrid's measured, shallower ABL gives LCL/ABL **≈ 1.2** → the overshooting branch. Same weather,
opposite class, decided by the ABL depth.

**Which is right?** The independent evidence (`scripts/validation/README.md`) is unambiguous:
inland Madrid radiosondes and ERA5 over Catalonia both put the real inland-summer LCL/ABL at
**≈ 1.0–1.2** with a deep ABL — the **hybrid's** value, not icon2i's 0.79. So the hybrid's
overshooting call is the physically correct one, and icon2i's deep-heavy midday is an artifact of
its over-deep proxy ABL. ERA5's own class balance (~58% overshoot) matches the hybrid, not icon2i.

## 5. Limitations

- The `rh_top_moist = 60` gate and the 1.1e-3 K/m stability threshold sit inside their validated
  error bars, so class *counts* are indicative, not calibrated.
- The hybrid trades horizontal resolution (6.5 km atmosphere vs 2.2 km) for vertical resolution;
  the fuel gate stays 2.2 km, so fire-power detail is preserved.
- Neither is validated against *observed* pyroCu in Tuscany; the diagnostics are validated
  against radiosondes and ERA5 (see `scripts/validation/`).

## 6. Which to use

- **Default is `hybrid`** — its ABL, mixed-layer stability and LCL/ABL are the measured,
  independently-validated ones; icon2i's proxy ABL mislocates the overshoot/deep boundary.
- icon2i remains the automatic fallback when the ICON-EU run is unpublished, and is 2.2 km
  throughout.
- Read the gated map for decisions and the diagnostic rasters (ABL, LCL/ABL, ML dθ/dz) to
  understand a class; do not read potential-map class totals quantitatively from either.

## 7. Reproduction

```bash
PYTHONPATH=src python tests/pyroconv_daily.py 2026-07-15 0                    # hybrid (default)
PYROCONV_SOURCE=icon2i PYTHONPATH=src python tests/pyroconv_daily.py 2026-07-15 0
python scripts/validation/compare_sources.py 2026-07-15                       # this report's tables
```

Products: `pyroconv_tuscany_{icon2i,hybrid}_{potential,gated}_2026-07-15.png`, the two PDFs, and
the per-hour diagnostic GeoTIFFs under `rasters_{icon2i,hybrid}_2026-07-15/` (gitignored). Method
and validation: `scripts/validation/README.md`.
