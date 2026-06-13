"""Validating pyflam against real FlamMap raster output (roadmap step 6).

The golden-master tests only lock pyflam's own numbers. The real acceptance test
is a cell-by-cell diff against a FlamMap run on the same landscape and inputs.
This module provides the generic machinery for that diff — it is deliberately
data-agnostic (it compares two aligned arrays), so it works for any FlamMap
output band, not just one dataset. The data wiring for a specific landscape lives
in a script (see ``tests/validate_flammap_ros.py``).

What it does:

* :func:`compare_fields` — robust comparison metrics between a pyflam field and a
  reference field over a shared mask (bias, MAE, RMSE, ratios, correlation, an
  OLS fit, and — because fire ROS spans orders of magnitude — log-space stats and
  "within X%" fractions). Also reports the burn/no-burn *classification*
  agreement, which catches moisture-of-extinction and nonburnable mismatches that
  a numbers-only diff hides.
* :func:`scan_parameter` — sweep a single run input (typically the wind speed,
  which FlamMap project files don't store in an easily parsed form) and report
  which value best reproduces the reference. This is a *diagnostic* to pin down
  an unknown run setting, not a substitute for validation: fitting the wind and
  then reporting the resulting agreement would be circular, so the two uses are
  kept distinct and labelled.

Nothing here imports rasterio; callers pass NumPy arrays already on a common grid.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np


@dataclass
class FieldComparison:
    """Comparison metrics between a pyflam field and a reference field."""

    n: int                      # cells compared (in the mask, both finite)
    bias: float                 # mean(pyflam - reference), same units
    mae: float                  # mean abs error
    rmse: float                 # root mean square error
    median_ratio: float         # median(pyflam / reference)
    pearson_r: float            # linear correlation
    log_rmse: float             # RMSE of log10 ratio (dimensionless)
    fit_slope: float            # OLS pyflam ~ a*reference + b
    fit_intercept: float
    within_10pct: float         # fraction with |pyflam/ref - 1| <= 0.10
    within_25pct: float
    # Burn / no-burn classification over the *whole* mask (incl. zero cells).
    n_mask: int
    both_burn: int
    both_zero: int
    pyflam_only: int            # pyflam burns, reference doesn't
    reference_only: int         # reference burns, pyflam doesn't

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self, units: str = "ft/min") -> str:
        """A human-readable multi-line report."""
        agree = self.both_burn + self.both_zero
        pct_agree = 100.0 * agree / self.n_mask if self.n_mask else float("nan")
        return "\n".join([
            f"cells compared (both burning): {self.n:,}",
            f"  bias (pyflam - FlamMap):  {self.bias:+.3f} {units}",
            f"  MAE:                      {self.mae:.3f} {units}",
            f"  RMSE:                     {self.rmse:.3f} {units}",
            f"  median ratio pyflam/ref:  {self.median_ratio:.3f}",
            f"  Pearson r:                {self.pearson_r:.4f}",
            f"  log10-ratio RMSE:         {self.log_rmse:.4f}  "
            f"(x{10**self.log_rmse:.2f})",
            f"  OLS fit pyflam = {self.fit_slope:.3f}*ref + {self.fit_intercept:.3f}",
            f"  within 10% / 25%:         "
            f"{100*self.within_10pct:.1f}% / {100*self.within_25pct:.1f}%",
            f"burn/no-burn classification over {self.n_mask:,} cells:",
            f"  agree:        {agree:,} ({pct_agree:.2f}%)",
            f"  pyflam-only burning:   {self.pyflam_only:,}",
            f"  FlamMap-only burning:  {self.reference_only:,}",
        ])


def compare_fields(
    pyflam: np.ndarray,
    reference: np.ndarray,
    *,
    mask: np.ndarray | None = None,
    burn_threshold: float = 1e-6,
) -> FieldComparison:
    """Compare a pyflam field to a reference (FlamMap) field, cell by cell.

    Both arrays must be on the same grid. ``mask`` (bool, same shape) selects the
    cells to consider — typically "valid in both rasters" (not nodata). Cells are
    classed as burning where the value exceeds ``burn_threshold``. The numeric
    metrics (bias, RMSE, ratios, correlation, fit) are computed over cells that
    burn in *both*; the classification counts use the whole mask.
    """
    a = np.asarray(pyflam, dtype=float)
    b = np.asarray(reference, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: pyflam {a.shape} vs reference {b.shape}")

    if mask is None:
        mask = np.ones(a.shape, dtype=bool)
    mask = mask & np.isfinite(a) & np.isfinite(b)
    n_mask = int(mask.sum())

    a_burn = mask & (a > burn_threshold)
    b_burn = mask & (b > burn_threshold)
    both = a_burn & b_burn

    both_burn = int(both.sum())
    both_zero = int((mask & ~a_burn & ~b_burn).sum())
    pyflam_only = int((a_burn & ~b_burn).sum())
    reference_only = int((b_burn & ~a_burn).sum())

    av = a[both]
    bv = b[both]
    n = av.size
    if n == 0:
        nan = float("nan")
        return FieldComparison(
            n=0, bias=nan, mae=nan, rmse=nan, median_ratio=nan, pearson_r=nan,
            log_rmse=nan, fit_slope=nan, fit_intercept=nan,
            within_10pct=nan, within_25pct=nan, n_mask=n_mask,
            both_burn=both_burn, both_zero=both_zero,
            pyflam_only=pyflam_only, reference_only=reference_only,
        )

    diff = av - bv
    ratio = av / bv
    log_ratio = np.log10(ratio)
    if n > 1 and av.std() > 0.0 and bv.std() > 0.0:
        pearson = float(np.corrcoef(av, bv)[0, 1])
        slope, intercept = np.polyfit(bv, av, 1)
    else:
        pearson = slope = intercept = float("nan")

    return FieldComparison(
        n=n,
        bias=float(diff.mean()),
        mae=float(np.abs(diff).mean()),
        rmse=float(np.sqrt((diff ** 2).mean())),
        median_ratio=float(np.median(ratio)),
        pearson_r=pearson,
        log_rmse=float(np.sqrt((log_ratio ** 2).mean())),
        fit_slope=float(slope),
        fit_intercept=float(intercept),
        within_10pct=float((np.abs(ratio - 1.0) <= 0.10).mean()),
        within_25pct=float((np.abs(ratio - 1.0) <= 0.25).mean()),
        n_mask=n_mask,
        both_burn=both_burn,
        both_zero=both_zero,
        pyflam_only=pyflam_only,
        reference_only=reference_only,
    )


@dataclass
class DirectionComparison:
    """Angular agreement between two compass-direction fields (degrees)."""

    n: int
    mean_abs_deg: float         # mean |wrapped difference|
    median_abs_deg: float
    rmse_deg: float
    mean_cos: float             # mean cos(difference): 1 = identical, 0 = orthogonal
    within_5deg: float
    within_10deg: float
    within_30deg: float

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        return "\n".join([
            f"cells compared: {self.n:,}",
            f"  mean |angular error|:   {self.mean_abs_deg:.2f} deg",
            f"  median |angular error|: {self.median_abs_deg:.2f} deg",
            f"  RMSE:                   {self.rmse_deg:.2f} deg",
            f"  mean cos(error):        {self.mean_cos:.4f}",
            f"  within 5 / 10 / 30 deg: "
            f"{100*self.within_5deg:.1f}% / {100*self.within_10deg:.1f}% / "
            f"{100*self.within_30deg:.1f}%",
        ])


def compare_directions(
    pyflam_deg: np.ndarray,
    reference_deg: np.ndarray,
    *,
    mask: np.ndarray | None = None,
) -> DirectionComparison:
    """Compare two compass-direction fields (degrees) with circular statistics.

    Differences are wrapped to [-180, 180] before any averaging — you can't take
    a plain mean of angles. Both arrays must share a grid; ``mask`` selects cells.
    """
    a = np.asarray(pyflam_deg, dtype=float)
    b = np.asarray(reference_deg, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    if mask is None:
        mask = np.ones(a.shape, dtype=bool)
    mask = mask & np.isfinite(a) & np.isfinite(b)
    av, bv = a[mask], b[mask]
    n = av.size
    if n == 0:
        nan = float("nan")
        return DirectionComparison(0, nan, nan, nan, nan, nan, nan, nan)

    diff = (av - bv + 180.0) % 360.0 - 180.0
    absd = np.abs(diff)
    return DirectionComparison(
        n=n,
        mean_abs_deg=float(absd.mean()),
        median_abs_deg=float(np.median(absd)),
        rmse_deg=float(np.sqrt((diff ** 2).mean())),
        mean_cos=float(np.cos(np.radians(diff)).mean()),
        within_5deg=float((absd <= 5.0).mean()),
        within_10deg=float((absd <= 10.0).mean()),
        within_30deg=float((absd <= 30.0).mean()),
    )


@dataclass
class PerimeterComparison:
    """Overlap and distance agreement between two burned areas / perimeters."""

    n_pyflam: int
    n_reference: int
    n_intersection: int
    jaccard: float               # |A & B| / |A | B|  (IoU)
    dice: float                  # 2|A & B| / (|A| + |B|)
    area_ratio: float            # |A| / |B|
    mean_perimeter_distance: float   # symmetric, map units
    hausdorff: float                 # symmetric max, map units

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self, units: str = "m") -> str:
        return "\n".join([
            f"burned area: pyflam {self.n_pyflam:,} vs FlamMap "
            f"{self.n_reference:,} cells (ratio {self.area_ratio:.3f})",
            f"  Jaccard (IoU):            {self.jaccard:.3f}",
            f"  Dice:                     {self.dice:.3f}",
            f"  mean perimeter distance:  {self.mean_perimeter_distance:.1f} {units}",
            f"  Hausdorff distance:       {self.hausdorff:.1f} {units}",
        ])


def compare_perimeters(pyflam_burned, reference_burned, *,
                       cellsize_x=1.0, cellsize_y=None) -> PerimeterComparison:
    """Compare two burned-area masks: overlap (Jaccard/Dice) + perimeter distance.

    ``pyflam_burned`` / ``reference_burned`` are boolean arrays on the same grid.
    Distances use the cell size (anisotropic via ``cellsize_y``) and are the
    symmetric mean / max nearest-neighbour distance between the two perimeters
    (Hausdorff). The acceptance test for the MTT growth engine vs a FlamMap
    single-fire (spotting-off) perimeter.
    """
    from scipy.ndimage import binary_erosion, distance_transform_edt
    a = np.asarray(pyflam_burned, dtype=bool)
    b = np.asarray(reference_burned, dtype=bool)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    csy = cellsize_x if cellsize_y is None else cellsize_y
    inter = int((a & b).sum())
    union = int((a | b).sum())
    na, nb = int(a.sum()), int(b.sum())

    pa = a & ~binary_erosion(a) if na else a    # perimeter (edge) cells
    pb = b & ~binary_erosion(b) if nb else b
    sampling = (csy, cellsize_x)
    if pa.any() and pb.any():
        d_to_b = distance_transform_edt(~pb, sampling=sampling)
        d_to_a = distance_transform_edt(~pa, sampling=sampling)
        da, db = d_to_b[pa], d_to_a[pb]
        mean_dist = 0.5 * (float(da.mean()) + float(db.mean()))
        hausdorff = max(float(da.max()), float(db.max()))
    else:
        mean_dist = hausdorff = float("nan")

    return PerimeterComparison(
        n_pyflam=na, n_reference=nb, n_intersection=inter,
        jaccard=inter / union if union else float("nan"),
        dice=2 * inter / (na + nb) if (na + nb) else float("nan"),
        area_ratio=na / nb if nb else float("nan"),
        mean_perimeter_distance=mean_dist, hausdorff=hausdorff)


def compare_arrival_times(pyflam_arrival, reference_arrival, *,
                          max_time=float("inf"), mask=None) -> FieldComparison:
    """Compare two fire arrival-time rasters over the cells both reach in time.

    Wraps :func:`compare_fields` on cells that burn in *both* within ``max_time``
    (so the metrics describe timing agreement where the fires overlap, not the
    extent disagreement — use :func:`compare_perimeters` for extent).
    """
    a = np.asarray(pyflam_arrival, dtype=float)
    b = np.asarray(reference_arrival, dtype=float)
    both = np.isfinite(a) & np.isfinite(b) & (a <= max_time) & (b <= max_time)
    if mask is not None:
        both = both & mask
    return compare_fields(a, b, mask=both, burn_threshold=-np.inf)


def scan_parameter(run, values, reference, *, mask=None, metric="log_rmse"):
    """Sweep one run input and report which value best matches the reference.

    ``run`` is a callable ``value -> pyflam_field`` (e.g. a closure that runs the
    spread model at a given wind speed). Returns ``(best_value, best_comparison,
    table)`` where ``table`` is a list of ``(value, FieldComparison)``. ``metric``
    is the attribute of :class:`FieldComparison` to minimize (``"log_rmse"`` by
    default — scale-free, the right choice for a quantity spanning decades).

    This is a diagnostic for recovering an unknown setting (e.g. the wind speed a
    FlamMap project used); it is not itself a validation result.
    """
    table = []
    best_value = None
    best_cmp = None
    best_score = np.inf
    for v in values:
        cmp = compare_fields(run(v), reference, mask=mask)
        table.append((v, cmp))
        score = abs(getattr(cmp, metric))
        if np.isfinite(score) and score < best_score:
            best_score, best_value, best_cmp = score, v, cmp
    return best_value, best_cmp, table
