"""Does the coarse-resolution stability metric reproduce the full-resolution truth?

For every JJA 12Z sounding in the period of record:
  TRUTH   : mixed-layer depth by the parcel method + dtheta/dz *inside* that layer,
            both at full sounding resolution. "pyroCu-capable" = truth gradient <= 1e-3.
  COARSE  : the same column subsampled to ICON-2I's 5 pressure levels, then run through
            pyflam's shipped Rib ABL and the surface->ABL bulk gradient.

Then: (a) quantify the ABL and gradient bias, (b) ask how well the coarse metric
reproduces the truth classification at the inherited 1.0e-3 threshold, and (c) find the
threshold on the coarse metric that best reproduces it (Youden J). Also tests a
parcel-based coarse ABL as an alternative estimator.
"""
import sys, numpy as np
from sounding_degradation import (parse_igra, prep, abl_rib, parcel_abl, ml_grad_inside,
                     ml_grad_surface_to_abl, degrade, ELEV)

THR = 1.0e-3


def rates(truth, pred):
    tp = int((truth & pred).sum()); tn = int((~truth & ~pred).sum())
    fp = int((~truth & pred).sum()); fn = int((truth & ~pred).sum())
    sens = tp / max(tp + fn, 1); spec = tn / max(tn + fp, 1)
    return sens, spec, (tp + tn) / max(len(truth), 1), sens + spec - 1


def run(stn, name, months=(6, 7, 8)):
    snds = parse_igra(f"por_{stn}/{stn}-data.txt", months=months, hour=12)
    R = []
    for s in snds:
        try:
            p, z, t, rh, u, v = prep(s, ELEV[stn])
        except Exception:
            continue
        if len(p) < 30 or z.max() < 4000:
            continue
        ml_par = parcel_abl(p, z, t)                     # truth mixing depth
        if ml_par < 300:                                 # not a convective column
            continue
        g_truth = ml_grad_inside(p, z, t, ml_par)
        a_full = abl_rib(p, z, t, rh, u, v)
        try:
            pc, zc, tc, rhc, uc, vc = degrade(p, z, t, rh, u, v)
        except Exception:
            continue
        if len(pc) < 3:
            continue
        a_crs = abl_rib(pc, zc, tc, rhc, uc, vc)         # Rib on 5 levels (what we ship)
        a_crs_par = parcel_abl(pc, zc, tc)               # parcel on 5 levels (alternative)
        g_crs = ml_grad_surface_to_abl(pc, zc, tc, a_crs)
        g_crs_par = ml_grad_surface_to_abl(pc, zc, tc, a_crs_par)
        g_full = ml_grad_surface_to_abl(p, z, t, a_full)
        if not np.isfinite([g_truth, g_crs, g_crs_par, g_full]).all():
            continue
        R.append((ml_par, a_full, a_crs, a_crs_par, g_truth, g_full, g_crs, g_crs_par))

    if len(R) < 20:
        print(f"\n{name}: only {len(R)} usable soundings -- too few\n"); return
    r = np.array(R)
    ml_par, a_full, a_crs, a_crs_par = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
    g_truth, g_full, g_crs, g_crs_par = r[:, 4], r[:, 5], r[:, 6], r[:, 7]
    med = np.nanmedian

    print(f"\n{'='*92}\n{name} -- {len(r)} convective JJA 12Z soundings\n{'='*92}")
    print("MIXING DEPTH (m, median)")
    print(f"  parcel, full-res (truth)      : {med(ml_par):7.0f}")
    print(f"  Rib,    full-res              : {med(a_full):7.0f}   ({med(a_full-ml_par):+.0f} vs truth)")
    print(f"  Rib,    ICON 5 levels         : {med(a_crs):7.0f}   ({med(a_crs-ml_par):+.0f} vs truth)")
    print(f"  parcel, ICON 5 levels         : {med(a_crs_par):7.0f}   ({med(a_crs_par-ml_par):+.0f} vs truth)")

    print("\nMIXED-LAYER dtheta/dz (K/m, median)")
    print(f"  inside ML, full-res (TRUTH)   : {med(g_truth):.2e}")
    print(f"  sfc->ABL,  full-res           : {med(g_full):.2e}")
    print(f"  sfc->ABL,  ICON 5 lev (Rib)   : {med(g_crs):.2e}   <- what pyflam computes")
    print(f"  sfc->ABL,  ICON 5 lev (parcel): {med(g_crs_par):.2e}")

    truth = g_truth <= THR
    print(f"\nCLASSIFICATION vs TRUTH  (pyroCu-capable = dtheta/dz <= {THR:.0e})")
    print(f"  truth base rate: {100*truth.mean():.1f}% of convective soundings are capable")
    print(f"  {'metric':<28} {'flagged':>8} {'sens':>6} {'spec':>6} {'acc':>6} {'J':>6}")
    for lbl, gv in [("sfc->ABL 5lev Rib  @1e-3", g_crs),
                    ("sfc->ABL 5lev parcel@1e-3", g_crs_par),
                    ("sfc->ABL full-res  @1e-3", g_full)]:
        pred = gv <= THR
        se, sp, ac, J = rates(truth, pred)
        print(f"  {lbl:<28} {100*pred.mean():7.1f}% {se:6.2f} {sp:6.2f} {ac:6.2f} {J:6.2f}")

    # Recalibrate: which threshold on the shipped coarse metric best reproduces truth?
    best = max(((t, rates(truth, g_crs <= t)) for t in np.arange(0.2e-3, 6.01e-3, 0.05e-3)),
               key=lambda x: x[1][3])
    t, (se, sp, ac, J) = best
    print(f"\n  RECALIBRATED threshold for the shipped coarse metric (max Youden J):")
    print(f"    dtheta/dz <= {t:.2e} K/m  -> sens {se:.2f}  spec {sp:.2f}  acc {ac:.2f}  J {J:.2f}"
          f"  (flags {100*np.mean(g_crs<=t):.1f}%)")

    return r


if __name__ == "__main__":
    run("ITM00016144", "S. Pietro Capofiume (INLAND Po valley -- Tuscany analogue)")
    run("ITM00016245", "Pratica di Mare (COASTAL Rome -- sea-breeze regime)")
