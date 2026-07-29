# SPDX-License-Identifier: AGPL-3.0-or-later
"""Per-moment firepower from the Wildfire Data Portal's isochrone perimeters.

Tory & Kepert's appendix D gives firepower as ``FP = alpha * h * w_a * dA/dt``. Everything
except ``dA/dt`` is measured (fuel load from the 10 m FBFM40 map, alpha and h from the paper),
and ``dA/dt`` had been taken from the portal's *fire-level* peak burn ratio -- one number for
the whole event, applied to sondes launched at arbitrary moments within it.

The portal also publishes isochrone perimeters per fire, which give ``dA/dt`` *at the moment a
sonde was launched* rather than at the fire's peak hour. This extracts them, raising per-moment
coverage of the campaign sondes from 11/26 to 22/27 and growth rates from 18 fires to 24.

What that bought is documented in docs/graf_vs_pyflam_2026-07-26.md sec. 21.4-21.5, and it is
not what was expected: per-moment firepower is the physically correct quantity but scores
*worse* than the fire-level peak, because it demotes Guissona -- the campaign's one observed
PyroCb, whose sonde was launched ~4 h before the blow-up, into an interval growing at ~460 ha/h
against the fire's 6000 ha/h peak. The gate is right about the moment and is marked wrong
because the label describes the fire's whole life. A moment-level gate cannot be validated
against fire-level labels; that needs per-launch observed classes, which the published table
does not carry.

Four things the files require care with, all found by inspection rather than documentation:

* **Placemark timestamps come in at least six formats** -- ``15:00``, ``20210724_1547``,
  ``20220717 1630``, ``2025/07/07 13:07``, ``21/6/2025 14:00:00``, ``01/07/2025 17:13`` -- and
  are irregular in spacing, not hourly despite the file names.
* **Some files carry no timestamp in ``<name>`` at all**, holding it instead in the ArcGIS
  attribute table rendered into ``<description>`` HTML, under a field named ``DiHo`` or ``FeHo``
  (dia/fecha + hora). Those files parse to zero isochrones unless the table is read.
* **Some files hold cumulative perimeters, others per-interval bands.** Martorell's polygons
  sum to 203 ha against a reported 197 ha total, so they are bands; a file whose areas increase
  monotonically is cumulative. Decided per fire from the geometry, since nothing declares it.
* **Placemarks carry an ``id`` attribute**, so a ``<Placemark>`` pattern silently matches
  nothing.

**Timestamps are local, not UTC** (``TIME_NOTE``). The files do not say so, but the campaign
sondes settle it: Guissona's sonde was launched at 15:59 UTC into an already-convecting plume,
while its isochrones run 17:13-19:29 -- under a UTC reading the sonde precedes every mapped
perimeter of a fire it was sampling. Reading them as local (CEST) puts it inside the span, and
does the same for Patagual, Junquillos and Vega Honda. Converted here per fire via ``zoneinfo``,
which handles DST from the date; zones assigned from the portal's coordinates, which fall in
four well-separated clusters (Catalonia, Netherlands, Greece, central Chile).

Usage:  PYTHONPATH=src python scripts/isochrone_firepower.py [--kmz-dir DIR] [--out PATH]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import zipfile
from datetime import datetime

import numpy as np

TIME_NOTE = ("isochrone timestamps are published in local time and converted to UTC here, per "
             "the launch-time evidence in the module docstring")

_TS_PATTERNS = [
    ("%Y%m%d_%H%M", r"^\d{8}_\d{4}$"),
    ("%Y%m%d %H%M", r"^\d{8} \d{4}$"),
    ("%Y/%m/%d %H:%M", r"^\d{4}/\d{1,2}/\d{1,2} \d{1,2}:\d{2}$"),
    ("%d/%m/%Y %H:%M:%S", r"^\d{1,2}/\d{1,2}/\d{4} \d{1,2}:\d{2}:\d{2}$"),
    ("%d/%m/%Y %H:%M", r"^\d{1,2}/\d{1,2}/\d{4} \d{1,2}:\d{2}$"),
]

# Attribute-table fields that carry the perimeter time when <name> does not.
_TIME_FIELDS = ("DiHo", "FeHo", "Fecha", "Data", "Hora")

# Portal coordinates fall in four well-separated clusters; box -> IANA zone, DST from the date.
_ZONE_BOXES = [
    ((-56.0, -17.0, -76.0, -66.0), "America/Santiago"),
    ((35.0, 44.0, -10.0, 4.0), "Europe/Madrid"),
    ((49.0, 54.0, 3.0, 8.0), "Europe/Amsterdam"),
    ((34.0, 42.0, 19.0, 30.0), "Europe/Athens"),
]


def zone_for(lat, lon):
    """Portal coordinates -> IANA zone name, or None if outside every known cluster."""
    if lat is None or lon is None:
        return None
    for (s, n, w, e), z in _ZONE_BOXES:
        if s <= lat <= n and w <= lon <= e:
            return z
    return None


def to_utc(when: datetime, zone: str | None) -> datetime:
    """Local isochrone timestamp -> naive UTC. Unknown zone: returned unchanged."""
    if zone is None:
        return when
    from zoneinfo import ZoneInfo
    return (when.replace(tzinfo=ZoneInfo(zone))
                .astimezone(ZoneInfo("UTC")).replace(tzinfo=None))


def parse_timestamp(name: str, day: datetime | None = None):
    """Placemark name -> datetime, across the formats the portal actually uses."""
    s = name.strip()
    for fmt, pat in _TS_PATTERNS:
        if re.match(pat, s):
            try:
                return datetime.strptime(s, fmt)
            except ValueError:
                continue
    m = re.match(r"^(\d{1,2})[:.](\d{2})$", s)        # bare time; needs the fire's day
    if m and day is not None:
        return day.replace(hour=int(m.group(1)) % 24, minute=int(m.group(2)),
                           second=0, microsecond=0)
    return None


def timestamp_from_description(desc: str, day: datetime | None = None):
    """Perimeter time from the ArcGIS attribute table rendered into <description> HTML.

    Some exports leave <name> empty or set it to a fire id, and carry the time only here.
    """
    kv = re.findall(r"<td>([^<]{1,20})</td>\s*<td>([^<]{0,40})</td>", desc)
    for key, val in kv:
        if key.strip() in _TIME_FIELDS:
            when = parse_timestamp(val, day)
            if when is not None:
                return when
    return None


def read_isochrones(kmz_path: str, day: datetime | None = None, zone: str | None = None):
    """-> [(datetime_utc, area_m2)] sorted in time, from a portal perimeter KMZ."""
    from pyproj import Transformer
    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    tr = Transformer.from_crs("EPSG:4326", "EPSG:3035", always_xy=True)
    with zipfile.ZipFile(kmz_path) as z:
        kml = z.read([n for n in z.namelist() if n.endswith(".kml")][0]).decode("utf-8", "ignore")

    out = []
    for m in re.finditer(r"<Placemark[^>]*>(.*?)</Placemark>", kml, re.S):
        desc = re.search(r"<description>(.*?)</description>", m.group(1), re.S)
        blk = re.sub(r"<description>.*?</description>", "", m.group(1), flags=re.S)
        nm = re.search(r"<name>([^<]*)</name>", blk)
        when = parse_timestamp(nm.group(1), day) if nm else None
        if when is None and desc is not None:          # fall back to the attribute table
            when = timestamp_from_description(desc.group(1), day)
        if when is None:
            continue
        when = to_utc(when, zone)
        polys = []
        for c in re.findall(r"<coordinates>(.*?)</coordinates>", blk, re.S):
            pts = []
            for tok in c.split():
                p = tok.split(",")
                if len(p) >= 2:
                    try:
                        pts.append(tr.transform(float(p[0]), float(p[1])))
                    except Exception:
                        pass
            if len(pts) >= 4:
                g = Polygon(pts).buffer(0)
                if g.is_valid and g.area > 0:
                    polys.append(g)
        if polys:
            out.append((when, unary_union(polys).area))
    return sorted(out, key=lambda r: r[0])


def growth_rate(iso, total_ha: float | None = None):
    """-> [(t_start, t_end, dA/dt m2/s)] from isochrones, cumulative or banded.

    Whether a file holds cumulative perimeters or per-interval bands is not declared, so it is
    decided from the geometry: monotonically increasing areas are cumulative, and a set whose
    areas *sum* to about the reported total is banded.
    """
    if len(iso) < 2:
        return []
    a = np.array([x[1] for x in iso])
    cumulative = bool(np.all(np.diff(a) >= -1e-6))
    if not cumulative and total_ha is not None:
        # corroborate the banded reading against the reported total where we have one
        summed = a.sum() / 1e4
        if not (0.5 * total_ha <= summed <= 2.0 * total_ha):
            cumulative = True                       # neither reading is clean; prefer differences
    out = []
    for i in range(1, len(iso)):
        dt = (iso[i][0] - iso[i - 1][0]).total_seconds()
        if dt <= 0:
            continue
        da = (a[i] - a[i - 1]) if cumulative else a[i]
        out.append((iso[i - 1][0], iso[i][0], max(da, 0.0) / dt))
    return out


# Growth rates are means over a mapping interval, so coarsely mapped fires under-report their
# peak. Measured, not assumed: taking the 8 fires mapped at <=30 min and progressively
# coarsening their own cumulative-area curves, the peak rate falls as a clean power law
# peak(W) ~ W^-beta over W = native..180 min, with r2 0.76-0.99 per fire and a median
# beta = 0.44 (range 0.08-0.61). A 60 min window costs a factor 2.21 against a 10 min
# reference; the median 60->native loss across those fires is 1.54x.
#
# Linear interpolation inside each interval hides sub-interval bursts, so beta is a LOWER bound.
RESOLUTION_BETA = 0.44
# Reference window: the plume overturning timescale, which is what the convective column
# responds to. Not 1 min -- firepower averaged over a minute is not what lifts a plume.
RESOLUTION_REF_MIN = 10.0
# Beyond ~2 h of mapping gap the power law is pure extrapolation, so the factor is capped.
RESOLUTION_MAX_FACTOR = 3.0


def resolution_factor(median_gap_min: float) -> float:
    """Correction from a fire's mapping cadence to the ``RESOLUTION_REF_MIN`` reference.

    Applies to the *estimate*, and is worth having for that reason alone -- a fire mapped hourly
    genuinely did burn faster than its hourly means say. It does **not** improve
    classification: applied across all 30 labelled columns it shifts every margin up and crosses
    no decision boundary (docs/graf_vs_pyflam_2026-07-26.md sec. 21.11). The gap it closes is
    near-uniform across observed classes, so it is calibration, not discrimination.
    """
    if not np.isfinite(median_gap_min) or median_gap_min <= 0:
        return 1.0
    g = max(float(median_gap_min), RESOLUTION_REF_MIN)
    return float(min((g / RESOLUTION_REF_MIN) ** RESOLUTION_BETA, RESOLUTION_MAX_FACTOR))


def median_gap_minutes(rates) -> float:
    """Median mapping interval (minutes) of a fire's isochrone series."""
    if not rates:
        return float("nan")
    return float(np.median([(t1 - t0).total_seconds() / 60.0 for t0, t1, _ in rates]))


def firepower_at(rates, when, *, w_a=2.0, alpha=0.7, heat=15.0e6, tolerance_s=1800.0):
    """Firepower (GW) at ``when``, from the isochrone interval containing it.

    ``tolerance_s`` carries the last interval's rate forward past the final mapped perimeter.
    Three campaign sondes were launched 14-23 min after their fire's last isochrone while it was
    still burning; the alternative to a bounded carry-forward is discarding them. Held to 30 min
    because mapping often stops near containment, when the rate is falling.
    """
    for t0, t1, dadt in rates:
        if t0 <= when <= t1:
            return alpha * heat * w_a * dadt / 1e9
    if rates and 0 < (when - rates[-1][1]).total_seconds() <= tolerance_s:
        return alpha * heat * w_a * rates[-1][2] / 1e9
    return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmz-dir", required=True)
    ap.add_argument("--out", default="docs/isochrone_firepower.json")
    a = ap.parse_args()
    port = {p["slug"]: p for p in json.load(open("/tmp/portal_fires.json"))}
    res = {}
    for fn in sorted(os.listdir(a.kmz_dir)):
        if not fn.endswith(".kmz"):
            continue
        slug = fn[:-4]
        p = port.get(slug, {})
        day = None
        if p.get("start"):
            try:
                day = datetime.strptime(p["start"][:10], "%d/%m/%Y")
            except ValueError:
                day = None
        zone = zone_for(p.get("lat"), p.get("lon"))
        try:
            iso = read_isochrones(os.path.join(a.kmz_dir, fn), day, zone)
        except Exception as e:
            print(f"  {slug}: {str(e)[:60]}")
            continue
        rates = growth_rate(iso, p.get("area"))
        if not rates:
            continue
        peak = max(r[2] for r in rates) * 3600 / 1e4
        gap = median_gap_minutes(rates)
        rf = resolution_factor(gap)
        res[slug] = dict(n=len(iso), zone=zone,
                         median_gap_min=round(gap, 1), resolution_factor=round(rf, 2),
                         peak_ha_per_h_res_corrected=round(peak * rf, 1),
                         first=iso[0][0].isoformat(), last=iso[-1][0].isoformat(),
                         peak_ha_per_h=round(peak, 1),
                         portal_br_max=p.get("br_max"), portal_area_ha=p.get("area"),
                         rates=[(t0.isoformat(), t1.isoformat(), r) for t0, t1, r in rates])
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"wrote {a.out}: {len(res)} fires with isochrone growth rates")
    print(f"{'fire':26s} {'n':>3s} {'peak ha/h':>10s} {'portal BRmax':>13s}")
    for s, r in sorted(res.items()):
        print(f"{s[:26]:26s} {r['n']:3d} {r['peak_ha_per_h']:10.1f} "
              f"{(r['portal_br_max'] if r['portal_br_max'] else float('nan')):13.0f}")


if __name__ == "__main__":
    main()
