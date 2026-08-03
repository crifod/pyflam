"""Fetch and filter the MISR Plume Height Project catalogue -- the escape rung's test set.

The `escape` rung of the cost ladder (atmosphere._escape_cost_gw) is Briggs' penetration
criterion, whose functional form is sourced but whose *level* rests on a prefactor no one has
measured: Briggs' own geometry gives 0.179 (T/theta), Tory & Kepert eq 25 expanded on the same
group gives 0.060, and nothing in either paper adjudicates. This script pulls the data that can.

The MISR Plume Height Project (Nelson et al. 2013) digitised every smoke plume MISR could
resolve for 2008-2011 plus the summers of 2017-2018, each with a stereo-retrieved plume-top
height, the terrain beneath it, and the summed MODIS FRP of the fire pixels feeding it. That is
both sides of the escape test -- how much fire, and did the plume get out -- for ~72k plumes,
against the 22 columns with 2 observed pyroCb that docs/graf_vs_pyflam_2026-07-26.md sec. 21.7
had to work with.

Served by MERLIN (the ASDC front end) from a single unauthenticated query endpoint. The whole
catalogue is one ~88 MB request, so it is fetched once and cached.

    python scripts/misr_plume_catalogue.py [--out DIR] [--refresh]

Writes `misr_plumes.json` (raw) and `misr_plumes_usable.csv` (filtered, with the derived
columns the calibration needs).
"""
import argparse, collections, datetime as dt, json, os, sys, urllib.parse, urllib.request

import numpy as np

MERLIN = "https://l0dup05.larc.nasa.gov/merlin/merlin/query/"
# The query is positional: biomes, 'eb', regions, 'er', then max/min pairs for plume height,
# FRP, AOD and single-scattering albedo, then the date range. Passing every biome and region
# with the sliders at their full range (newview.js `ranges`) returns the entire catalogue.
BIOMES, REGIONS = range(18), range(13)
SLIDERS = "9000,1,50000,1,3,0,1,0"
DATE_RANGE = ("2008-01-01", "2018-12-31")

REGION_NAME = {0: "Africa", 1: "Australia", 2: "Europe/W. Asia", 3: "North America",
               4: "N. Asia", 5: "South America", 6: "S./E. Asia"}


def fetch(out_dir, refresh=False):
    path = os.path.join(out_dir, "misr_plumes.json")
    if os.path.exists(path) and not refresh:
        print(f"[misr] cached {path} ({os.path.getsize(path)/1e6:.1f} MB)")
        return path
    q = (",".join(str(b) for b in BIOMES) + ",eb," + ",".join(str(r) for r in REGIONS)
         + ",er," + SLIDERS + "," + ",".join(DATE_RANGE))
    url = MERLIN + "?" + urllib.parse.urlencode({"q": q})
    print(f"[misr] fetching the full catalogue (~88 MB) ...")
    with urllib.request.urlopen(url, timeout=600) as r, open(path, "wb") as f:
        f.write(r.read())
    print(f"[misr] wrote {path} ({os.path.getsize(path)/1e6:.1f} MB)")
    return path


def usable(records):
    """Filter to plumes that can carry the escape test, and derive its two observables.

    Rejected: cloud-contaminated retrievals, plumes with no MODIS fire pixel (no firepower to
    test), and any plume whose top does not sit above its own terrain (a failed stereo
    retrieval, not a plume that stayed on the ground).

    `plume_top_agl` is the stereo height minus the mean terrain under the digitised plume.
    MISR heights are above the ellipsoid, so subtracting the *co-located* terrain is what makes
    them comparable to a reanalysis boundary-layer depth, which is above ground.
    """
    out = []
    for r in records:
        if r.get("p_clouds"):
            continue
        frp = r.get("p_total_frp")
        if not frp or float(frp) <= 0 or not r.get("p_num_fire_pts"):
            continue
        top, terr = r.get("p_max_ht"), r.get("p_terr_mean_ht")
        med = r.get("p_med_ht")
        if top is None or terr is None or med is None:
            continue
        agl_top, agl_med = float(top) - float(terr), float(med) - float(terr)
        if not (np.isfinite(agl_top) and agl_top > 0):
            continue
        t = dt.datetime.strptime(r["p_date"], "%Y-%m-%dT%H:%M:%SZ")
        out.append(dict(
            name=r["p_name"], date=t.isoformat(), year=t.year, month=t.month, hour=t.hour,
            lat=float(r["p_src_lat"]), lon=float(r["p_src_long"]),
            region=r["p_region_id"], biome=r["p_biome_id"],
            frp_mw=float(frp), frp_max_mw=float(r.get("p_max_frp") or np.nan),
            n_fire_pts=int(r["p_num_fire_pts"]),
            plume_top_agl_m=agl_top, plume_med_agl_m=agl_med,
            terrain_m=float(terr), area_km2=r.get("p_area"),
            n_heights=r.get("p_num_hts"), url=r.get("p_url")))
    return out


def describe(rows):
    print(f"\n[misr] usable plumes: {len(rows)}")
    frp = np.array([r["frp_mw"] for r in rows])
    top = np.array([r["plume_top_agl_m"] for r in rows])
    for name, a, unit in [("total FRP", frp, "MW"), ("plume top AGL", top, "m")]:
        q = np.percentile(a, [1, 10, 25, 50, 75, 90, 99])
        print(f"  {name:14s} ({unit:2s}) p01..p99: " + " ".join(f"{x:8.0f}" for x in q))
    by_region = collections.Counter(r["region"] for r in rows)
    print("  by region: " + ", ".join(
        f"{REGION_NAME.get(k, k)} {v}" for k, v in by_region.most_common()))
    print("  by year:   " + ", ".join(
        f"{k} {v}" for k, v in sorted(collections.Counter(r['year'] for r in rows).items())))
    # The UTC overpass window per region decides how much reanalysis the calibration must pull:
    # MISR is sun-synchronous, so a region occupies only a few hours of the day.
    for reg, n in by_region.most_common(3):
        hrs = sorted({r["hour"] for r in rows if r["region"] == reg})
        print(f"  {REGION_NAME.get(reg, reg):15s} n={n:6d}  UTC hours {hrs}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "..", "data", "misr"))
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)

    rows = usable(json.load(open(fetch(out, a.refresh))))
    describe(rows)

    csv = os.path.join(out, "misr_plumes_usable.csv")
    cols = list(rows[0].keys())
    with open(csv, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    print(f"\n[misr] wrote {csv} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
