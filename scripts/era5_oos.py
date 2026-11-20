"""ERA5 at the peak-growth hour of the 8 portal fires that carry a label but no sonde.

Out-of-sample arm of the capability-margin test (docs/graf_vs_pyflam_2026-07-26.md sec. 21.9):
these fires have a published pyroconvection class and perimeter isochrones but no sounding, so
the atmosphere comes from reanalysis alone -- which is also how an operational run would work.

Inputs (built by the campaign-side analysis, not by this script):
  docs/pyroconv_validation/oos_fires.json   the 8 fires, from the portal catalogue
  docs/isochrone_firepower.json  perimeter growth rates (scripts/isochrone_firepower.py)
Writes era5o_<slug>_<date>_<hour>.nc into the corpus era5/ dir, plus
era5_oos_index.json beside it.
Needs a configured cdsapi client.
"""
import json

import valdata, os, warnings
warnings.simplefilter("ignore")
from datetime import datetime
import cdsapi

S = os.path.dirname(os.path.abspath(__file__))
oos = json.load(open(os.path.join(valdata.DATA,"oos_fires.json")))
REPO = os.environ.get("PYFLAM_REPO",
    "/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam")
iso = json.load(open(os.path.join(REPO, "docs/isochrone_firepower.json")))

req = {}
for o in oos:
    r = iso.get(o["slug"])
    if r and r["rates"]:
        t0, t1, best = max(r["rates"], key=lambda x: x[2])
        when = datetime.fromisoformat(t0)                     # start of the peak interval
    elif o.get("start"):
        d, m, y = o["start"][:10].split("/")
        when = datetime(int(y), int(m), int(d), int(o["start"][11:13]))
    else:
        continue
    req[o["slug"]] = (when, o["lat"], o["lon"])

c = cdsapi.Client(quiet=True, progress=False)
LEV = ["1000","975","950","925","900","875","850","825","800","775","750","700",
       "650","600","550","500","450","400","350","300","250","200"]
print(f"{len(req)} ERA5 requests at peak-growth hours", flush=True)
idx = {}
for i, (slug, (when, lat, lon)) in enumerate(sorted(req.items()), 1):
    dst = os.path.join(valdata.ERA5, f"era5o_{slug}_{when:%Y%m%d}_{when.hour:02d}.nc")
    idx[slug] = dict(path=dst, when=when.isoformat())
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        print(f"  [{i}/{len(req)}] {slug} {when} cached", flush=True); continue
    try:
        c.retrieve("reanalysis-era5-pressure-levels", {
            "product_type": "reanalysis", "format": "netcdf",
            "variable": ["temperature","relative_humidity","geopotential",
                         "u_component_of_wind","v_component_of_wind","specific_humidity"],
            "pressure_level": LEV, "year": f"{when:%Y}", "month": f"{when:%m}",
            "day": f"{when:%d}", "time": f"{when:%H}:00",
            "area": [lat + 0.5, lon - 0.5, lat - 0.5, lon + 0.5]}, dst)
        print(f"  [{i}/{len(req)}] {slug} {when} OK", flush=True)
    except Exception as e:
        print(f"  [{i}/{len(req)}] {slug} FAILED: {str(e)[:90]}", flush=True)
json.dump(idx, open(os.path.join(valdata.DATA,"era5_oos_index.json"), "w"), indent=1)
print("done", flush=True)
