import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("m", sys.argv[1] + "/tests/pyroconv_icon2i_fullitaly.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import pyflam_gui.core.pyroconv as pc
from pyflam.atmosphere import relative_humidity_from_dewpoint
d = m.load(); RH = relative_humidity_from_dewpoint(d["T2m"], d["Td2m"])
for hour in (12, 15):
    si = d["idx"][m.HOURS.index(hour)]
    Tl = {p: d["T"][p][si] for p in (850, 700, 500)}
    old = pc.classify(d["T2m"][si], RH[si], Tl)                       # legacy: every cell at 1000 hPa
    new = pc.classify(d["T2m"][si], RH[si], Tl, ps_hpa=d["PS"][si])
    PS = d["PS"][si]; ok = new >= 0
    print(f"{hour:02d}Z  old classes {np.bincount(old.ravel(), minlength=5).tolist()}  "
          f"new {np.bincount(new[ok].ravel(), minlength=5).tolist()} + nodata {int((~ok).sum())}")
    for lo_, hi_ in ((975, 1100), (925, 975), (850, 925)):
        sel = ok & (PS > lo_) & (PS <= hi_)
        print(f"   ps {lo_}-{hi_}: {sel.sum():5d} cells, changed {(old[sel] != new[sel]).sum():4d}, "
              f"old mean {old[sel].mean():.2f} -> new {new[sel].mean():.2f}")
