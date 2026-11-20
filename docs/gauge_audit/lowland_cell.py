import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("m", sys.argv[1] + "/tests/pyroconv_icon2i_fullitaly.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import pyflam.atmosphere as A
import pyflam_gui.core.pyroconv as pc
from pyflam.atmosphere import relative_humidity_from_dewpoint, theta_kelvin, lcl_height_m
d = m.load(); RH = relative_humidity_from_dewpoint(d["T2m"], d["Td2m"])
si = d["idx"][m.HOURS.index(12)]
Tl = {p: d["T"][p][si] for p in (850, 700, 500)}; T2m, PS = d["T2m"][si], d["PS"][si]
seen = {}
real = A.pyroconvection_type
def spy(**k):
    seen.setdefault(tag, []).append(k); return real(**k)
A.pyroconvection_type = spy
tag = "old"; old = pc.classify(T2m, RH[si], Tl)
tag = "new"; new = pc.classify(T2m, RH[si], Tl, ps_hpa=PS)
A.pyroconvection_type = real
print("legacy ABL_MIN skip count:", 13920 - len(seen["old"]), " new calls:", len(seen["new"]))
# a lowland cell that changed
idx = np.argwhere((PS > 990) & (old != new) & (new >= 0))[len(np.argwhere((PS > 990) & (old != new) & (new >= 0)))//2]
a, b = idx
print(f"cell ps {PS[a,b]:.0f} hPa T2m {T2m[a,b]:.1f} RH {RH[si][a,b]:.0f}%: old class {old[a,b]} new class {new[a,b]}")
th = {p: theta_kelvin(Tl[p], float(p)) for p in Tl}
print(f"  th_sfc(1000) {theta_kelvin(T2m[a,b],1000.0):.2f}  th_sfc(ps) {theta_kelvin(T2m[a,b],PS[a,b]):.2f}  th850 {th[850][a,b]:.2f} th700 {th[700][a,b]:.2f}  LCL {lcl_height_m(T2m, RH[si])[a,b]:.0f}")
_, zs, g, abl, ok = pc.surface_layer_terms(T2m, th, PS, pc.STD_LEVEL_HEIGHT_M, abl_levels=(850, 700))
print(f"  new: z_sfc {zs[a,b]:.0f} ML grad {g[a,b]:.2e} ABL {abl[a,b]:.0f}")
# find the spy kwargs for this cell: calls are in row-major order over cells that were classified
