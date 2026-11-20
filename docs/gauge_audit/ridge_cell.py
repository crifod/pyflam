import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("m", sys.argv[1] + "/tests/pyroconv_icon2i_native.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import pyflam_gui.core.pyroconv as pc
from pyflam.atmosphere import theta_kelvin, lcl_height_m, relative_humidity_from_dewpoint, _ML_STABLE, _GAMMA_DEEP
d = m.load(); si = d["idx"][m.HOURS.index(12)]
T2m, Tlev, PS = d["T2m"][si], d["Tlev"][si], d["PS"][si]
RH = relative_humidity_from_dewpoint(T2m, d["Td2m"][si])
cls = m.classify(T2m, RH, Tlev, PS)
th = theta_kelvin(Tlev, np.array(m.LEVELS)[:, None, None])
_, zs, g, abl, ok = pc.surface_layer_terms(T2m, {p: th[i] for i, p in enumerate(m.LEVELS)}, PS, dict(zip(m.LEVELS, m.PLEV_BARO)))
lcl = lcl_height_m(T2m, RH); z = m.PLEV_BARO
gamma = (th[4] - th[3]) / (z[4] - z[3])
print(f"gates: ML stable {_ML_STABLE:.1e} K/m, gamma deep < {_GAMMA_DEEP:.1e} K/m")
for c in (4, 3, 0):
    a, b = [x[len(x)//2] for x in np.where(cls == c)]
    print(f"class {c} @ lat {d['lat'][a]:.2f} lon {d['lon'][b]:.2f}: ps {PS[a,b]:.0f} hPa (z~{zs[a,b]:.0f} m), "
          f"T2m {T2m[a,b]:.1f} C RH {RH[a,b]:.0f}%  th_sfc {theta_kelvin(T2m[a,b], PS[a,b]):.1f}  "
          f"th850 {th[2][a,b]:.1f} th700 {th[3][a,b]:.1f} th500 {th[4][a,b]:.1f}  ML grad {g[a,b]:.1e}  "
          f"ABL {abl[a,b]:.0f} m  LCL {lcl[a,b]:.0f} m  LCL/ABL {lcl[a,b]/max(abl[a,b],100):.2f}  gamma {gamma[a,b]:.1e}")
