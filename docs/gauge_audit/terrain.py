import sys, runpy
sys.argv = ["x"]
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(sys.path[0] + "/audit.py")
A = g["A"]; install = g["install"]; pft = g["pft"]; elevated = g["elevated"]; orig = g["orig"]
for name, c0 in g["cases"].items():
    for ps in (100000.0, 92500.0, 85000.0, 75000.0):
        c = c0 if ps == 100000.0 else elevated(c0, ps)
        out = []
        for kw in ({}, dict(p0_pa=ps), dict(dtheta_mode="T")):
            install(**kw); out.append(pft(c, beta_steps=1200)[0])
        A.theta_kelvin, A._mixed_layer_state, A._saturation_point = orig["theta"], orig["ml"], orig["sp"]
        f = (1000/(ps/100))**0.286
        print(f"{name:26s} ps={ps/100:5.0f}  (p0/ps)^k={f:.3f}  as-coded={out[0]:7.0f}  p0=ps:{100*(out[1]/out[0]-1):+6.1f}%  moist-per-T:{100*(out[2]/out[0]-1):+6.1f}%")
