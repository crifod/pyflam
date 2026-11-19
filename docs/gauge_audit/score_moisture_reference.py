"""Score PFT moisture_reference="theta" (paper) vs "surface_t" on the GRAF capability-margin set.

Runs scripts/margin_test.py unchanged, once per variant, with the PFT default switched.
Run from the repository root:  python docs/gauge_audit/score_moisture_reference.py <variant>
"""
import functools, runpy, sys
sys.path.insert(0, "src")
import pyflam.atmosphere as A

variant = sys.argv[1] if len(sys.argv) > 1 else "theta"
A.pyrocb_firepower_threshold_grid = functools.partial(A.pyrocb_firepower_threshold_grid,
                                                      moisture_reference=variant)
print(f"### moisture_reference = {variant}")
runpy.run_path("scripts/margin_test.py", run_name="__main__")
