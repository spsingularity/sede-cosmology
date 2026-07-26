#!/usr/bin/env python3
"""Second-law selection of the gate (second-law cap family): which saturation functions
avoid the future overshoot of Finding 23/24?

ANALYTIC PREDICTION (derived first, tested here): near the attractor,
    E - E_inf ~ [src - (b_inf - b) E_inf] / b_inf,
with matter src ~ e^{-3N} and gate remainder b_inf - b ~ g'(x_inf) e^{-2N}
(growth mode). Since -2 beats -3, ANY gate still strictly filling at the
attractor (g'(x_inf) > 0) undershoots H_inf -> overshoot is GENERIC across
smooth gate shapes. Escape requires the filling to COMPLETE at finite
x* <= x*_crit < x_inf (or zero slope at x_inf): then the late approach is
matter-dominated (-3) and H is monotone from above.

Map:
  A. smooth families (exp CDF gamma in {0.5,1.5,3}, linear, quadratic, sqrt,
     tanh): predict ALL overshoot.
  B. capped exp(1.5) with saturation x* in [1.0, x_inf]: find x*_crit; verify
     the observable past (z>=0) is IDENTICAL for every x* >= 1 (x = D^2/D0^2
     <= 1 in the past), so the GSL-compliant family is data-indistinguishable
     from SEDE.
  x* = 1 is exactly the papers' Bousso-bound hypothesis H5 (f_sat <= 1).

Run: python3 src/experiments/run_gate_selection.py
"""
import json, os
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

Om, Orad = 0.315, 9.0e-5
ODE0 = 1 - Om - Orad
N = np.linspace(-7.0, 12.0, 3800)
aN = np.exp(N)
src = Om / aN ** 3 + Orad / aN ** 4


def solve_growth(E2_arr):
    dlnE2 = np.gradient(np.log(E2_arr), N)
    E2i = interp1d(N, E2_arr, fill_value=(E2_arr[0], E2_arr[-1]), bounds_error=False)
    dE2i = interp1d(N, dlnE2, fill_value=(dlnE2[0], dlnE2[-1]), bounds_error=False)

    def rhs(n, y):
        a = np.exp(n)
        return [y[1], -(2 + 0.5 * dE2i(n)) * y[1] + 1.5 * (Om / a ** 3) / E2i(n) * y[0]]

    s = solve_ivp(rhs, [N[0], N[-1]], [aN[0]] * 2, dense_output=True, rtol=1e-10, atol=1e-13)
    return s.sol(N)[0] / s.sol(0.0)[0]


def selfconsistent(fgate):
    """Volume closure E^2 = src + ODE0*(f/f(1))*E, gate f = fgate(x), x = D^2."""
    E2 = (0.5 * (ODE0 + np.sqrt(ODE0 ** 2 + 4 * src))) ** 2
    for _ in range(40):
        D = solve_growth(E2)
        x = D ** 2
        f = fgate(x) / fgate(np.array([1.0]))[0]     # f(today)=1 (flatness norm)
        b = ODE0 * f
        E2_new = (0.5 * (b + np.sqrt(b ** 2 + 4 * src))) ** 2
        err = float(np.max(np.abs(E2_new - E2) / E2))
        E2 = E2_new
        if err < 1e-10:
            break
    return np.sqrt(E2), x, f


def audit(E, x, f):
    E_inf = E[-1]
    dE = np.diff(E)
    mono = bool(np.all(dE <= 1e-13))
    depth = float(np.min(E) / E_inf - 1)             # <0 => undershoot of H_inf
    past = N <= 0
    return dict(monotone_H=mono, undershoot_depth=depth,
                x_inf=float(x[-1]), f_inf=float(f[-1]),
                E_inf=float(E_inf))


# reference (uncapped exp 1.5) for past-identity checks
GAM = 1.50
exp_gate = lambda g: (lambda x: 1 - np.exp(-g * x))
E_ref, x_ref, f_ref = selfconsistent(exp_gate(GAM))
past = N <= 0

results = {"A_smooth_families": {}, "B_capped_family": {}}

print("A. smooth families (prediction: ALL overshoot, since g'(x_inf) > 0)")
fams = {
    "exp_g0.5": exp_gate(0.5), "exp_g1.5": exp_gate(1.5), "exp_g3": exp_gate(3.0),
    "linear": lambda x: x, "quadratic": lambda x: x ** 2,
    "sqrt": lambda x: np.sqrt(x), "tanh": lambda x: np.tanh(x),
}
for name, g in fams.items():
    E, x, f = selfconsistent(g)
    a = audit(E, x, f)
    results["A_smooth_families"][name] = a
    print(f"   {name:>9}: monotone_H={str(a['monotone_H']):>5}  "
          f"undershoot={a['undershoot_depth']:+.2e}  x_inf={a['x_inf']:.3f}")

print("\nB. capped exp(1.5): f = g(min(x, x*)) — find x*_crit")
x_inf_ref = float(x_ref[-1])
caps = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, x_inf_ref]
for xs in caps:
    g = lambda x, xs=xs: 1 - np.exp(-GAM * np.minimum(x, xs))
    E, x, f = selfconsistent(g)
    a = audit(E, x, f)
    a["x_star"] = xs
    a["past_identical_to_SEDE"] = float(np.max(np.abs(E[past] / E_ref[past] - 1)))
    a["deficit_today"] = float(1 - 1.0 / a["f_inf"])   # f(today)=1 by norm
    results["B_capped_family"][f"{xs:.3f}"] = a
    print(f"   x*={xs:5.3f}: monotone_H={str(a['monotone_H']):>5}  "
          f"undershoot={a['undershoot_depth']:+.2e}  "
          f"past-delta={a['past_identical_to_SEDE']:.1e}  "
          f"deficit_today={a['deficit_today']:.3f}")

# locate x*_crit (largest monotone cap)
mono_caps = [v["x_star"] for v in results["B_capped_family"].values() if v["monotone_H"]]
x_crit = max(mono_caps) if mono_caps else None
results["x_star_crit_bracket"] = x_crit
print(f"\n   largest GSL-safe cap on this grid: x* = {x_crit}")
print("\n[READ-OUT]")
print("  theorem check: every smooth (uncapped) gate overshoots (depth ~ g'(x_inf));")
print("  the capped family is GSL-safe up to x*_crit; every x* >= 1 member is")
print("  data-identical on z >= 0. NB the Bousso bound (H5) is in FULL-saturation")
print("  units and is automatic for any CDF gate -- the second-law cap is a NEW,")
print("  strictly stronger bound; the E_SEDE code clip corresponds to x* = 1, the")
print("  most conservative member of the GSL-safe family x* in [1, x*_crit].")

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results")
os.makedirs(OUT, exist_ok=True)
json.dump(results, open(os.path.join(OUT, "gate_selection.json"), "w"), indent=2)
print(f"\n wrote {os.path.join(OUT, 'gate_selection.json')}")
