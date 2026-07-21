#!/usr/bin/env python3
"""
Generate every pre-registerable number from the CANONICAL parameter-free SEDE-H
(Barrow Δ=1 → λ=0.5, γ=Sheth-Tormen=1.4964). No fitting — these are committed predictions.
Output feeds PREREGISTRATION.md.  Run: python gen_predictions.py
"""
import numpy as np
# NumPy 2.0 renamed np.trapz -> np.trapezoid (np.trapz removed in later 2.x). Bind whichever exists so
# this generator reproduces the registered predictions.json byte-for-byte on any NumPy (the two are
# numerically identical). See PREREGISTRATION.md cryptographic lock.
_trapz = getattr(np, "trapezoid", None) or np.trapz
from sede import friedmann as fr
from sede import perturbations as pert
from sede import barrow_bh as bh

Om, Or, GAM, LAM, DEL, H0 = 0.30, 9.0e-5, 1.4964, 0.5, 1.0, 67.5
def E(z): return fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, LAM)
def Elcdm(z): z=np.atleast_1d(z); return np.sqrt(Om*(1+z)**3+Or*(1+z)**4+(1-Om-Or))

P = {}

# --- EOS: w(z), w0, wa (CPL fit 0<z<1), crossing, min ---
# NOTE: grid is z in [0, 2] BY DESIGN (the observable window). 'w_min'/'z_wmin' below are therefore the
# minimum of w over z<=2, NOT a global extremum: w(z) keeps descending to ~-1.16 near z~20 (-> -7/6 matter-era
# asymptote), but Omega_DE <~ 1e-4 there so it is inert. Do NOT widen this grid to "find the true min" — that
# would change predictions.json and break the registered sha256 lock (see PREREGISTRATION.md A5 footnote).
z = np.linspace(0, 2, 200); Ez = E(z); rho = np.maximum(Ez**2 - Om*(1+z)**3 - Or*(1+z)**4, 1e-12)
w = -1 + (1/3.)*(1+z)*np.gradient(np.log(rho), z)
zc = z[z <= 1]; A = np.column_stack([np.ones_like(zc), zc/(1+zc)]); c, *_ = np.linalg.lstsq(A, w[z <= 1], rcond=None)
P['w0'] = float(c[0]); P['wa'] = float(c[1]); P['w0_today'] = float(w[0])
P['w_min'] = float(w.min()); P['z_wmin'] = float(z[np.argmin(w)])
P['crosses_m1'] = bool((w < -1).any()); P['z_cross'] = float(z[np.argmax(w < -1)]) if (w<-1).any() else None

# --- Barrow Δ, λ ---
P['Delta'] = DEL; P['lambda'] = LAM; P['sigma_Delta_DR3Euclid'] = 0.087

# --- early-DE fraction ---
for zz in (3.0, 1100.0):
    Ezz = E(zz)[0]; P[f'OmegaDE_z{int(zz)}'] = float(max(Ezz**2 - Om*(1+zz)**3 - Or*(1+zz)**4, 0)/Ezz**2)

# --- growth: sigma8, S8, growth index, fsigma8 anchors ---
zg = np.array([0.0, 0.5, 1.0]); D, f = fr.compute_growth_model(zg, Om, lambda z: E(z))
s8 = 0.760; P['sigma8'] = s8; P['S8'] = float(s8*np.sqrt(Om/0.3))
Om_z = lambda zz: Om*(1+zz)**3/E(zz)[0]**2
P['gamma_growth_0'] = float(np.log(f[0])/np.log(Om_z(0.0)))   # zg[0]=0.0
P['gamma_growth_1'] = float(np.log(f[2])/np.log(Om_z(1.0)))   # zg[2]=1.0 (f[1] is z=0.5)
zf = np.array([0.3, 0.5, 0.7, 1.0, 1.5]); Df, ff = fr.compute_growth_model(zf, Om, lambda z: E(z))
P['fsigma8'] = {float(zz): float(ff[i]*s8*Df[i]) for i, zz in enumerate(zf)}

# --- ISW relative amplitude ---
P['ISW_over_LCDM'] = float(pert.isw_ratio_cousin(Om, GAM, E_run=lambda zz: E(zz))[0])

# --- age of universe ---
zz = np.concatenate([[0.0], np.logspace(-6, 6, 4000)])
P['t0_Gyr'] = float((9.778/(H0/100.)) * _trapz(1.0/((1+zz)*E(zz)), zz))

# --- BBN null ---
P['BBN_speedup'] = float(E(4e8)[0]/Elcdm(4e8)[0]); P['BBN_dYp'] = 0.16*(P['BBN_speedup']-1)

# --- cross-horizon: BH Barrow signature ---
P['BH_logS_enhancement_60Msun'] = float(bh.entropy_enhancement_log10(60.0, 1.0, 0.7))
P['BH_S_scaling'] = "S ∝ A^{3/2} (Δ=1), not A/4"

# --- statefinder / Om diagnostic fingerprint ---
zz2 = np.linspace(0, 1.5, 6000); Ez2 = E(zz2); Ep = np.gradient(Ez2, zz2); Epp = np.gradient(Ep, zz2)
r = 1 - 2*(1+zz2)*Ep/Ez2 + (1+zz2)**2*(Ep**2 + Ez2*Epp)/Ez2**2
P['statefinder_max_r_dev'] = float(np.max(np.abs(r-1)[40:5960]))

if __name__ == "__main__":
    import json, os
    for k, v in P.items():
        print(f"{k:28s} = {v}")
    _root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    _out = os.path.join(_root, "results", "predictions.json")
    with open(_out, "w") as fh:
        json.dump(P, fh, indent=2, default=str)
    print("\nwrote results/predictions.json")
