#!/usr/bin/env python3
"""E_G consistency check for Paper III §5 (spec: submissions/CHECK_EG_consistency.md).

THE CLAIM UNDER TEST
    "SEDE is GR with smooth dark energy, so it predicts no gravitational slip,
     E_G(z) = Omega_m,0/f(z), within <1% of LCDM ... at 0.3% ... not a discriminator."

WHY IT IS DOUBTED
    No slip constrains phi = psi, i.e. Sigma = mu -- NOT Sigma = 1. Paper IV finds
    mu_inf = 1.051 and +2-3% fsigma8. E_G = Omega_m,0 Sigma / f then moves unless the
    growth enhancement exactly cancels mu.

METHOD (ratios only, so normalisations cancel)
    Integrated in SYNCHRONOUS gauge (mochi-class does not implement newtonian gauge with
    smg), but phi/psi are output as the Bardeen potentials regardless. At k = 0.1 h/Mpc,
    deep sub-horizon, the synchronous d_m matches the comoving Delta_m to well under a
    percent, so the ratios below are unaffected; the k = 0.01 row is shown for context only.
      psi        -> Newtonian potential governing clustering   ->  mu
      (phi+psi)/2-> lensing potential                          ->  Sigma
      slip        -> phi/psi  (== 1 means no slip)
    mu_ratio    = [psi/d_m]_SEDE      / [psi/d_m]_LCDM
    Sigma_ratio = [(phi+psi)/d_m]_SEDE/ [(phi+psi)/d_m]_LCDM
    f           = -(1+z)/2 dlnP/dz
    E_G ratio   = Sigma_ratio / (f_SEDE/f_LCDM)      [Omega_m,0 matched by construction]
"""
import os, sys
import numpy as np
from scipy.interpolate import CubicSpline

MOCHI = os.environ.get("MOCHI_DIR") or os.path.expanduser("~/Projects/usc-program/tools/mochi-class")
import glob
for _b in glob.glob(os.path.join(MOCHI, "build", "lib*")):
    sys.path.insert(0, _b)
os.chdir(MOCHI)                      # stable_params_input/ paths are relative
from classy import Class

H = 0.68
COS = dict(h=H, omega_b=0.02237, omega_cdm=0.30*H**2 - 0.02237, A_s=2.10e-9, n_s=0.965,
           N_ur=3.046, output='mPk,dTk,vTk', gauge='synchronous')
ZS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.0, 1.1]
ZP = ",".join(str(z) for z in ZS)


def run(sede):
    c = Class(); p = dict(COS); p['P_k_max_1/Mpc'] = 2.0; p['z_pk'] = ZP; p['z_max_pk'] = 1.2
    if sede:
        p.update({'Omega_Lambda': 0, 'Omega_fld': 0, 'Omega_smg': -1,
                  'gravity_model': 'stable_params',
                  'smg_file_name': 'stable_params_input/fpab_sede_b020_stable.dat',
                  'parameters_smg': '1.39982e-1', 'method_gr_smg': 'on', 'z_gr_smg': 5,
                  'expansion_model': 'rho_de', 'expansion_smg': 0.69991,
                  'expansion_file_name': 'stable_params_input/fpab_sede_b020_rho.dat',
                  'output_background_smg': 3, 'skip_stability_tests_smg': 'no',
                  'method_qs_smg': 'fully_dynamic', 'pert_initial_conditions_smg': 'zero'})
    c.set(p); c.compute(); return c


def at_k(c, z, kh):
    """phi, psi, d_m interpolated to k (h/Mpc) at redshift z."""
    t = c.get_transfer(z=z)
    k = t['k (h/Mpc)']
    out = {}
    for key in ('phi', 'psi', 'd_m'):
        out[key] = float(np.interp(kh, k, t[key]))
    return out


def growth_f(c, kh, z):
    k = kh * H
    lnP = np.array([np.log(c.pk(k, zz)) for zz in ZS])
    sp = CubicSpline(np.array(ZS), lnP)
    return float(-(1 + z) / 2 * sp(z, 1))


L, S = run(False), run(True)
print("=" * 78)
print("E_G CONSISTENCY CHECK  —  Paper III §5")
print("=" * 78)
print(f"  Omega_m: LCDM={L.Omega_m():.5f}  SEDE={S.Omega_m():.5f}   (matched by construction)")
print()

for kh in (0.1, 0.01):
    print(f"--- k = {kh} h/Mpc " + "-" * 58)
    print(f"  {'z':>5} {'slip_S':>9} {'slip_L':>9} {'mu':>9} {'Sigma':>9} {'Sig/mu':>8} "
          f"{'f_S/f_L':>9} {'E_G ratio':>10}")
    worst = 0.0
    for z in [0.0, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]:
        s, l = at_k(S, z, kh), at_k(L, z, kh)
        slip_S = s['phi'] / s['psi']
        slip_L = l['phi'] / l['psi']
        mu = (s['psi'] / s['d_m']) / (l['psi'] / l['d_m'])
        Sig = ((s['phi'] + s['psi']) / s['d_m']) / ((l['phi'] + l['psi']) / l['d_m'])
        fr = growth_f(S, kh, z) / growth_f(L, kh, z)
        eg = Sig / fr
        worst = max(worst, abs(eg - 1))
        print(f"  {z:5.2f} {slip_S:9.5f} {slip_L:9.5f} {mu:9.5f} {Sig:9.5f} "
              f"{Sig/mu:8.5f} {fr:9.5f} {eg:10.5f}")
    print(f"  --> max |E_G(SEDE)/E_G(LCDM) - 1| over z<=1 : {100*worst:.3f} %")
    print()

print("SPEC DECISION RULE: if the residual exceeds ~1%, the 'shared GR null' sentence")
print("and the 0.3% figure must be revised, and E_G becomes a (weak) discriminator.")
