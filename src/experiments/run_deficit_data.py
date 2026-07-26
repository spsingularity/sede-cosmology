#!/usr/bin/env python3
"""The deficit from data : what does model-independent H(z) say
about the approach to the de Sitter asymptote?

Data: the standard 31-point cosmic-chronometer compilation (0.07 <= z <= 1.965;
Moresco-type differential-age measurements; the widely used table, fetched and
cross-checked 2026-07-26). Diagonal errors (systematic covariance not included
-- caveat in the finding).

Fits (H0 profiled analytically at each shape point -- chi2 is quadratic in H0):
  1. flat LCDM  (Om)                 -- the rate -3 class
  2. SEDE volume, gamma = 1.5 FIXED  (Om)  -- the growth-gated rate -2 class,
     same parameter count as flat LCDM: a direct chi2 comparison
  3. oLCDM      (Om, Ok)             -- exposes the deficit/curvature degeneracy
Plus: the effective curvature of SEDE (fit oLCDM to SEDE's own best-fit curve),
and the deficit trajectories Delta(z) = ln(H/H_inf) of the best fits, including
the x* family degeneracy of the asymptote (second-law cap family).

Run: python3 src/experiments/run_deficit_data.py
"""
import json, os
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from sede.friedmann import E_SEDE_volume

CC = np.array([  # z, H [km/s/Mpc], sigma
    [0.07, 69.0, 19.6], [0.09, 69.0, 12.0], [0.12, 68.6, 26.2], [0.17, 83.0, 8.0],
    [0.179, 75.0, 4.0], [0.199, 75.0, 5.0], [0.2, 72.9, 29.6], [0.27, 77.0, 14.0],
    [0.28, 88.8, 36.6], [0.352, 83.0, 14.0], [0.3802, 83.0, 13.5], [0.4, 95.0, 17.0],
    [0.4004, 77.0, 10.2], [0.4247, 87.1, 11.2], [0.4497, 92.8, 12.9], [0.47, 89.0, 50.0],
    [0.4783, 80.9, 9.0], [0.48, 97.0, 62.0], [0.593, 104.0, 13.0], [0.68, 92.0, 8.0],
    [0.78, 105.0, 12.0], [0.875, 125.0, 17.0], [0.88, 90.0, 40.0], [0.9, 117.0, 23.0],
    [1.037, 154.0, 20.0], [1.3, 168.0, 17.0], [1.363, 160.0, 33.6], [1.43, 177.0, 18.0],
    [1.53, 140.0, 14.0], [1.75, 202.0, 40.0], [1.965, 186.5, 50.4]])
z, Hd, sig = CC.T
GAM = 1.50


def fit_shape(E_of_z):
    """Profile H0 analytically: chi2(H0) = sum ((Hd - H0 E)/sig)^2."""
    E = E_of_z
    H0 = np.sum(Hd * E / sig ** 2) / np.sum(E ** 2 / sig ** 2)
    chi2 = float(np.sum(((Hd - H0 * E) / sig) ** 2))
    return H0, chi2


def E_lcdm(zz, Om, Ok=0.0):
    return np.sqrt(Om * (1 + zz) ** 3 + Ok * (1 + zz) ** 2 + (1 - Om - Ok))


def scan_lcdm():
    best = (1e9,)
    for Om in np.arange(0.15, 0.55, 0.0025):
        H0, c2 = fit_shape(E_lcdm(z, Om))
        if c2 < best[0]:
            best = (c2, Om, H0)
    return dict(chi2=best[0], Om=best[1], H0=best[2], k=2)


def scan_sede():
    best = (1e9,)
    for Om in np.arange(0.15, 0.55, 0.005):
        E = np.asarray(E_SEDE_volume(z, Om, GAM), float)
        H0, c2 = fit_shape(E)
        if c2 < best[0]:
            best = (c2, Om, H0)
    return dict(chi2=best[0], Om=best[1], H0=best[2], k=2)


def scan_olcdm(zz=None, Hh=None, ss=None):
    zz = z if zz is None else zz
    Hh = Hd if Hh is None else Hh
    ss = sig if ss is None else ss
    best = (1e9,)
    for Om in np.arange(0.10, 0.60, 0.005):
        for Ok in np.arange(-0.5, 0.5, 0.01):
            E = np.sqrt(np.maximum(
                Om * (1 + zz) ** 3 + Ok * (1 + zz) ** 2 + (1 - Om - Ok), 1e-10))
            H0 = np.sum(Hh * E / ss ** 2) / np.sum(E ** 2 / ss ** 2)
            c2 = float(np.sum(((Hh - H0 * E) / ss) ** 2))
            if c2 < best[0]:
                best = (c2, Om, Ok, H0)
    return dict(chi2=best[0], Om=best[1], Ok=best[2], H0=best[3], k=3)


def main():
    L = scan_lcdm()
    S = scan_sede()
    O = scan_olcdm()
    print(f"flat LCDM : Om={L['Om']:.3f} H0={L['H0']:.1f}  chi2={L['chi2']:.2f} (k=2)")
    print(f"SEDE vol  : Om={S['Om']:.3f} H0={S['H0']:.1f}  chi2={S['chi2']:.2f} (k=2, gamma fixed)")
    print(f"oLCDM     : Om={O['Om']:.3f} Ok={O['Ok']:+.2f} H0={O['H0']:.1f}  chi2={O['chi2']:.2f} (k=3)")
    dchi = S['chi2'] - L['chi2']
    print(f"\n Delta chi2 (SEDE - LCDM, same k): {dchi:+.2f}  on 31 points")

    # effective curvature of SEDE: fit oLCDM to SEDE's own best-fit curve with CC-like errors
    zg = np.linspace(0.05, 1.95, 60)
    E_S = np.asarray(E_SEDE_volume(zg, S['Om'], GAM), float) * S['H0']
    sg = 0.05 * E_S                                   # 5% mock errors
    Ofit = scan_olcdm(zg, E_S, sg)
    print(f" SEDE's effective curvature (oLCDM fit to SEDE curve): Ok_eff = {Ofit['Ok']:+.3f} "
          f"(Om={Ofit['Om']:.3f}), residual chi2={Ofit['chi2']:.2f}/60pts at 5% errors")

    # deficit trajectories: Delta(z) = ln(H/H_inf); LCDM vs SEDE x* family
    ODE0 = 1 - S['Om']
    g = lambda x: 1 - np.exp(-GAM * x)
    x_inf = 1.985                                     # self-consistent attractor (the deficit analysis)
    members = {}
    for xs, label in [(1.0, "x*=1 (clip)"), (1.8, "x*=1.8 (GSL edge)"),
                      (x_inf, "x*=x_inf (uncapped)")]:
        f_inf = g(min(xs, x_inf)) / g(1.0)
        E_inf = ODE0 * f_inf                          # volume closure: E_inf = b_inf
        members[label] = dict(x_star=xs, f_inf=f_inf, E_inf=E_inf,
                              deficit_today=float(np.log(1.0 / E_inf)))
    E_inf_L = np.sqrt(1 - L['Om'])
    print("\n deficit today Delta0 = ln(H0/H_inf):")
    print(f"   LCDM: {np.log(1/E_inf_L):.3f}")
    for lab, m in members.items():
        print(f"   SEDE {lab}: {m['deficit_today']:.3f}")
    print("\n NB all SEDE members fit the z>=0 data IDENTICALLY (second-law cap family): the deficit")
    print(" NORMALIZATION is x*-degenerate; only the SHAPE of H(z) is measured.")

    json.dump(dict(lcdm=L, sede=S, olcdm=O, delta_chi2_sede_minus_lcdm=dchi,
                   sede_effective_curvature=Ofit, deficit_members=members,
                   deficit_today_lcdm=float(np.log(1 / E_inf_L)),
                   n_points=len(z)),
              open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
                  os.path.abspath(__file__)))), "results", "deficit_data.json"), "w"),
              indent=2)
    print("\n wrote results/deficit_data.json")


if __name__ == "__main__":
    main()
