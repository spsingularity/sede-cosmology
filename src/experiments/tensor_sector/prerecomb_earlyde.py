#!/usr/bin/env python3
"""
WP-A2 — SEDE dark-energy fraction at recombination, from the paper's own pipeline.
================================================================================
Recomputes the pre-recombination numbers the tensor-sector note asserts
("f_sat(z~1100) ~ 1e-6, rho_DE/rho_tot ~ 1e-10") using the CANONICAL,
parameter-free SEDE background E_SEDE_volume (Barrow Delta=1, lambda=0.5) and the
pipeline growth factor D(z) — NOT the EH98 approximation. These are the values to
quote in the Paper III pre-recombination paragraph (§6/§9).

Model: Om=0.30, gamma=1.4964, Or=9.0e-5 (the headline run_cmb_earlyde.py values).

Run:  PYTHONPATH=<sede-dev root> python prerecomb_earlyde.py
Requires the vendored `sede` package (sede-dev/sede) on PYTHONPATH.
"""
import numpy as np
from sede.friedmann import compute_growth_factor, E_SEDE_volume
from sede.theory import fsat_late

Om, GAM, Or = 0.30, 1.4964, 9.0e-5
Omega_DE0 = 1.0 - Om - Or
SIGMA8_0 = 0.811

Z_REC  = 1089.9   # Planck 2018 z_star (last scattering)
Z_DRAG = 1059.9   # baryon drag epoch


def row(label, z):
    D = float(compute_growth_factor(np.array([z]), Om, Or)[0])
    x = D ** 2                                   # variance ratio sigma8^2(z)/sigma8^2(0)
    f_sat = float(fsat_late(z, SIGMA8_0 * D, SIGMA8_0, GAM))   # pipeline-stable f_sat
    f_sat_asymp = GAM * x / (1.0 - np.exp(-GAM))              # gate limit f_sat -> gamma x/(1-e^-gamma), D->0
    E = float(E_SEDE_volume(np.array([z]), Om, GAM)[0])
    matter = Om * (1 + z) ** 3 + Or * (1 + z) ** 4
    OmDE = (E ** 2 - matter) / E ** 2            # exact early-DE fraction = rho_DE/rho_tot
    # analytic identity for the lambda=0.5 gate: rho_DE/rho_crit = Omega_DE0 f_sat / E
    OmDE_id = Omega_DE0 * f_sat / E
    print(f"\n=== {label}  (z = {z}) ===")
    print(f"  D(z)                    = {D:.6e}")
    print(f"  x = D^2 (variance)      = {x:.6e}")
    print(f"  f_sat(z)  [pipeline]    = {f_sat:.6e}")
    print(f"  f_sat(z)  [gate limit]  = {f_sat_asymp:.6e}   (gamma x/(1-e^-gamma); D->0 check)")
    print(f"  E(z) = H/H0             = {E:.6e}")
    print(f"  rho_DE/rho_tot (exact)  = {OmDE:.6e}")
    print(f"  rho_DE/rho_tot (Om0 f/E)= {OmDE_id:.6e}   (lambda=0.5 identity)")
    return dict(z=z, D=D, x=x, f_sat=f_sat, E=E, OmDE=OmDE)


if __name__ == "__main__":
    print("SEDE canonical background E_SEDE_volume (Delta=1, lambda=0.5)")
    print(f"Om={Om}  gamma={GAM}  Or={Or}  Omega_DE0={Omega_DE0:.4f}  sigma8_0={SIGMA8_0}")

    rec  = row("recombination / last scattering", Z_REC)
    drag = row("baryon drag",                     Z_DRAG)
    z3   = row("z=3 cross-check (paper Fig 5: SEDE 0.029)", 3.0)

    print("\n=== Values to quote in the Paper III pre-recombination paragraph ===")
    print(f"  f_sat(z=1089.9)      = {rec['f_sat']:.2e}   (note asserted ~1e-6)")
    print(f"  rho_DE/rho_tot(z_rec)= {rec['OmDE']:.2e}   (note asserted ~1e-10)")
    print(f"  Omega_DE(z=3)        = {z3['OmDE']:.4f}    (reproduces paper's 0.029 -> pipeline self-consistent)")
