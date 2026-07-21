#!/usr/bin/env python3
"""
T2.2 — Information criteria (AIC / BIC / DIC) alongside ΔDIC  [REVISION_PLAN.md]
================================================================================
GP "E": DIC alone is distrusted for non-Gaussian posteriors; report AIC/BIC too.

The key fact: SEDE's dark sector adds NO fitted parameter, so k_SEDE = k_LCDM = 5.
With Δk = 0, AIC = χ²+2k, BIC = χ²+k·lnN and DIC ≈ χ²+p_D all reduce to the SAME
Δ = Δχ²_min in the model comparison — the preference is penalty-free, not bought with
complexity. This script reports the absolute AIC/BIC and confirms ΔAIC = ΔBIC = Δχ².

χ²_min from the canonical compressed-CMB joint (run_probe_decomposition.json);
full-CMB joint Δχ² from run_joint_fullcmb.json.  Run: python run_info_criteria.py
"""
import json
import numpy as np

K = 5                      # Ω_m, H0, ω_b, M_B, σ8 — same for SEDE and ΛCDM (dark sector fixed)
N = 1628                   # BAO 13 + SN 1580 + CC 15 + fσ8 16 + CMB 2 + SH0ES 1 + ω_b 1


def crit(chi2, k=K, n=N):
    return dict(chi2=chi2, AIC=chi2 + 2 * k, BIC=chi2 + k * np.log(n))


def main():
    pd = json.load(open('results/probe_decomposition.json'))
    chiL, chiS = pd['lcdm']['TOTAL'], pd['sede']['TOTAL']
    cL, cS = crit(chiL), crit(chiS)
    print("=== Information criteria (compressed-CMB joint; k=5 both, N=1628) ===\n")
    print(f"{'':6s} {'χ²_min':>10s} {'AIC':>10s} {'BIC':>10s}")
    print(f"{'ΛCDM':6s} {cL['chi2']:10.2f} {cL['AIC']:10.2f} {cL['BIC']:10.2f}")
    print(f"{'SEDE':6s} {cS['chi2']:10.2f} {cS['AIC']:10.2f} {cS['BIC']:10.2f}")
    dchi = cS['chi2'] - cL['chi2']
    print(f"\n  ΔAIC = {cS['AIC']-cL['AIC']:+.2f}   ΔBIC = {cS['BIC']-cL['BIC']:+.2f}   "
          f"(Δχ²_min = {dchi:+.2f}; ΔDIC ≈ {dchi:+.2f}, headline)")
    print("  Δk = 0 ⟹ all information criteria coincide at Δχ²: the preference is penalty-free.")
    try:
        jf = json.load(open('results/joint_fullcmb.json'))
        d = jf['delta_chi2_sede_minus_lcdm_fullcmb']
        print(f"\n  Full-CMB joint: ΔAIC = ΔBIC = Δχ² = {d:+.2f} (same Δk=0 argument).")
    except Exception:
        pass
    json.dump({'lcdm': cL, 'sede': cS, 'dAIC': cS['AIC']-cL['AIC'],
               'dBIC': cS['BIC']-cL['BIC'], 'N': N, 'k': K},
              open('results/info_criteria.json', 'w'), indent=2)
    print("  saved -> results/info_criteria.json")


if __name__ == '__main__':
    main()
