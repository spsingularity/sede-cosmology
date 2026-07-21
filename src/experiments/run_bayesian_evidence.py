#!/usr/bin/env python3
"""
Bayesian evidence (nested sampling) for SEDE vs ΛCDM
==========================================================================
R6 (statistics) asked for a genuine Bayesian evidence rather than DIC/AIC/BIC, which
coincide here only because Δk = 0. Nested sampling (dynesty) integrates the joint
likelihood over the prior for each model, returning ln Z; the Bayes factor ln B =
ln Z_SEDE − ln Z_ΛCDM folds in the Occam factor automatically.

Both models are sampled over the SAME 5 priors (Ω_m, H0, ω_b, M_B, σ8) on the canonical
CAMB-in-the-loop compressed-CMB joint (run_probe_decomposition.chi2_byprobe). SEDE's dark
sector is fixed (λ=0.5, γ=1.4964), so the two models have identical prior volume — the
Bayes factor is then a clean test of whether the fit improvement survives the Occam penalty.

Run: python run_bayesian_evidence.py [--nlive 400]
"""
from __future__ import annotations
import argparse, json
import numpy as np
import dynesty

import run_probe_decomposition as PD

# uniform priors (identical for both models — dark sector fixed)
PRIORS = [("Om", 0.20, 0.45), ("H0", 60.0, 78.0), ("ombh2", 0.019, 0.025),
          ("MB", -20.5, -18.5), ("s8", 0.62, 0.90)]
LO = np.array([p[1] for p in PRIORS]); HI = np.array([p[2] for p in PRIORS])


def ptform(u):
    return LO + u * (HI - LO)


def make_loglike(model):
    def loglike(theta):
        try:
            return -0.5 * PD.chi2_byprobe(theta, model)['TOTAL']
        except Exception:
            return -1e10
    return loglike


def run_model(model, nlive):
    sampler = dynesty.NestedSampler(make_loglike(model), ptform, ndim=len(PRIORS),
                                    nlive=nlive, bound='multi', sample='rwalk')
    sampler.run_nested(print_progress=False, dlogz=0.1)
    res = sampler.results
    return float(res.logz[-1]), float(res.logzerr[-1])


def main(nlive):
    PD.DATA = PD.load()
    print(f"=== Bayesian evidence (dynesty, nlive={nlive}) — SEDE vs ΛCDM, compressed-CMB joint ===\n")
    lnZ = {}
    for model in ('lcdm', 'sede'):
        z, ze = run_model(model, nlive)
        lnZ[model] = (z, ze)
        print(f"  {model.upper():5s}  ln Z = {z:8.2f} ± {ze:.2f}")
    lnB = lnZ['sede'][0] - lnZ['lcdm'][0]
    err = np.hypot(lnZ['sede'][1], lnZ['lcdm'][1])
    print(f"\n  >>> ln B (SEDE − ΛCDM) = {lnB:+.2f} ± {err:.2f} <<<")
    # Jeffreys/Kass-Raftery scale
    a = abs(lnB)
    scale = ("inconclusive (|lnB|<1)" if a < 1 else
             "positive/weak (1–2.5)" if a < 2.5 else
             "strong (2.5–5)" if a < 5 else "very strong (>5)")
    fav = "SEDE" if lnB > 0 else "ΛCDM"
    print(f"  Jeffreys scale: {scale}, favouring {fav}.")
    print(f"  (Δk=0 ⟹ equal prior volume; ln B ≈ ½Δχ² minus the Occam/prior-shape factor."
          f"  Compare ½Δχ² = {0.5*-2.95:.2f}.)")
    json.dump({'lnZ_lcdm': lnZ['lcdm'], 'lnZ_sede': lnZ['sede'], 'lnB': lnB, 'lnB_err': err,
               'nlive': nlive}, open('results/bayesian_evidence.json', 'w'), indent=2)
    print("  saved -> results/bayesian_evidence.json")


if __name__ == '__main__':
    import os; os.makedirs('results', exist_ok=True)
    ap = argparse.ArgumentParser(); ap.add_argument('--nlive', type=int, default=400)
    main(ap.parse_args().nlive)
