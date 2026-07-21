#!/usr/bin/env python3
"""
T2.3 — Per-probe Δχ² decomposition at the joint best-fit  [REVISION_PLAN.md]
============================================================================
CL (round 2): the leave-one-probe-out range is NOT the decomposition — LOO is blind to a
*paired* tension (two probes in mutual tension that SEDE relieves jointly look stable under
single-probe removal). The right diagnostic is the DIRECT per-probe Δχ² split at the joint
best-fit: which probes' χ² actually went down, and by how much.

This fits SEDE (canonical Barrow λ=0.5, γ=1.4964) and ΛCDM on the SAME canonical CAMB-in-the-loop
joint as the headline (run_lambda_verify; compressed R,ℓ_A CMB), then decomposes the total
Δχ²(SEDE−ΛCDM) into per-probe contributions {DESI BAO, Pantheon+ SN, CMB(R,ℓ_A), ω_b prior,
SH0ES, Moresco CC, fσ8}. Reports whether the preference is spread (consensus) or concentrated
in one/two probes (which would qualify the "consensus" claim).

Run: python run_probe_decomposition.py
"""
from __future__ import annotations
import json
import numpy as np
import camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve

from sede import friedmann as fr
from run_lambda_verify import w_of_a
from run_camb_joint import load, R_PL, LA_PL, CMB_CINV, SHOES, OMBH2_PRIOR

C = 299792.458
GAMMA, LAM = 1.4964, 0.5
DATA = None


def chi2_byprobe(p, model):
    """Per-probe χ² dict at params p=(Om,H0,ombh2,MB,s8) for 'sede' or 'lcdm'."""
    Om, H0, ombh2, MB, s8 = p
    h = H0 / 100.
    pa = camb.CAMBparams()
    pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om * h ** 2 - ombh2, mnu=0.06)
    if model == 'sede':
        a, w = w_of_a(Om, GAMMA, LAM)
        de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False
    bg = camb.get_background(pa)
    d = bg.get_derived_params(); rd, zs, rs = d['rdrag'], d['zstar'], d['rstar']
    o = {}
    # DESI BAO
    z, t, m, icov = DATA['desi']
    pred = np.array([bg.comoving_radial_distance(zz) / rd if tp == 'DM/rd'
                     else (C / bg.hubble_parameter(zz)) / rd if tp == 'DH/rd'
                     else ((zz * bg.comoving_radial_distance(zz) ** 2 * (C / bg.hubble_parameter(zz))) ** (1 / 3.)) / rd
                     for zz, tp in zip(z, t)])
    dd = m - pred; o['DESI_BAO'] = float(dd @ icov @ dd)
    # Pantheon+ SN
    zp, mu, chol = DATA['pan']
    dmu = mu - (5 * np.log10((1 + zp) * bg.comoving_radial_distance(zp)) + 25 + MB)
    o['Pantheon_SN'] = float(dmu @ cho_solve(chol, dmu))
    # compressed CMB (R, ℓ_A)
    DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om) * H0 * DMz / C; lA = np.pi * DMz / rs
    v = np.array([R - R_PL, lA - LA_PL]); o['CMB_RlA'] = float(v @ CMB_CINV @ v)
    # ω_b BBN prior
    o['omega_b_prior'] = float(((ombh2 - OMBH2_PRIOR[0]) / OMBH2_PRIOR[1]) ** 2)
    # SH0ES
    o['SH0ES'] = float(((H0 - SHOES[0]) / SHOES[1]) ** 2)
    # Moresco CC
    zc, H, icc = DATA['cc']; dH = H - bg.hubble_parameter(zc); o['Moresco_CC'] = float(dH @ icc @ dH)
    # fσ8
    zf, fo, fe = DATA['fs8']
    Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz)) / H0)
    o['fsigma8'] = float(np.sum(((fo - fd * s8 * Dd) / fe) ** 2))
    o['TOTAL'] = float(sum(o.values()))
    return o


def fit(model):
    def obj(p):
        Om, H0, ob, MB, s8 = p
        if not (0.2 < Om < 0.45 and 60 < H0 < 78 and 0.019 < ob < 0.025
                and -20.5 < MB < -18.5 and 0.62 < s8 < 0.90):
            return 1e9
        try:
            return chi2_byprobe(p, model)['TOTAL']
        except Exception:
            return 1e9
    r = minimize(obj, [0.305, 68.4, 0.02237, -19.40, 0.78], method='Nelder-Mead',
                 options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000))
    return r.x


def main():
    print("=== T2.3 — per-probe Δχ² decomposition at the canonical joint best-fit ===\n")
    xS = fit('sede'); xL = fit('lcdm')
    bS = chi2_byprobe(xS, 'sede'); bL = chi2_byprobe(xL, 'lcdm')
    probes = ['DESI_BAO', 'Pantheon_SN', 'CMB_RlA', 'omega_b_prior', 'SH0ES', 'Moresco_CC', 'fsigma8', 'TOTAL']
    print(f"{'probe':16s} {'ΛCDM χ²':>10s} {'SEDE χ²':>10s} {'Δχ²':>9s}")
    print("-" * 48)
    deltas = {}
    for k in probes:
        dlt = bS[k] - bL[k]; deltas[k] = dlt
        print(f"{k:16s} {bL[k]:10.2f} {bS[k]:10.2f} {dlt:+9.2f}")
    tot = deltas['TOTAL']
    # which probe(s) carry the preference — test for single AND paired concentration
    contrib = {k: v for k, v in deltas.items() if k != 'TOTAL'}
    neg = {k: v for k, v in contrib.items() if v < 0}
    ranked = sorted(contrib.items(), key=lambda x: x[1])      # most-negative first
    top1, top2 = ranked[0], ranked[1]
    frac1 = top1[1] / tot if tot != 0 else 0
    frac2 = (top1[1] + top2[1]) / tot if tot != 0 else 0
    print("\n  Interpretation:")
    print(f"   total Δχ²(SEDE−ΛCDM) = {tot:+.2f}")
    print(f"   probes favouring SEDE (Δχ²<0): {', '.join(f'{k}({v:+.2f})' for k,v in sorted(neg.items(), key=lambda x:x[1]))}")
    print(f"   largest single  : {top1[0]} ({top1[1]:+.2f} = {100*frac1:.0f}% of total)")
    print(f"   largest pair    : {top1[0]}+{top2[0]} ({top1[1]+top2[1]:+.2f} = {100*frac2:.0f}% of total)")
    if frac1 >= 0.6:
        verdict = f"CONCENTRATED in a single probe ({top1[0]}) — consensus claim must be dropped"
    elif frac2 >= 0.85:
        verdict = (f"CONCENTRATED in the {top1[0]}+{top2[0]} PAIR ({100*frac2:.0f}%) — this is the "
                   f"paired-tension LOO is blind to; the 'multi-probe consensus' framing must be qualified")
    else:
        verdict = "SPREAD across ≥3 probes (genuine multi-probe consensus)"
    print(f"   ⟹ {verdict}.")
    json.dump({'lcdm': bL, 'sede': bS, 'delta': deltas, 'bestfit_sede': list(map(float, xS)),
               'bestfit_lcdm': list(map(float, xL))}, open('results/probe_decomposition.json', 'w'), indent=2)
    print("\n  saved -> results/probe_decomposition.json")


if __name__ == '__main__':
    import os
    os.makedirs('results', exist_ok=True)
    DATA = load()
    main()
