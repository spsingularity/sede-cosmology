#!/usr/bin/env python3
"""
Profiled JOINT fit with the FULL primary CMB in the loop  [REVISION_PLAN.md T1.1, step 1]
=========================================================================================
The headline ΔDIC≈−2.9 used the COMPRESSED (R, ℓ_A) CMB inside the joint. E1
(run_full_cmb_mcmc.py) showed the full primary CMB *alone* does not break SEDE (Δχ²=−0.43).
This closes the loop: it folds the full primary CMB (plik_lite TTTEEE + lowl TT/EE) INTO the
joint, replacing the compressed (R, ℓ_A) term, and profiles the joint χ² for SEDE vs ΛCDM.

Canonical background: E_SEDE_lambda(λ=0.5, γ=1.4964) — the SAME model as run_lambda_verify
(the headline ΔDIC pipeline) and as E1; NOT the E_SEDE_cousin cousin in run_camb_joint.

Late-time terms (DESI BAO, Pantheon+, Moresco CC, fσ8, SH0ES, ω_b BBN prior) are replicated
from run_lambda_verify.chi2 on a CAMB-in-the-loop background (rd = r_drag, no free ruler).
CMB term: full plik_lite + lowl via the cobaya model of run_full_cmb_mcmc (H0-parametrised,
σ8 derived). w(z) → CPL (w0,wa) recomputed per Ω_m (CPL-degenerate to ≈0.02). σ8 is DERIVED
from CAMB (As), not sampled — so the dark sector adds no fitted parameter.

Params (both models): Om, H0, ombh2, M_B, logA, ns, tau, A_planck   (8; γ,λ fixed for SEDE).
Profiled with Nelder-Mead (multi-seed). Reports joint Δχ²_min(SEDE−ΛCDM) with full CMB.

Usage:  python run_joint_fullcmb.py --smoke      # one eval per model (validate)
        python run_joint_fullcmb.py --minimize   # profiled joint fit (the deliverable)
"""
from __future__ import annotations
import argparse, json, time
import numpy as np
import camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve

from sede import friedmann as fr
from run_camb_joint import load, R_PL, LA_PL, CMB_CINV, SHOES, OMBH2_PRIOR  # data + consts
from run_full_cmb_mcmc import sede_cpl

C = 299792.458
GAMMA, LAM = 1.4964, 0.5          # canonical fixed inputs (Result 4C; Barrow Δ=1)
OMNU_H2 = 0.06 / 93.14            # mnu=0.06 eV
LMAX = 2700
DATA = None
_CMB = {}                         # cached cobaya CMB models per model


def w_of_a_canonical(Om, n=250):
    z = np.linspace(0, 8, n)
    E = fr.E_SEDE_lambda(z, Om, GAMMA, LAM)
    rho = np.maximum(E ** 2 - Om * (1 + z) ** 3, 1e-8)
    w = -1 + (1 / 3.) * (1 + z) * np.gradient(np.log(rho), z)
    a = 1 / (1 + z); i = np.argsort(a)
    return a[i], w[i]


def cmb_model(model):
    """Prebuilt cobaya CMB model (full plik_lite + lowl), H0-parametrised, σ8 derived.
    All cosmological params are INPUTS passed per-eval via logposterior()."""
    if model in _CMB:
        return _CMB[model]
    from cobaya.model import get_model
    extra = {'lens_potential_accuracy': 1}
    def rng(a, b, ref):
        return {'prior': {'min': a, 'max': b}, 'ref': ref}
    params = {
        'ombh2': rng(0.018, 0.026, 0.02237), 'omch2': rng(0.08, 0.16, 0.12),
        'H0': rng(55, 80, 68.0), 'As': rng(1.0e-9, 3.5e-9, 2.1e-9),
        'ns': rng(0.90, 1.02, 0.965), 'tau': rng(0.01, 0.12, 0.054),
        'A_planck': rng(0.97, 1.03, 1.0),
        'sigma8': None,   # derived
    }
    if model == 'sede':
        extra['dark_energy_model'] = 'ppf'
        params['w'] = rng(-1.6, -0.4, -1.0)
        params['wa'] = rng(-1.5, 1.5, 0.0)
    info = {'likelihood': {'planck_2018_highl_plik.TTTEEE_lite_native': None,
                           'planck_2018_lowl.TT': None, 'planck_2018_lowl.EE': None},
            'theory': {'camb': {'extra_args': extra}},
            'params': params, 'packages_path': './packages'}
    _CMB[model] = get_model(info)
    return _CMB[model]


def cmb_chi2_and_sigma8(model, Om, H0, ombh2, logA, ns, tau, A_planck):
    """Full primary-CMB χ² and derived σ8 at the given cosmology."""
    h = H0 / 100.
    omch2 = Om * h ** 2 - ombh2                       # match run_lambda_verify convention
    point = {'ombh2': ombh2, 'omch2': omch2, 'H0': H0, 'As': 1e-10 * np.exp(logA),
             'ns': ns, 'tau': tau, 'A_planck': A_planck}
    if model == 'sede':
        w0, wa, _ = sede_cpl(Om)
        point['w'] = w0; point['wa'] = wa
    m = cmb_model(model)
    lp = m.logposterior(point)
    chi2 = -2.0 * float(sum(lp.loglikes))
    sig8 = float(lp.derived[list(m.parameterization.derived_params()).index('sigma8')])
    return chi2, sig8


def latetime_chi2(model, Om, H0, ombh2, MB, sigma8):
    """DESI BAO + Pantheon+ + Moresco + fσ8 + SH0ES + ω_b prior on a CAMB background."""
    h = H0 / 100.
    pa = camb.CAMBparams()
    pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om * h ** 2 - ombh2, mnu=0.06)
    if model == 'sede':
        a, w = w_of_a_canonical(Om)
        de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False
    bg = camb.get_background(pa)
    rd = bg.get_derived_params()['rdrag']
    c = 0.0
    # BAO
    z, t, m, icov = DATA['desi']
    pred = np.array([bg.comoving_radial_distance(zz) / rd if tp == 'DM/rd'
                     else (C / bg.hubble_parameter(zz)) / rd if tp == 'DH/rd'
                     else ((zz * bg.comoving_radial_distance(zz) ** 2 * (C / bg.hubble_parameter(zz))) ** (1 / 3.)) / rd
                     for zz, tp in zip(z, t)])
    dd = m - pred; c += float(dd @ icov @ dd)
    # SN
    zp, mu, chol = DATA['pan']
    dmu = mu - (5 * np.log10((1 + zp) * bg.comoving_radial_distance(zp)) + 25 + MB)
    c += float(dmu @ cho_solve(chol, dmu))
    # ω_b prior, SH0ES
    c += ((ombh2 - OMBH2_PRIOR[0]) / OMBH2_PRIOR[1]) ** 2
    c += ((H0 - SHOES[0]) / SHOES[1]) ** 2
    # Moresco CC
    zc, H, icc = DATA['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    # fσ8 (smooth-DE growth; σ8 DERIVED from CAMB)
    zf, fo, fe = DATA['fs8']
    Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz)) / H0)
    c += float(np.sum(((fo - fd * sigma8 * Dd) / fe) ** 2))
    return c


def joint_chi2(p, model):
    Om, H0, ombh2, MB, logA, ns, tau, A_planck = p
    if not (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.019 < ombh2 < 0.025
            and -20.5 < MB < -18.5 and 2.9 < logA < 3.2 and 0.92 < ns < 1.0
            and 0.02 < tau < 0.10 and 0.98 < A_planck < 1.02):
        return 1e9
    try:
        cmb, sig8 = cmb_chi2_and_sigma8(model, Om, H0, ombh2, logA, ns, tau, A_planck)
        late = latetime_chi2(model, Om, H0, ombh2, MB, sig8)
        # A_planck Gaussian prior (1±0.0025)
        cmb += ((A_planck - 1.0) / 0.0025) ** 2
        return cmb + late
    except Exception as e:
        return 1e9


def fit(model, seeds):
    best = None
    for s in seeds:
        r = minimize(lambda p: joint_chi2(p, model), s, method='Nelder-Mead',
                     options=dict(xatol=1e-4, fatol=1e-3, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return best


SEED = [0.305, 68.4, 0.02237, -19.40, 3.044, 0.965, 0.054, 1.0]


def run_smoke():
    print("=== SMOKE: one joint eval per model (full CMB in the loop) ===")
    for model in ('lcdm', 'sede'):
        t0 = time.time()
        c = joint_chi2(SEED, model)
        print(f"  [{model:4s}] joint χ²(seed) = {c:.2f}   ({time.time()-t0:.1f}s/eval)")


def run_minimize():
    print("=== PROFILED JOINT with FULL primary CMB — SEDE vs ΛCDM ===")
    out = {}
    for model in ('lcdm', 'sede'):
        t0 = time.time()
        b = fit(model, [SEED])
        out[model] = {'chi2': float(b.fun), 'x': list(map(float, b.x)), 'sec': time.time() - t0}
        Om, H0, ob, MB, lA, ns, tau, Ap = b.x
        print(f"\n[{model.upper()}] joint χ²_min = {b.fun:.2f}   ({out[model]['sec']:.0f}s)")
        print(f"   Om={Om:.3f} H0={H0:.2f} 100ωb={ob*100:.3f} M_B={MB:.3f} "
              f"logA={lA:.3f} ns={ns:.4f} τ={tau:.3f} A_planck={Ap:.4f}")
    d = out['sede']['chi2'] - out['lcdm']['chi2']
    out['delta_chi2_sede_minus_lcdm_fullcmb'] = d
    print(f"\n  >>> JOINT (full primary CMB) Δχ²(SEDE−ΛCDM) = {d:+.2f} <<<")
    print(f"  (compressed-CMB joint headline was ≈ −2.9; CMB-only full-primary was −0.43)")
    json.dump(out, open('results/joint_fullcmb.json', 'w'), indent=2)
    print("  saved -> results/joint_fullcmb.json")


if __name__ == '__main__':
    import os
    os.makedirs('results', exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--minimize', action='store_true')
    a = ap.parse_args()
    DATA = load()
    if not (a.smoke or a.minimize):
        a.smoke = True
    if a.smoke:
        run_smoke()
    if a.minimize:
        run_minimize()
