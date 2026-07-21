#!/usr/bin/env python3
"""
E1 — FULL (uncompressed) primary-CMB likelihood test of SEDE  [REVISION_PLAN.md T1.1]
====================================================================================
The headline ΔDIC used the COMPRESSED (R, ℓ_A) Planck priors — exactly the compression
that hid the early–late tension which kills standard HDE (Wu et al. 2509.02945). This
script upgrades to the FULL primary-CMB likelihood, MARGINALISED over the cosmological
parameters — the step run_plik_check.py did not take (it fixed θ* and everything but w).

Likelihood stack (installed under ./packages):
  - planck_2018_highl_plik.TTTEEE_lite_native   (high-ℓ TT/TE/EE peak structure)
  - planck_2018_lowl.TT                          (low-ℓ TT)
  - planck_2018_lowl.EE                          (low-ℓ EE, pins τ)
Marginalised params (BOTH models, same set): ombh2, omch2, cosmomc_theta, logA(→As),
ns, tau, A_planck.  → the early–late tension can hide in correlated shifts of these, so
varying them is the whole point.

SEDE dark sector: w(z) enters CAMB as PPF.  We use a CPL (w0,wa) DERIVED from the SEDE
background (sede.friedmann.E_SEDE_cousin); this is justified by the paper's own result that
SEDE's w(z) is CPL-degenerate to max|w−w_CPL| ≈ 0.027 — far below CMB w-sensitivity —
and the CPL coefficients are recomputed at the best-fit Ω_m as a robustness check.

HONEST SCOPE (state in the paper):
  * plik_lite has foregrounds PRE-marginalised; the full non-lite plik (sampled
    foregrounds) and ACT DR6 PRIMARY are NOT installed → this is plan step (a), the
    full-lite primary + lowl.  ACT DR6 LENSING (Δχ²=−0.32) is run separately
    (run_act_lensing.py); add --lensing-note to print it alongside.
  * The genuinely-untested-before number is the MARGINALISED Δχ² on the primary peaks.

Usage:
  python run_full_cmb_mcmc.py --smoke      # evaluate both models at fiducial (fast)
  python run_full_cmb_mcmc.py --minimize   # profiled Δχ²_min  (the E1 headline; ~30–80 min)
  python run_full_cmb_mcmc.py --mcmc       # full MCMC → ΔDIC   (long; optional)
"""
from __future__ import annotations
import argparse, json, time, sys
import numpy as np

PKG = './packages'
THETA0 = 0.0104109                 # Planck theta*_MC fiducial
LMAX_EXTRA = {'lens_potential_accuracy': 1}

# ----------------------------------------------------------------------------- SEDE w(z)→CPL
def sede_cpl(Om, gamma=1.5, zmax=3.0, n=300):
    """CPL (w0,wa) least-squares fit to the CANONICAL SEDE w(z) over a∈[1/(1+zmax),1].
    Canonical = E_SEDE_volume (Barrow Δ=1 = λ=0.5), γ=1.5 (Result 4C): w0≈−1.0, crosses −1,
    dips to ≈−1.08 by z≈2 (paper §3.2). DE-density-weighted so the fit tracks w where dark
    energy actually matters. NOT E_SEDE_cousin (that is the λ=1 cousin, w0≈−0.86)."""
    from sede.friedmann import E_SEDE_volume
    z = np.linspace(0.0, zmax, n)
    E = E_SEDE_volume(z, Om, gamma)
    rho = np.maximum(E**2 - Om * (1 + z) ** 3, 1e-8)        # ρ_DE/ρ_crit0 (Ω_r negligible here)
    w = -1.0 + (1.0 / 3.0) * (1 + z) * np.gradient(np.log(rho), z)
    a = 1.0 / (1.0 + z)
    wgt = rho / rho[0]                                       # weight by DE density
    # w = w0 + wa (1-a)  → linear LSQ in [1, (1-a)]
    Xc = np.vstack([np.ones_like(a), (1.0 - a)]).T
    W = np.diag(wgt)
    beta = np.linalg.solve(Xc.T @ W @ Xc, Xc.T @ W @ w)
    w0, wa = float(beta[0]), float(beta[1])
    resid = np.max(np.abs(w - (w0 + wa * (1 - a))))
    return w0, wa, resid

# ----------------------------------------------------------------------------- cobaya model
def build_info(model, Om_fid=0.30, gamma=1.5, sampled=True, fullplik=False):  # γ=1.5 canonical (Result 4C)
    """cobaya info for 'lcdm' or 'sede'. sampled=True → priors for minimize/mcmc;
    sampled=False → all fixed at fiducial (smoke). fullplik=True → the full non-lite plik
    (TTTEEE with ~15 foreground nuisance parameters sampled), the gold-standard primary CMB."""
    extra = dict(LMAX_EXTRA)
    de = {}
    meta = {}
    if model == 'sede':
        w0, wa, resid = sede_cpl(Om_fid, gamma)
        extra['dark_energy_model'] = 'ppf'
        de = {'w': w0, 'wa': wa}
        meta = {'w0': w0, 'wa': wa, 'cpl_resid': resid}
    highl = 'planck_2018_highl_plik.TTTEEE' if fullplik else 'planck_2018_highl_plik.TTTEEE_lite_native'
    likelihood = {
        highl: None,
        'planck_2018_lowl.TT': None,
        'planck_2018_lowl.EE': None,
    }
    if sampled:
        params = {
            'ombh2':         {'prior': {'min': 0.019,  'max': 0.025},  'ref': 0.02237, 'proposal': 0.0001},
            'omch2':         {'prior': {'min': 0.10,   'max': 0.14},   'ref': 0.1200,  'proposal': 0.001},
            'cosmomc_theta': {'prior': {'min': 0.0103, 'max': 0.0105}, 'ref': THETA0,  'proposal': 1e-6},
            'logA':          {'prior': {'min': 2.9,    'max': 3.2},    'ref': 3.044,   'proposal': 0.01, 'drop': True},
            'As':            {'value': 'lambda logA: 1e-10*np.exp(logA)', 'derived': False},
            'ns':            {'prior': {'min': 0.92,   'max': 1.00},   'ref': 0.9649,  'proposal': 0.004},
            'tau':           {'prior': {'min': 0.02,   'max': 0.10},   'ref': 0.0544,  'proposal': 0.008},
            'A_planck':      {'prior': {'dist': 'norm', 'loc': 1.0, 'scale': 0.0025}, 'ref': 1.0, 'proposal': 0.002},
            'H0':            {'latex': 'H_0'},        # derived
            'omegam':        {'latex': '\\Omega_m'},  # derived
        }
    else:
        params = {'ombh2': 0.02237, 'omch2': 0.1200, 'cosmomc_theta': THETA0,
                  'logA': {'value': 3.044, 'drop': True},
                  'As': {'value': 'lambda logA: 1e-10*np.exp(logA)', 'derived': False},
                  'ns': 0.9649, 'tau': 0.0544, 'A_planck': 1.0}
    params.update(de)
    info = {'likelihood': likelihood,
            'theory': {'camb': {'extra_args': extra}},
            'params': params, 'packages_path': PKG}
    return info, meta

def chi2_breakdown(model_obj, point):
    """per-likelihood χ² at a parameter point."""
    ll = model_obj.loglikes(point, as_dict=True)[0]
    out = {k: -2.0 * float(v) for k, v in ll.items()}
    out['TOTAL'] = sum(out.values())
    return out


def act_lensing_chi2_at(point, meta, lmax=4000):
    """ACT DR6 CMB-lensing χ² evaluated AT a given (best-fit) cosmology — reuses the
    act_dr6_lenslike machinery of run_act_lensing.py, but at the marginalised primary
    best-fit rather than a fixed fiducial. point: dict of sampled params; meta: SEDE (w0,wa) or {}.
    Returns χ² or None if the package/data are unavailable."""
    try:
        import act_dr6_lenslike as alike
    except Exception:
        return None
    ddir = "packages/data/act_dr6_lens/v1.2/"
    pars = camb.CAMBparams()
    pars.set_cosmology(cosmomc_theta=point['cosmomc_theta'], ombh2=point['ombh2'],
                       omch2=point['omch2'], mnu=0.06, tau=point['tau'])
    As = 1e-10 * np.exp(point['logA'])
    pars.InitPower.set_params(As=As, ns=point['ns'])
    if meta:  # SEDE: PPF w(a) via the CPL coefficients used in the fit
        pars.set_dark_energy(w=meta['w0'], wa=meta['wa'], dark_energy_model='ppf')
    pars.set_for_lmax(lmax, lens_potential_accuracy=4)
    res = camb.get_results(pars)
    pot = res.get_lens_potential_cls(lmax=lmax)
    tot = res.get_cmb_power_spectra(pars, lmax=lmax, CMB_unit='muK', raw_cl=False)['lensed_scalar']
    ell = np.arange(pot.shape[0])
    cl_kk = pot[:, 0] / 4. * 2 * np.pi
    with np.errstate(divide='ignore', invalid='ignore'):
        pf = 2 * np.pi / ell / (ell + 1.)
    pf[:2] = 0.0
    cl_tt, cl_ee, cl_bb, cl_te = tot[:, 0] * pf, tot[:, 1] * pf, tot[:, 2] * pf, tot[:, 3] * pf
    d = alike.load_data('act_baseline', ddir=ddir, lens_only=True, like_corrections=False)
    return -2 * alike.generic_lnlike(d, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb, trim_lmax=2998)


# channels that carry the ISW (late-time DE) signal: low-ℓ TT
ISW_CHANNEL = 'planck_2018_lowl.TT'

import camb  # noqa: E402  (used by act_lensing_chi2_at)

# ----------------------------------------------------------------------------- modes
def run_smoke():
    from cobaya.model import get_model
    print("=== SMOKE: full primary-CMB χ² at the Planck fiducial (no marginalisation) ===")
    res = {}
    for model in ('lcdm', 'sede'):
        info, meta = build_info(model, sampled=False)
        m = get_model(info)
        bd = chi2_breakdown(m, {})
        res[model] = {'chi2': bd, 'meta': meta}
        tag = f"  ({meta.get('w0'):.4f},{meta.get('wa'):.4f}) CPL-resid {meta.get('cpl_resid'):.4f}" if meta else ""
        print(f"\n[{model.upper()}]{tag}")
        for k, v in bd.items():
            print(f"    {k:42s} {v:10.2f}")
    d = res['sede']['chi2']['TOTAL'] - res['lcdm']['chi2']['TOTAL']
    print(f"\n  Δχ²(SEDE−ΛCDM) at fiducial (NOT marginalised) = {d:+.2f}")
    print("  (marginalised number requires --minimize; this only checks the PPF path works)")
    json.dump(res, open('results/e1_smoke.json', 'w'), indent=2)
    print("  saved -> results/e1_smoke.json")
    return res

def run_minimize(fullplik=False):
    from cobaya.model import get_model
    from cobaya.sampler import get_sampler  # noqa
    from cobaya import run as cobaya_run
    tag = "FULL non-lite plik (foregrounds sampled)" if fullplik else "plik_lite"
    print(f"=== MINIMIZE: profiled Δχ²_min on full primary CMB [{tag}] (MARGINALISED) ===")
    out = {}
    for model in ('lcdm', 'sede'):
        info, meta = build_info(model, sampled=True, fullplik=fullplik)
        info['sampler'] = {'minimize': {'method': 'bobyqa', 'max_evals': 1500 if fullplik else 800}}
        info['output'] = f'results/e1_min_{"fp_" if fullplik else ""}{model}'
        info['force'] = True
        t0 = time.time()
        updated, sampler = cobaya_run(info)
        prod = sampler.products()
        mfull = prod['minimum']
        chi2min = -2.0 * float(mfull['minuslogpost'])  # note: includes A_planck prior
        # recompute clean per-likelihood χ² at the minimum
        m = get_model(build_info(model, sampled=True, fullplik=fullplik)[0])
        bestpoint = {p: float(mfull[p]) for p in m.parameterization.sampled_params()}
        bd = chi2_breakdown(m, bestpoint)
        derived = {k: float(mfull[k]) for k in ('H0', 'omegam') if k in mfull}
        lens = act_lensing_chi2_at(bestpoint, meta)   # ACT DR6 lensing AT the primary best-fit
        out[model] = {'chi2_breakdown': bd, 'bestfit': bestpoint, 'derived': derived, 'meta': meta,
                      'act_lensing_chi2': lens, 'minuslogpost_chi2': chi2min,
                      'seconds': time.time() - t0}
        print(f"\n[{model.upper()}] χ²(likes only) = {bd['TOTAL']:.2f}   ({out[model]['seconds']:.0f}s)")
        for k, v in bd.items():
            tag = "   (← low-ℓ ISW channel)" if k == ISW_CHANNEL else ""
            print(f"    {k:42s} {v:10.2f}{tag}")
        if lens is not None:
            print(f"    {'ACT DR6 lensing (at best-fit, not in joint)':42s} {lens:10.2f}")
    dchi2 = out['sede']['chi2_breakdown']['TOTAL'] - out['lcdm']['chi2_breakdown']['TOTAL']
    out['delta_chi2_sede_minus_lcdm'] = dchi2
    # per-channel split
    print("\n=== Δχ²(SEDE−ΛCDM), per channel (marginalised) ===")
    for k in out['lcdm']['chi2_breakdown']:
        dd = out['sede']['chi2_breakdown'][k] - out['lcdm']['chi2_breakdown'][k]
        tag = "   (← low-ℓ ISW channel)" if k == ISW_CHANNEL else ""
        print(f"    {k:42s} {dd:+8.2f}{tag}")
    lL, lS = out['lcdm'].get('act_lensing_chi2'), out['sede'].get('act_lensing_chi2')
    if lL is not None and lS is not None:
        out['delta_chi2_act_lensing_at_bestfit'] = lS - lL
        print(f"    {'ACT DR6 lensing (at best-fit, not in joint)':42s} {lS - lL:+8.2f}")
    print(f"\n  >>> FULL PRIMARY-CMB (marginalised) Δχ²(SEDE−ΛCDM) = {dchi2:+.2f} <<<")
    print("  (negative → SEDE favoured; |Δχ²|<~1 → compressed CMB captured it; the HDE")
    print("   early–late tension would show as a LARGE positive Δχ² here.)")
    # robustness: recompute SEDE (w0,wa) at its best-fit Om vs the fiducial used in the fit
    Om_bf = out['sede'].get('derived', {}).get('omegam')
    if Om_bf:
        w0b, wab, residb = sede_cpl(Om_bf)
        out['sede']['cpl_at_bestfit'] = {'Om': Om_bf, 'w0': w0b, 'wa': wab, 'resid': residb}
        print(f"\n  Robustness — SEDE CPL recomputed at best-fit Ω_m={Om_bf:.4f}: "
              f"(w0,wa)=({w0b:.4f},{wab:.4f}) vs fiducial ({out['sede']['meta']['w0']:.4f},"
              f"{out['sede']['meta']['wa']:.4f}); shift Δw0={w0b-out['sede']['meta']['w0']:+.4f}")
        print("  (small shift ⟹ fixing (w0,wa) at fiducial during the fit is a sound approximation)")
    fn = 'results/e1_minimize_fullplik.json' if fullplik else 'results/e1_minimize.json'
    json.dump(out, open(fn, 'w'), indent=2, default=float)
    print(f"  saved -> {fn}")
    return out

def run_mcmc():
    from cobaya import run as cobaya_run
    print("=== MCMC: ΔDIC on full primary CMB (long) ===")
    for model in ('lcdm', 'sede'):
        info, meta = build_info(model, sampled=True)
        info['sampler'] = {'mcmc': {'max_tries': 1000, 'Rminus1_stop': 0.05, 'burn_in': 200}}
        info['output'] = f'results/e1_mcmc_{model}'
        info['force'] = True
        print(f"\n[{model.upper()}] launching MCMC -> results/e1_mcmc_{model}")
        cobaya_run(info)
    print("\nChains done. Compute DIC from results/e1_mcmc_*.{1.txt} with getdist (see --dic).")

if __name__ == '__main__':
    import os
    os.makedirs('results', exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--minimize', action='store_true')
    ap.add_argument('--mcmc', action='store_true')
    ap.add_argument('--fullplik', action='store_true', help='use the full non-lite plik (foregrounds sampled)')
    a = ap.parse_args()
    if not (a.smoke or a.minimize or a.mcmc):
        a.smoke = True
    if a.smoke:
        run_smoke()
    if a.minimize:
        run_minimize(fullplik=a.fullplik)
    if a.mcmc:
        run_mcmc()
