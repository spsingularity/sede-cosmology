#!/usr/bin/env python3
"""
Model comparisons on the CAMB-in-the-loop joint likelihood.
Adds the two comparisons the ΛCDM-only headline lacked:

  (1) w0waCDM (CPL) — the generic evolving-DE alternative. Gives ΔDIC(SEDE − w0waCDM):
      does SEDE beat *generic* evolving DE, or only ΛCDM?  (7 params: +w0, wa)
  (2) SEDE free-shape — γ and λ are SAMPLED, not frozen at theory. Gives the hidden-
      parameter penalty  DIC(SEDE_free) − DIC(SEDE_frozen): how much of the −4.7 is
      "unearned" because the two structure-gate shape parameters were fixed?  (7 params)

All four models (ΛCDM, SEDE-frozen, w0waCDM, SEDE-free) run on IDENTICAL chains/settings
and the SAME likelihood as the headline (DESI DR2 BAO + Pantheon+ SN + compressed CMB
R/l_A + BBN ω_b prior + SH0ES + cosmic chronometers + fσ8), with the same smooth-DE GR
growth. Marginalised DIC = <χ²> + p_D, p_D = <χ²> − χ²_min (effective # parameters).

Reuses run_lambda_verify (SEDE chi2 + data + constants) and run_barrow_mcmc (ΛCDM +
SEDE-frozen log-probs and MCMC settings), so nothing about the headline likelihood is
re-implemented — only the CPL background and the free-shape priors are added here.

Usage: python run_model_comparison.py [--steps N] [--burn N] [--walkers N] [--workers N]
"""
import argparse, os, time
import multiprocessing as mp
import numpy as np
import camb
from camb import dark_energy
from scipy.linalg import cho_solve
import run_lambda_verify as V   # chi2(Om,H0,ombh2,MB,s8,gamma,lam), w_of_a, data, consts
import run_barrow_mcmc as B     # logp_lcdm, logp_barrow (SEDE-frozen), MCMC pattern

C = V.C
data = V.data


def cpl_table(w0, wa, n=250):
    """CPL w(a) = w0 + wa (1 − a), on the same z-grid convention as V.w_of_a."""
    z = np.linspace(0, 8, n); a = 1.0 / (1.0 + z); w = w0 + wa * (1.0 - a)
    i = np.argsort(a); return a[i], w[i]


def chi2_cpl(Om, H0, ombh2, MB, s8, w0, wa, use_shoes=True):
    """Identical likelihood to V.chi2, with a CPL (w0,wa) background in place of SEDE's."""
    h = H0 / 100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2-ombh2, mnu=0.06)
    a, w = cpl_table(w0, wa); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; bg = camb.get_background(pa)
    d = bg.get_derived_params(); rd = d['rdrag']; zs = d['zstar']; rs = d['rstar']; c = 0.0
    z, t, m, icov = data['desi']
    pred = np.array([bg.comoving_radial_distance(zz)/rd if tp == 'DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp == 'DH/rd'
        else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz, tp in zip(z, t)])
    dd = m - pred; c += float(dd @ icov @ dd)
    zp, mu, chol = data['pan']; dmu = mu - (5*np.log10((1+zp)*bg.comoving_radial_distance(zp)) + 25 + MB); c += float(dmu @ cho_solve(chol, dmu))
    DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om)*H0*DMz/C; lA = np.pi*DMz/rs
    v = np.array([R-V.R_PL, lA-V.LA_PL]); c += float(v @ V.CMB_CINV @ v)
    c += ((ombh2-V.OMBH2_PRIOR[0])/V.OMBH2_PRIOR[1])**2
    if use_shoes: c += ((H0-V.SHOES[0])/V.SHOES[1])**2
    zc, H, icc = data['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    zf, fo, fe = data['fs8']; Dd, fd = V.fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz))/H0)
    c += float(np.sum(((fo - fd*s8*Dd)/fe)**2)); return c


def logp_w0wa(t):
    Om, H0, ob, MB, s8, w0, wa = t
    if not (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.0195 < ob < 0.025 and -20.5 < MB < -18.5
            and 0.62 < s8 < 0.90 and -2.0 < w0 < -0.2 and -3.0 < wa < 2.0):
        return -np.inf
    try:
        c = chi2_cpl(Om, H0, ob, MB, s8, w0, wa)
    except Exception:
        return -np.inf
    return -0.5 * c


def logp_sede_free(t):
    Om, H0, ob, MB, s8, g, lam = t
    if not (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.0195 < ob < 0.025 and -20.5 < MB < -18.5
            and 0.62 < s8 < 0.90 and 0.3 < g < 6.0 and 0.30 < lam < 1.00):
        return -np.inf
    try:
        c = V.chi2(Om, H0, ob, MB, s8, g, lam)[0]
    except Exception:
        return -np.inf
    return -0.5 * c


def run_model(name, lp, theta0, scales, labels, args):
    import emcee
    ndim = len(theta0); rng = np.random.default_rng(args.seed)
    p0 = np.array(theta0) + np.array(scales) * rng.standard_normal((args.walkers, ndim))
    pool = mp.Pool(args.workers) if args.workers > 1 else None
    s = emcee.EnsembleSampler(args.walkers, ndim, lp, pool=pool)
    t0 = time.time(); print(f"[{name}] burn {args.burn}...", flush=True)
    st = s.run_mcmc(p0, args.burn, progress=False); s.reset()
    print(f"[{name}] prod {args.steps}...", flush=True)
    s.run_mcmc(st, args.steps, progress=False)
    if pool: pool.close()
    flat = s.get_chain(flat=True); chi2s = -2 * s.get_log_prob(flat=True)
    ok = np.isfinite(chi2s); flat, chi2s = flat[ok], chi2s[ok]
    print(f"[{name}] done {time.time()-t0:.0f}s acc={np.mean(s.acceptance_fraction):.2f} nsamp={len(chi2s)}", flush=True)
    med = np.median(flat, axis=0); lo, hi = np.percentile(flat, [16, 84], axis=0)
    for nlab, mm, l, h in zip(labels, med, lo, hi):
        print(f"    {nlab:6s}= {mm:.4f} +{h-mm:.4f} -{mm-l:.4f}", flush=True)
    cmin, cmean = float(chi2s.min()), float(chi2s.mean()); pD = cmean - cmin; DIC = cmean + pD
    print(f"    χ²_min={cmin:.2f} <χ²>={cmean:.2f} p_D≈{pD:.1f} DIC={DIC:.2f}", flush=True)
    os.makedirs('results', exist_ok=True); np.save(f'results/model_chain_{name}.npy', flat)
    return dict(cmin=cmin, cmean=cmean, pD=pD, DIC=DIC, med=[float(x) for x in med])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=1200); ap.add_argument('--burn', type=int, default=400)
    ap.add_argument('--walkers', type=int, default=32); ap.add_argument('--workers', type=int, default=7)
    ap.add_argument('--seed', type=int, default=11)
    args = ap.parse_args()
    try: mp.set_start_method('fork')
    except RuntimeError: pass
    print("=" * 72); print("REFEREE COMPARISONS: SEDE vs w0waCDM, and SEDE frozen vs free-shape"); print("=" * 72)
    b5 = ([0.305, 68.4, 0.02237, -19.40, 0.78], [0.008, 0.4, 0.0003, 0.02, 0.015], ['Om', 'H0', 'ombh2', 'M_B', 's8'])
    rL = run_model('lcdm', B.logp_lcdm, *b5, args)
    rS = run_model('sede_frozen', B.logp_barrow, *b5, args)
    rW = run_model('w0wa', logp_w0wa, b5[0] + [-0.95, -0.1], b5[1] + [0.05, 0.2], b5[2] + ['w0', 'wa'], args)
    rF = run_model('sede_free', logp_sede_free, b5[0] + [1.4964, 0.5], b5[1] + [0.3, 0.06], b5[2] + ['gamma', 'lam'], args)
    print("\n" + "=" * 72 + "\nDIC MATRIX\n" + "=" * 72)
    print(f"  ΛCDM          DIC={rL['DIC']:.2f}  p_D={rL['pD']:.1f}")
    print(f"  SEDE(frozen)  DIC={rS['DIC']:.2f}  p_D={rS['pD']:.1f}   ΔDIC(SEDE−ΛCDM)    ={rS['DIC']-rL['DIC']:+.2f}")
    print(f"  w0waCDM       DIC={rW['DIC']:.2f}  p_D={rW['pD']:.1f}   ΔDIC(w0wa−ΛCDM)    ={rW['DIC']-rL['DIC']:+.2f}")
    print(f"  SEDE(free)    DIC={rF['DIC']:.2f}  p_D={rF['pD']:.1f}   ΔDIC(SEDEfree−ΛCDM)={rF['DIC']-rL['DIC']:+.2f}")
    print("-" * 72)
    tag = 'SEDE' if rS['DIC'] < rW['DIC'] else 'w0waCDM'
    print(f"  >>> ΔDIC(SEDE − w0waCDM) = {rS['DIC']-rW['DIC']:+.2f}  ({tag} preferred; SEDE uses {rW['pD']-rS['pD']:+.1f} fewer effective params)")
    print(f"  >>> hidden-param penalty = DIC(SEDEfree) − DIC(SEDEfrozen) = {rF['DIC']-rS['DIC']:+.2f}")
    print(f"  >>> SEDE-free median (γ, λ) = ({rF['med'][5]:.3f}, {rF['med'][6]:.3f})  [theory: 1.496, 0.500]")
    import json
    os.makedirs('results', exist_ok=True)
    json.dump({'lcdm': rL, 'sede_frozen': rS, 'w0wa': rW, 'sede_free': rF,
               'dDIC_SEDE_minus_LCDM': rS['DIC']-rL['DIC'], 'dDIC_w0wa_minus_LCDM': rW['DIC']-rL['DIC'],
               'dDIC_SEDE_minus_w0wa': rS['DIC']-rW['DIC'], 'hidden_param_penalty': rF['DIC']-rS['DIC']},
              open('results/model_comparison.json', 'w'), indent=2)
    print("\n  wrote results/model_comparison.json")


if __name__ == "__main__":
    main()
