#!/usr/bin/env python3
"""
Marginalised CAMB-in-the-loop MCMC for PARAMETER-FREE Barrow SEDE-H.
Δ=1 (maximal Barrow deformation -> λ=0.5) and γ=γ_theory(1.4964) both FIXED, so the
model has the SAME 5 params as ΛCDM: (Ω_m, H0, ω_b, M_B, σ8). Reports marginalised
χ²_min, <χ²>, DIC vs ΛCDM.  Usage: python run_barrow_mcmc.py [--steps N] [--burn N]
"""
import argparse, time
import multiprocessing as mp
import numpy as np
import camb
from scipy.linalg import cho_solve
import run_lambda_verify as V   # provides chi2(Om,H0,ombh2,MB,s8,gamma,lam), data, consts

GAMMA_TH, LAM = 1.4964, 0.5     # Barrow Δ=1 -> λ=0.5; γ at halo-statistics theory
C = V.C; data = V.data

def lcdm_chi2(Om, H0, ombh2, MB, s8):
    h = H0/100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2-ombh2, mnu=0.06)
    pa.WantCls = False; pa.WantTransfer = False; bg = camb.get_background(pa)
    d = bg.get_derived_params(); rd = d['rdrag']; zs = d['zstar']; rs = d['rstar']; c = 0.0
    z, t, m, icov = data['desi']
    pred = np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
        else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz,tp in zip(z,t)])
    dd = m-pred; c += float(dd@icov@dd)
    zp, mu, chol = data['pan']; dmu = mu-(5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB); c += float(dmu@cho_solve(chol, dmu))
    DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om)*H0*DMz/C; lA = np.pi*DMz/rs
    v = np.array([R-V.R_PL, lA-V.LA_PL]); c += float(v@V.CMB_CINV@v)
    c += ((ombh2-V.OMBH2_PRIOR[0])/V.OMBH2_PRIOR[1])**2 + ((H0-V.SHOES[0])/V.SHOES[1])**2
    zc, H, icc = data['cc']; dH = H-bg.hubble_parameter(zc); c += float(dH@icc@dH)
    zf, fo, fe = data['fs8']; Dd, fd = V.fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz))/H0)
    c += float(np.sum(((fo-fd*s8*Dd)/fe)**2)); return c


def logp(theta, model):
    Om, H0, ob, MB, s8 = theta
    if not (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.0195 < ob < 0.025
            and -20.5 < MB < -18.5 and 0.62 < s8 < 0.90):
        return -np.inf
    try:
        c = lcdm_chi2(Om, H0, ob, MB, s8) if model == 'lcdm' else V.chi2(Om, H0, ob, MB, s8, GAMMA_TH, LAM)[0]
    except Exception:
        return -np.inf
    return -0.5 * c

def logp_barrow(t): return logp(t, 'barrow')
def logp_lcdm(t):   return logp(t, 'lcdm')


def run(model, lp, args):
    import emcee
    theta0 = np.array([0.305, 68.4, 0.02237, -19.40, 0.78])
    scales = np.array([0.008, 0.4, 0.0003, 0.02, 0.015])
    rng = np.random.default_rng(11)
    p0 = theta0 + scales * rng.standard_normal((args.walkers, 5))
    pool = mp.Pool(args.workers) if args.workers > 1 else None
    s = emcee.EnsembleSampler(args.walkers, 5, lp, pool=pool)
    t0 = time.time(); print(f"[{model}] burn {args.burn}...")
    st = s.run_mcmc(p0, args.burn, progress=False); s.reset()
    print(f"[{model}] prod {args.steps}...")
    s.run_mcmc(st, args.steps, progress=False)
    if pool: pool.close()
    flat = s.get_chain(flat=True); chi2s = -2*s.get_log_prob(flat=True)
    print(f"[{model}] done {time.time()-t0:.0f}s acc={np.mean(s.acceptance_fraction):.2f}")
    med = np.median(flat, axis=0); lo, hi = np.percentile(flat, [16, 84], axis=0)
    for n, mm, l, h in zip(['Om','H0','ombh2','M_B','s8'], med, lo, hi):
        print(f"    {n:6s}= {mm:.4f} +{h-mm:.4f} -{mm-l:.4f}")
    cmin, cmean = chi2s.min(), chi2s.mean(); pD = cmean-cmin; DIC = cmean+pD
    print(f"    χ²_min={cmin:.2f} <χ²>={cmean:.2f} p_D≈{pD:.1f} DIC={DIC:.2f}")
    np.save(f'results/barrow_chain_{model}.npy', flat)
    return dict(cmin=cmin, cmean=cmean, DIC=DIC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=1500); ap.add_argument('--burn', type=int, default=500)
    ap.add_argument('--walkers', type=int, default=32); ap.add_argument('--workers', type=int, default=7)
    args = ap.parse_args()
    try: mp.set_start_method('fork')
    except RuntimeError: pass
    print("="*68); print(f"PARAMETER-FREE Barrow SEDE-H (Δ=1, λ=0.5, γ=theory) vs ΛCDM — marginalised"); print("="*68)
    rB = run('barrow', logp_barrow, args)
    rL = run('lcdm', logp_lcdm, args)
    print("\n"+"="*68+"\nMARGINALISED COMPARISON\n"+"="*68)
    print(f"  χ²_min : Barrow {rB['cmin']:.2f}  ΛCDM {rL['cmin']:.2f}  Δ={rB['cmin']-rL['cmin']:+.2f}")
    print(f"  <χ²>   : Barrow {rB['cmean']:.2f}  ΛCDM {rL['cmean']:.2f}  Δ={rB['cmean']-rL['cmean']:+.2f}")
    print(f"  DIC    : Barrow {rB['DIC']:.2f}  ΛCDM {rL['DIC']:.2f}  ΔDIC={rB['DIC']-rL['DIC']:+.2f}  "
          f"({'Barrow' if rB['DIC']<rL['DIC'] else 'ΛCDM'} preferred)")


if __name__ == "__main__":
    main()
