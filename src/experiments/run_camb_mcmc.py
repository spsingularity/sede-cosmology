#!/usr/bin/env python3
"""
⚠️ DEPRECATED / SUPERSEDED — uses the λ=1 "temperature-factor" cousin (E_SEDE_cousin via
camb_background / run_camb_joint, w0≈−0.86), NOT canonical Barrow λ=0.5 (E_SEDE_lambda,
w0≈−1.0, crosses −1). Canonical headline: run_barrow_mcmc.py / run_lambda_verify.py
(ΔDIC≈−2.9); full-CMB joint: run_joint_fullcmb.py (Δχ²=−3.17). Kept for the cousin contrast
(paper §3.2) + history; NOT in reproduce_all.py. Do not use for headline results.

PHASE 2: marginalised joint MCMC with CAMB in the loop.
=======================================================
emcee over the CAMB-background joint (sede/camb_background.py): BAO (rd=r_drag),
Pantheon, compressed CMB (R,l_A) CAMB-exact, SH0ES, Moresco, fσ8. Produces
MARGINALISED posteriors + model comparison (mean χ², DIC, min χ²) for SEDE-H vs
LCDM — the proper metric (best-fit χ² overstates tensions). Answers whether the
SEDE-H DESI–CMB concordance tension survives marginalisation.

Usage:  python run_camb_mcmc.py [--steps N] [--burn N] [--walkers N] [--workers N]
"""
import argparse, time
import multiprocessing as mp
import numpy as np

from sede import friedmann as fr
from sede.camb_background import Background
from sede.mcmc_pipeline import build_or_load_table
from run_camb_joint import load, chi2

# module-level state (inherited by forked workers)
DATA = load()
_tH = build_or_load_table(model='cousin'); IDH, IFH = fr.make_growth_interpolators(_tH)
_tL = build_or_load_table(model='lcdm');  IDL, IFL = fr.make_growth_interpolators(_tL)


def logp(theta, model):
    if model == 'cousin':
        Om, H0, ob, g, MB, s8 = theta
        iD, iF = IDH, IFH
    else:
        Om, H0, ob, MB, s8 = theta; g = 1.5
        iD, iF = IDL, IFL
    if not (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.0195 < ob < 0.0250
            and -20.5 < MB < -18.5 and 0.55 < s8 < 0.95 and 0.5 < g < 6):
        return -np.inf
    try:
        bg = Background(Om, H0, g, ob, model=model)
    except Exception:
        return -np.inf
    bg._ombh2 = ob
    return -0.5 * chi2(bg, DATA, MB, s8, iD, iF, g)


def logp_cousin(t): return logp(t, 'cousin')
def logp_lcdm(t):  return logp(t, 'lcdm')


def run(model, lp, p0names, theta0, scales, args):
    import emcee
    ndim = len(theta0)
    rng = np.random.default_rng(7)
    p0 = theta0 + scales * rng.standard_normal((args.walkers, ndim))
    pool = mp.Pool(args.workers) if args.workers > 1 else None
    sampler = emcee.EnsembleSampler(args.walkers, ndim, lp, pool=pool)
    t0 = time.time()
    print(f"[{model}] burn {args.burn}...")
    st = sampler.run_mcmc(p0, args.burn, progress=False)
    sampler.reset()
    print(f"[{model}] prod {args.steps}...")
    sampler.run_mcmc(st, args.steps, progress=False)
    if pool: pool.close()
    flat = sampler.get_chain(flat=True)
    lps = sampler.get_log_prob(flat=True)
    chi2s = -2 * lps
    print(f"[{model}] done {time.time()-t0:.0f}s  acc={np.mean(sampler.acceptance_fraction):.2f}")
    # marginalised summary
    med = np.median(flat, axis=0); lo, hi = np.percentile(flat, [16, 84], axis=0)
    for n, m, l, h in zip(p0names, med, lo, hi):
        print(f"    {n:8s}= {m:.4f}  +{h-m:.4f} -{m-l:.4f}")
    chi2_min = chi2s.min(); chi2_mean = chi2s.mean()
    pD = chi2_mean - chi2_min          # effective #params (approx)
    DIC = chi2_mean + pD
    print(f"    χ²_min={chi2_min:.2f}  <χ²>={chi2_mean:.2f}  p_D≈{pD:.1f}  DIC={DIC:.2f}")
    np.save(f'results/camb_chain_{model}.npy', flat)
    return dict(flat=flat, names=p0names, chi2_min=chi2_min, chi2_mean=chi2_mean, DIC=DIC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=1200)
    ap.add_argument('--burn', type=int, default=400)
    ap.add_argument('--walkers', type=int, default=32)
    ap.add_argument('--workers', type=int, default=max(1, mp.cpu_count() - 1))
    args = ap.parse_args()
    try:
        mp.set_start_method('fork')
    except RuntimeError:
        pass
    print("=" * 70)
    print(f"PHASE 2: marginalised CAMB-in-the-loop MCMC ({args.workers} workers)")
    print("=" * 70)

    rH = run('cousin', logp_cousin,
             ['Om', 'H0', 'ombh2', 'gamma', 'M_B', 'sigma8'],
             np.array([0.30, 68.0, 0.02237, 1.0, -19.40, 0.74]),
             np.array([0.008, 0.4, 0.0003, 0.1, 0.02, 0.015]), args)
    rL = run('lcdm', logp_lcdm,
             ['Om', 'H0', 'ombh2', 'M_B', 'sigma8'],
             np.array([0.30, 68.0, 0.02237, -19.40, 0.74]),
             np.array([0.008, 0.4, 0.0003, 0.02, 0.015]), args)

    print("\n" + "=" * 70 + "\nMARGINALISED MODEL COMPARISON\n" + "=" * 70)
    print(f"  χ²_min : SEDE-H {rH['chi2_min']:.2f}  LCDM {rL['chi2_min']:.2f}  "
          f"Δ={rH['chi2_min']-rL['chi2_min']:+.2f}")
    print(f"  <χ²>   : SEDE-H {rH['chi2_mean']:.2f}  LCDM {rL['chi2_mean']:.2f}  "
          f"Δ={rH['chi2_mean']-rL['chi2_mean']:+.2f}")
    print(f"  DIC    : SEDE-H {rH['DIC']:.2f}  LCDM {rL['DIC']:.2f}  "
          f"ΔDIC={rH['DIC']-rL['DIC']:+.2f}  ({'SEDE-H' if rH['DIC']<rL['DIC'] else 'LCDM'} preferred)")


if __name__ == "__main__":
    main()
