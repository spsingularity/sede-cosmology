#!/usr/bin/env python3
"""
No-SH0ES robustness of the SEDE preference
===========================================================================
Round 4 (R6, statistics; seconded by R4, DESI) asked: the per-probe decomposition shows the
joint preference is carried ≈100% by the CMB-distance + SH0ES pair, and SH0ES is the single
most contested prior in the field (DESI's own baseline excludes it). So how much of the
preference survives with SH0ES dropped?

This recomputes the THREE headline statistics with the SH0ES H0 term removed, on the SAME
canonical CAMB-in-the-loop compressed-CMB joint as the headline (so the numbers are directly
comparable):

  (1) joint Δχ²(SEDE−ΛCDM)   — refit both models to the no-SH0ES joint
  (2) nested-sampling ln B    — dynesty, identical 5 priors, loglike = −½·(TOTAL − SH0ES)
  (3) false-preference p      — ΛCDM-truth mocks (all probes re-drawn EXCEPT SH0ES), locate
                                the no-SH0ES REAL Δχ² in the null

If the preference largely survives (CMB-distance carrying it), the result is robust to the H₀
controversy; if it collapses, the honest headline becomes "favoured by the CMB-distance + SH0ES
combination." Either way the reader needs the number.

Run: python run_no_shoes_robustness.py [--nlive 400] [--nmock 300] [--skip-calib]
"""
from __future__ import annotations
import argparse, json
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import cho_solve, cholesky

import run_probe_decomposition as PD
from run_lambda_verify import w_of_a, data as LV_DATA, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES

C = 299792.458
RNG = np.random.default_rng(20260628)
X0 = [0.305, 68.4, 0.02237, -19.40, 0.78]
BOUNDS = lambda Om, H0, ob, MB, s8: (0.2 < Om < 0.45 and 60 < H0 < 78 and 0.019 < ob < 0.025
                                     and -20.5 < MB < -18.5 and 0.62 < s8 < 0.90)


# ───────────────────────── (1) joint Δχ² with SH0ES removed ─────────────────────────
def total_noshoes(p, model):
    o = PD.chi2_byprobe(p, model)
    return o['TOTAL'] - o['SH0ES'], o          # drop the SH0ES datum from the sum


def fit_noshoes(model):
    def obj(p):
        if not BOUNDS(*p):
            return 1e9
        try:
            return total_noshoes(p, model)[0]
        except Exception:
            return 1e9
    r = minimize(obj, X0, method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000))
    return r.x


def part1_joint():
    print("=== (1) joint Δχ²(SEDE−ΛCDM) with SH0ES dropped ===\n")
    xS, xL = fit_noshoes('sede'), fit_noshoes('lcdm')
    tS, bS = total_noshoes(xS, 'sede'); tL, bL = total_noshoes(xL, 'lcdm')
    probes = ['DESI_BAO', 'Pantheon_SN', 'CMB_RlA', 'omega_b_prior', 'Moresco_CC', 'fsigma8']
    print(f"  {'probe':16s} {'ΛCDM χ²':>10s} {'SEDE χ²':>10s} {'Δχ²':>9s}")
    print("  " + "-" * 48)
    for k in probes:
        print(f"  {k:16s} {bL[k]:10.2f} {bS[k]:10.2f} {bS[k]-bL[k]:+9.2f}")
    dchi2 = tS - tL
    print(f"  {'JOINT (no SH0ES)':16s} {tL:10.2f} {tS:10.2f} {dchi2:+9.2f}")
    print(f"\n  >>> no-SH0ES joint Δχ²(SEDE−ΛCDM) = {dchi2:+.2f}  (headline WITH SH0ES = −2.95) <<<")
    print(f"      ΛCDM best-fit H0 = {xL[1]:.2f}, SEDE best-fit H0 = {xS[1]:.2f} "
          f"(no longer pulled toward 73.04)")
    return dchi2, tuple(xL), tuple(xS)


# ───────────────────────── (2) nested-sampling ln B, no SH0ES ─────────────────────────
PRIORS_LO = np.array([0.20, 60.0, 0.019, -20.5, 0.62])
PRIORS_HI = np.array([0.45, 78.0, 0.025, -18.5, 0.90])


def part2_evidence(nlive):
    import dynesty
    print(f"\n=== (2) nested-sampling ln B (dynesty, nlive={nlive}), SH0ES dropped ===\n")

    def ptform(u):
        return PRIORS_LO + u * (PRIORS_HI - PRIORS_LO)

    def make_loglike(model):
        def ll(theta):
            try:
                return -0.5 * total_noshoes(theta, model)[0]
            except Exception:
                return -1e10
        return ll

    lnZ = {}
    for model in ('lcdm', 'sede'):
        s = dynesty.NestedSampler(make_loglike(model), ptform, ndim=5, nlive=nlive,
                                  bound='multi', sample='rwalk')
        s.run_nested(print_progress=False, dlogz=0.1)
        lnZ[model] = (float(s.results.logz[-1]), float(s.results.logzerr[-1]))
        print(f"  {model.upper():5s}  ln Z = {lnZ[model][0]:8.2f} ± {lnZ[model][1]:.2f}")
    lnB = lnZ['sede'][0] - lnZ['lcdm'][0]
    err = float(np.hypot(lnZ['sede'][1], lnZ['lcdm'][1]))
    a = abs(lnB)
    scale = ("inconclusive (|lnB|<1)" if a < 1 else "positive/weak (1–2.5)" if a < 2.5
             else "strong (2.5–5)" if a < 5 else "very strong (>5)")
    print(f"\n  >>> no-SH0ES ln B (SEDE − ΛCDM) = {lnB:+.2f} ± {err:.2f} <<<  "
          f"[{scale}, favouring {'SEDE' if lnB > 0 else 'ΛCDM'}]  (headline WITH SH0ES = +1.89)")
    return lnB, err, lnZ


# ───────────────────────── (3) false-preference calibration, no SH0ES ─────────────────────────
def predict(Om, H0, ombh2, MB, s8, lcdm):
    import camb
    from camb import dark_energy
    h = H0 / 100.; pa = camb.CAMBparams()
    pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om * h ** 2 - ombh2, mnu=0.06)
    if not lcdm:
        a, w = w_of_a(Om, 1.4964, 0.5); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False
    bg = camb.get_background(pa); dd = bg.get_derived_params()
    rd, zs, rs = dd['rdrag'], dd['zstar'], dd['rstar']
    from sede import friedmann as fr
    z, t, _, _ = LV_DATA['desi']
    desi = np.array([bg.comoving_radial_distance(zz) / rd if tp == 'DM/rd'
                     else (C / bg.hubble_parameter(zz)) / rd if tp == 'DH/rd'
                     else ((zz * bg.comoving_radial_distance(zz) ** 2 * (C / bg.hubble_parameter(zz))) ** (1 / 3.)) / rd
                     for zz, tp in zip(z, t)])
    zp, _, _ = LV_DATA['pan']; mu = 5 * np.log10((1 + zp) * bg.comoving_radial_distance(zp)) + 25 + MB
    zc, _, _ = LV_DATA['cc']; H = bg.hubble_parameter(zc)
    zf, _, _ = LV_DATA['fs8']; Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz)) / H0)
    fs8 = fd * s8 * Dd
    R = np.sqrt(Om) * H0 * bg.comoving_radial_distance(zs) / C; lA = np.pi * bg.comoving_radial_distance(zs) / rs
    return desi, mu, H, fs8, R, lA


def chi2_mock_ns(Om, H0, ombh2, MB, s8, lcdm, mock):
    """no-SH0ES mock χ²: identical to the calibration but WITHOUT the SH0ES H0 term."""
    d_desi, d_mu, d_H, d_fs8, d_R, d_lA = mock
    desi, mu, H, fs8, R, lA = predict(Om, H0, ombh2, MB, s8, lcdm)
    _, _, _, icov = LV_DATA['desi']; c = float((d_desi - desi) @ icov @ (d_desi - desi))
    _, _, chol = LV_DATA['pan']; c += float((d_mu - mu) @ cho_solve(chol, (d_mu - mu)))
    _, _, icc = LV_DATA['cc']; c += float((d_H - H) @ icc @ (d_H - H))
    _, fo, fe = LV_DATA['fs8']; c += float(np.sum(((d_fs8 - fs8) / fe) ** 2))
    v = np.array([R - d_R, lA - d_lA]); c += float(v @ CMB_CINV @ v)
    c += ((ombh2 - OMBH2_PRIOR[0]) / OMBH2_PRIOR[1]) ** 2          # ω_b kept; SH0ES removed
    return c


def fit_mock_ns(lcdm, mock):
    def obj(v):
        if not BOUNDS(*v):
            return 1e9
        try:
            return chi2_mock_ns(*v, lcdm, mock)
        except Exception:
            return 1e9
    return minimize(obj, [0.30, 68.5, 0.02237, -19.40, 0.78], method='Nelder-Mead',
                    options=dict(xatol=2e-4, fatol=2e-3, maxiter=4000)).fun


def make_mock_ns(truth, lcdm):
    desi, mu, H, fs8, R, lA = predict(*truth, lcdm)
    _, _, _, icov = LV_DATA['desi']; Ld = cholesky(np.linalg.inv(icov), lower=True)
    d_desi = desi + Ld @ RNG.standard_normal(len(desi))
    d_mu = mu + LV_DATA['_panL'] @ RNG.standard_normal(len(mu))
    _, _, icc = LV_DATA['cc']; Lc = cholesky(np.linalg.inv(icc), lower=True)
    d_H = H + Lc @ RNG.standard_normal(len(H))
    _, _, fe = LV_DATA['fs8']; d_fs8 = fs8 + fe * RNG.standard_normal(len(fs8))
    cmb = np.array([R, lA]) + LV_DATA['_cmbL'] @ RNG.standard_normal(2)
    return (d_desi, d_mu, d_H, d_fs8, cmb[0], cmb[1])           # no d_H0


def part3_calibration(real_dchi2, truth_L, nmock):
    import sede.data_loader as dl
    _, _, covp = dl.load_pantheon_plus(); LV_DATA['_panL'] = cholesky(covp, lower=True)
    LV_DATA['_cmbL'] = cholesky(np.linalg.inv(CMB_CINV), lower=True)
    print(f"\n=== (3) false-preference calibration, SH0ES dropped (N={nmock} ΛCDM-truth mocks) ===\n")
    d = []
    for i in range(nmock):
        mk = make_mock_ns(truth_L, True)
        d.append(fit_mock_ns(False, mk) - fit_mock_ns(True, mk))
        if (i + 1) % 50 == 0:
            print(f"  [null] {i+1}/{nmock} ...", flush=True)
    d = np.array(d)
    k = int(np.sum(d <= real_dchi2)); p = k / nmock
    ub = 3.0 / nmock if k == 0 else (k + 1.6 * np.sqrt(k)) / nmock
    print(f"\n  null Δχ²(SEDE−ΛCDM): mean={d.mean():+.2f} sd={d.std():.2f} "
          f"[{np.percentile(d,5):+.2f},{np.percentile(d,95):+.2f}] min={d.min():+.2f}")
    print(f"  no-SH0ES REAL Δχ²={real_dchi2:+.2f} -> {k}/{nmock} beyond it, "
          f"p = {p:.4f} (95% UL {ub:.4f})  "
          f"[{'SIGNAL beyond the null' if p < 0.05 else 'within null scatter'}]")
    print(f"  (headline WITH SH0ES: null mean +0.42, p = 0.004)")
    return {'p': p, 'k': k, 'nmock': nmock, 'null_mean': float(d.mean()),
            'null_sd': float(d.std()), 'real': real_dchi2, 'ub': ub}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nlive', type=int, default=400)
    ap.add_argument('--nmock', type=int, default=300)
    ap.add_argument('--skip-calib', action='store_true')
    a = ap.parse_args()
    PD.DATA = PD.load()
    print("=" * 72)
    print("NO-SH0ES ROBUSTNESS OF THE SEDE PREFERENCE  (round-4 R6+R4)")
    print("=" * 72 + "\n")
    out = {}
    dchi2, truth_L5, truth_S5 = part1_joint()
    out['joint_dchi2_noshoes'] = dchi2
    out['bestfit_lcdm'] = list(map(float, truth_L5))
    out['bestfit_sede'] = list(map(float, truth_S5))
    lnB, lnBerr, lnZ = part2_evidence(a.nlive)
    out['lnB_noshoes'] = lnB; out['lnB_err'] = lnBerr; out['lnZ'] = lnZ
    if not a.skip_calib:
        out['calibration'] = part3_calibration(dchi2, truth_L5, a.nmock)
    json.dump(out, open('results/no_shoes_robustness.json', 'w'), indent=2)
    print("\n  saved -> results/no_shoes_robustness.json")
    print("\n" + "=" * 72)
    print("SUMMARY (no SH0ES):  joint Δχ² = %+.2f   ln B = %+.2f   %s" % (
        dchi2, lnB, ("p = %.4f" % out['calibration']['p']) if not a.skip_calib else "(calib skipped)"))
    print("HEADLINE (w/ SH0ES): joint Δχ² = −2.95   ln B = +1.89   p = 0.004")
    print("=" * 72)


if __name__ == '__main__':
    import os
    os.makedirs('results', exist_ok=True)
    main()
