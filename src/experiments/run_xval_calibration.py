#!/usr/bin/env python3
"""
#1 FALSE-PREFERENCE CALIBRATION (mock injection) + #6 Bayes factor.

Is ΔDIC≈−2.9 real signal, or the preference a flexible model wins by chance?
  (a) ΛCDM-TRUTH null: generate N mocks from the ΛCDM best-fit (real covariances),
      fit BOTH models to each, build the distribution of Δχ²(SEDE−ΛCDM). Locate the
      REAL −2.96 in that null → a frequentist p-value for the preference.
  (b) SEDE-TRUTH: inject SEDE mocks, confirm SEDE is recovered (detectability).
  (#6) Laplace Bayes factor: with equal dimensionality & priors, ln B ≈ ½Δχ²_min + Occam.

Mocks are FULLY SELF-CONSISTENT: every probe — DESI + SN + CC + fσ8 + CMB (R,l_A) + SH0ES
H0 — is re-drawn from the SAME truth + its covariance (a proper parametric bootstrap), so
the null has no cross-probe inconsistency. SEDE = canonical λ=0.5, γ=theory.
Run: python run_xval_calibration.py (~15 min)
"""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve, cholesky
from sede import friedmann as fr
from run_lambda_verify import w_of_a, data, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES
C = 299792.458
RNG = np.random.default_rng(20260624)
import sys as _sys
# usage: run_xval_calibration.py [n_null] [n_sede]  (defaults 100, 100)
N_NULL = int(_sys.argv[1]) if len(_sys.argv) > 1 else 100
N_SEDE = int(_sys.argv[2]) if len(_sys.argv) > 2 else 100

def bg_of(Om, H0, ombh2, lcdm):
    h = H0/100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2-ombh2, mnu=0.06)
    if not lcdm:
        a, w = w_of_a(Om, 1.4964, 0.5); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; return camb.get_background(pa)

def predict(Om, H0, ombh2, MB, s8, lcdm):
    """Model data vectors (DESI, SN mu, CC H, fσ8) + (R,l_A) for given params."""
    bg = bg_of(Om, H0, ombh2, lcdm); dd = bg.get_derived_params(); rd = dd['rdrag']; zs = dd['zstar']; rs = dd['rstar']
    z, t, _, _ = data['desi']
    desi = np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
        else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz, tp in zip(z, t)])
    zp, _, _ = data['pan']; mu = 5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB
    zc, _, _ = data['cc']; H = bg.hubble_parameter(zc)
    zf, _, _ = data['fs8']; Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz))/H0)
    fs8 = fd*s8*Dd
    R = np.sqrt(Om)*H0*bg.comoving_radial_distance(zs)/C; lA = np.pi*bg.comoving_radial_distance(zs)/rs
    return desi, mu, H, fs8, R, lA

def chi2_mock(Om, H0, ombh2, MB, s8, lcdm, mock):
    d_desi, d_mu, d_H, d_fs8, d_R, d_lA, d_H0 = mock
    desi, mu, H, fs8, R, lA = predict(Om, H0, ombh2, MB, s8, lcdm)
    z, t, m, icov = data['desi']; c = float((d_desi-desi) @ icov @ (d_desi-desi))
    zp, _, chol = data['pan']; c += float((d_mu-mu) @ cho_solve(chol, (d_mu-mu)))
    zc, _, icc = data['cc']; c += float((d_H-H) @ icc @ (d_H-H))
    zf, fo, fe = data['fs8']; c += float(np.sum(((d_fs8-fs8)/fe)**2))
    v = np.array([R-d_R, lA-d_lA]); c += float(v @ CMB_CINV @ v)            # CMB re-drawn from truth
    c += ((ombh2-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2 + ((H0-d_H0)/SHOES[1])**2  # SH0ES re-drawn from truth
    return c

def fit_mock(lcdm, mock):
    def obj(v):
        Om, H0, ob, MB, s8 = v
        if not (0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and -20.5<MB<-18.5 and 0.62<s8<0.90): return 1e9
        try: return chi2_mock(Om, H0, ob, MB, s8, lcdm, mock)
        except Exception: return 1e9
    return minimize(obj, [0.30, 68.5, 0.02237, -19.40, 0.78], method='Nelder-Mead',
                    options=dict(xatol=2e-4, fatol=2e-3, maxiter=4000)).fun

def make_mock(truth, lcdm):
    """Draw a fully self-consistent mock from `truth`: ALL probes (DESI, SN, CC, fσ8,
    CMB R/l_A, SH0ES H0) are re-drawn from the truth + their covariances."""
    desi, mu, H, fs8, R, lA = predict(*truth, lcdm)
    _, _, _, icov = data['desi']; Ld = cholesky(np.linalg.inv(icov), lower=True)
    d_desi = desi + Ld @ RNG.standard_normal(len(desi))
    d_mu = mu + data['_panL'] @ RNG.standard_normal(len(mu))
    _, _, icc = data['cc']; Lc2 = cholesky(np.linalg.inv(icc), lower=True)
    d_H = H + Lc2 @ RNG.standard_normal(len(H))
    _, _, fe = data['fs8']; d_fs8 = fs8 + fe*RNG.standard_normal(len(fs8))
    cmb = np.array([R, lA]) + data['_cmbL'] @ RNG.standard_normal(2)        # CMB from truth
    d_H0 = truth[1] + SHOES[1]*RNG.standard_normal()                        # SH0ES from truth H0
    return (d_desi, d_mu, d_H, d_fs8, cmb[0], cmb[1], d_H0)

if __name__ == "__main__":
    # cache the SN covariance (data['pan'] stores only the cholesky factor)
    import sede.data_loader as dl
    _, _, covp = dl.load_pantheon_plus(); data['_panL'] = cholesky(covp, lower=True)
    data['_cmbL'] = cholesky(np.linalg.inv(CMB_CINV), lower=True)   # CMB (R,l_A) covariance factor
    print("="*72); print("#1 FALSE-PREFERENCE CALIBRATION — is ΔDIC≈−2.9 real or a fluke?"); print("="*72)

    # best-fit truths (from the joint)
    truth_L = (0.299, 68.66, 0.02237, -19.40, 0.760)
    truth_S = (0.298, 68.81, 0.02237, -19.40, 0.760)
    REAL = -4.68   # canonical joint optimiser Δχ²(SEDE−ΛCDM); matches paper §5.1–§5.2

    draws, summary = {}, {"REAL": REAL, "n_null": N_NULL, "n_sede": N_SEDE, "seed": 20260624}
    for label, truth, nmock in [("ΛCDM-TRUTH (null)", truth_L, N_NULL), ("SEDE-TRUTH", truth_S, N_SEDE)]:
        lcdm_truth = label.startswith("ΛCDM")
        d = []
        for i in range(nmock):
            mk = make_mock(truth, lcdm_truth)
            dS = fit_mock(False, mk); dL = fit_mock(True, mk); d.append(dS - dL)
            if (i+1) % 50 == 0: print(f"  [{label}] {i+1}/{nmock} ...", flush=True)
        d = np.array(d)
        draws["null" if lcdm_truth else "sede"] = d
        if lcdm_truth:
            k = int(np.sum(d <= REAL)); p = k / nmock     # fraction at least as SEDE-preferring as real
            # one-sided 95% upper bound when k is small (rule-of-three style)
            ub = 3.0/nmock if k == 0 else (k + 1.6*np.sqrt(k)) / nmock
            # Gaussian read-off of the null (reported alongside the empirical bound, not instead of it)
            z_gauss = (REAL - float(d.mean())) / float(d.std())
            print(f"\n  {label} (N={nmock}): Δχ²(SEDE−ΛCDM) mean={d.mean():+.2f} sd={d.std():.2f} "
                  f"[{np.percentile(d,5):+.2f},{np.percentile(d,95):+.2f}]  min={d.min():+.2f}")
            print(f"    REAL Δχ²={REAL:+.2f} -> {k}/{nmock} beyond it, p-value = {p:.4f} (95% UL {ub:.4f})  "
                  f"({'SIGNAL: real beyond the null' if p<0.05 else 'within null scatter (flexibility)'})")
            print(f"    Gaussian read-off: z = {z_gauss:+.2f}  (~{abs(z_gauss):.1f}sigma)")
            summary.update(null_mean=float(d.mean()), null_sd=float(d.std()), null_min=float(d.min()),
                           null_p5=float(np.percentile(d, 5)), null_p95=float(np.percentile(d, 95)),
                           k_beyond_real=k, p_value=float(p), p_upper_95=float(ub), z_gauss=float(z_gauss))
        else:
            frac = float(np.mean(d < 0))
            print(f"\n  {label}: Δχ²(SEDE−ΛCDM) mean={d.mean():+.2f} sd={d.std():.2f}; "
                  f"SEDE recovered (Δχ²<0) in {100*frac:.0f}% of mocks (detectability)")
            summary.update(sede_mean=float(d.mean()), sede_sd=float(d.std()), sede_recovered_frac=frac)

    # Persist raw draws (for Figure 6 regeneration) and a machine-readable summary (reproducibility).
    import os, json
    resdir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(resdir, exist_ok=True)
    np.savez(os.path.join(resdir, "calibration_draws.npz"), **draws)
    with open(os.path.join(resdir, "calibration_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\n  saved draws -> results/calibration_draws.npz ; summary -> results/calibration_summary.json")

    # #6 Laplace Bayes factor (equal dimensionality & priors -> Occam ≈ 0)
    lnB = 0.5 * abs(REAL)
    print(f"\n[#6 Bayes factor] equal params (5=5) & priors -> ln B ≈ ½Δχ²_min = {lnB:+.2f} "
          f"(B≈{np.exp(lnB):.1f}), 'positive' on the Jeffreys scale, no Occam penalty.")
