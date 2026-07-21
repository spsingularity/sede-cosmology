#!/usr/bin/env python3
"""
⚠️ DEPRECATED / SUPERSEDED — uses the λ=1 "temperature-factor" cousin (E_SEDE_cousin via
camb_background, w0≈−0.86, does NOT cross −1), NOT canonical Barrow λ=0.5 (E_SEDE_lambda,
w0≈−1.0, crosses −1). Canonical headline joint: run_lambda_verify.py / run_barrow_mcmc.py
(ΔDIC≈−2.9); full-CMB joint: run_joint_fullcmb.py (Δχ²=−3.17). Kept for the cousin contrast
(paper §3.2) + history; NOT in reproduce_all.py. Do not use for headline results.

PHASE 1: full joint best-fit with CAMB in the loop.
====================================================
Single engine (CAMB background, the λ=1 COUSIN via w(a)) for ALL geometry + sound horizons:
  * BAO with rd = r_drag (CAMB; no free-rd fudge),
  * SN distance moduli from CAMB D_L,
  * compressed CMB (R, l_A) CAMB-exact (Planck-accurate; LCDM -> R -0.7σ, l_A +0.6σ),
  * Moresco H(z), SH0ES H0, fσ8 (smooth-DE growth, exact per Theorem 6).
ombh2 is sampled with a BBN prior; rd is NOT free (it is r_drag).

Compares the λ=1 COUSIN vs LCDM best-fit chi^2. Nelder-Mead (each eval ~0.05 s).
Usage:  python run_camb_joint.py
"""
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import cho_factor, cho_solve

from sede import friedmann as fr
from sede import data_loader as dl
from sede.camb_background import Background
from sede.mcmc_pipeline import build_or_load_table

C = 299792.458
# Planck 2018 distance priors (CHW2019; TT,TE,EE+lowE)
R_PL, LA_PL = 1.7502, 301.471
SR, SLA, RHO = 0.0046, 0.090, 0.46
CMB_CINV = np.linalg.inv(np.array([[SR**2, RHO*SR*SLA], [RHO*SR*SLA, SLA**2]]))
OMBH2_PRIOR = (0.02237, 0.00037)     # BBN-ish
SHOES = (73.04, 1.04)


def load():
    d = {}
    z, t, m, cov = dl.load_desi_dr2(); d['desi'] = (z, t, m, np.linalg.inv(cov))
    zp, mu, covp = dl.load_pantheon_plus()
    d['pan'] = (zp, mu, cho_factor(covp, lower=True)) if zp is not None else None
    zc, H, covc = dl.load_moresco(); d['cc'] = (zc, H, np.linalg.inv(covc))
    d['fs8'] = dl.load_fss8()
    return d


def chi2(bg, data, MB, sigma8, iD, iF, struct):
    # BAO
    z, types, mean, icov = data['desi']
    pred = np.empty(len(z))
    for k, (zz, tp) in enumerate(zip(z, types)):
        if tp == 'DM/rd':   pred[k] = bg.D_M(zz)[0] / bg.r_drag
        elif tp == 'DH/rd': pred[k] = bg.D_H(zz)[0] / bg.r_drag
        else:               pred[k] = bg.D_V(zz)[0] / bg.r_drag
    dd = mean - pred; c = float(dd @ icov @ dd)
    # SN
    if data['pan'] is not None:
        zp, mu, chol = data['pan']
        dmu = mu - (bg.distance_modulus(zp) + MB)
        c += float(dmu @ cho_solve(chol, dmu))
    # CMB (R, l_A)
    v = np.array([bg.R_shift() - R_PL, bg.l_A() - LA_PL]); c += float(v @ CMB_CINV @ v)
    # ombh2 BBN prior, SH0ES
    c += ((bg._ombh2 - OMBH2_PRIOR[0]) / OMBH2_PRIOR[1])**2
    c += ((bg.H0 - SHOES[0]) / SHOES[1])**2
    # Moresco CC
    zc, H, icc = data['cc']; dH = H - bg.hubble(zc); c += float(dH @ icc @ dH)
    # fsigma8 (smooth-DE growth interpolators)
    zf, fo, fe = data['fs8']
    pts = np.column_stack([np.full(len(zf), bg.Om), np.full(len(zf), struct), zf])
    fs8 = iF(pts) * sigma8 * iD(pts)
    c += float(np.sum(((fo - fs8) / fe)**2))
    return c


def make_obj(model, data, iD, iF, struct_is_gamma):
    def obj(p):
        if model == 'cousin':
            Om, H0, ombh2, g, MB, s8 = p
        else:
            Om, H0, ombh2, MB, s8 = p; g = 1.5
        if not (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.019 < ombh2 < 0.025
                and -20.5 < MB < -18.5 and 0.5 < s8 < 1.0 and 0.5 < g < 6):
            return 1e9
        try:
            bg = Background(Om, H0, g, ombh2, model=model)
        except Exception:
            return 1e9
        bg._ombh2 = ombh2
        struct = g if struct_is_gamma else g
        return chi2(bg, data, MB, s8, iD, iF, struct)
    return obj


def fit(model, data, iD, iF, seeds):
    best = None
    for s in seeds:
        r = minimize(make_obj(model, data, iD, iF, True), s, method='Nelder-Mead',
                     options=dict(xatol=1e-4, fatol=1e-3, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return best


def main():
    print("=" * 70)
    print("PHASE 1: full joint with CAMB in the loop — COUSIN (λ=1) vs LCDM")
    print("=" * 70)
    data = load()
    tH = build_or_load_table(model='cousin'); iDH, iFH = fr.make_growth_interpolators(tH)
    tL = build_or_load_table(model='lcdm');  iDL, iFL = fr.make_growth_interpolators(tL)

    print("\nFitting COUSIN (Om,H0,ombh2,gamma,M_B,sigma8)...")
    bH = fit('cousin', data, iDH, iFH,
             [[0.31, 68, 0.02237, g, -19.4, 0.78] for g in (1.0, 1.5)])
    print("Fitting LCDM (Om,H0,ombh2,M_B,sigma8)...")
    bL = fit('lcdm', data, iDL, iFL, [[0.31, 68, 0.02237, -19.4, 0.78]])

    OmH, H0H, obH, gH, MBH, s8H = bH.x
    bgH = Background(OmH, H0H, gH, obH, model='cousin')
    OmL, H0L, obL, MBL, s8L = bL.x
    bgL = Background(OmL, H0L, 1.5, obL, model='lcdm')
    print("\n" + "-" * 70)
    print(f"COUSIN: chi2={bH.fun:.2f}  Om={OmH:.3f} H0={H0H:.2f} 100ωb={obH*100:.3f} "
          f"γ={gH:.2f} σ8={s8H:.3f}")
    print(f"        r_drag={bgH.r_drag:.2f}  R={bgH.R_shift():.4f}({(bgH.R_shift()-R_PL)/SR:+.1f}σ) "
          f"l_A={bgH.l_A():.2f}({(bgH.l_A()-LA_PL)/SLA:+.1f}σ)  100θ*={bgH.theta_star:.5f}")
    print(f"LCDM:   chi2={bL.fun:.2f}  Om={OmL:.3f} H0={H0L:.2f} 100ωb={obL*100:.3f} σ8={s8L:.3f}")
    print(f"        r_drag={bgL.r_drag:.2f}  R={bgL.R_shift():.4f}({(bgL.R_shift()-R_PL)/SR:+.1f}σ) "
          f"l_A={bgL.l_A():.2f}({(bgL.l_A()-LA_PL)/SLA:+.1f}σ)  100θ*={bgL.theta_star:.5f}")
    kH, kL = 6, 5
    print("\n" + "=" * 70)
    print(f"Δχ² (COUSIN − LCDM) = {bH.fun - bL.fun:+.2f}   "
          f"ΔAIC = {(bH.fun+2*kH)-(bL.fun+2*kL):+.2f}   "
          f"({'COUSIN' if bH.fun < bL.fun else 'LCDM'} preferred)")


if __name__ == "__main__":
    main()
