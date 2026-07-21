#!/usr/bin/env python3
"""
OUT-OF-SAMPLE PREDICTION cross-validation (borrowed from the SEDE_V2 team's CV1/CV2 — the
strongest NON-FLEXIBILITY test).  Pure flexibility helps the in-sample fit but HURTS
held-out prediction; if SEDE predicts data it never saw better than ΛCDM, the preference is
signal, not over-fitting.  Canonical SEDE-H (Barrow λ=0.5, γ=theory) vs ΛCDM; geometry only
(fit Ω_m, H0, ω_b, M_B; r_d from CAMB).  Held-out χ² uses the marginal sub-covariance.

  CV2 redshift-split : train low-z BAO (z<0.8) + SN + CC + priors → predict held-out high-z BAO.
                       (and the reverse, high-z → low-z.)
  CV1 blind-DESI     : train SN + CC + CMB + priors, NO BAO → predict ALL DESI BAO.

Run:  python run_xval_oos.py
"""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve
from sede import friedmann as fr
from sede import data_loader as dl
from run_lambda_verify import w_of_a, data, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES
C = 299792.458
Z, T, M, ICOV = data['desi']; COV = np.linalg.inv(ICOV)
LOW = np.where(Z < 0.8)[0]; HIGH = np.where(Z >= 0.8)[0]

# eBOSS DR16 full-shape (real pre-DESI BAO); keep only geometry rows (DM/rd, DH/rd)
EZ, ET, EM, ECOV = dl.load_eboss_dr16_fs()
if EZ is not None:
    _g = [i for i, tp in enumerate(ET) if tp in ('DM/rd', 'DH/rd')]
    EBZ = EZ[_g]; EBT = [ET[i] for i in _g]; EBM = EM[_g]
    EBICOV = np.linalg.inv(ECOV[np.ix_(_g, _g)])

def eboss_chi2(bg, rd):
    pred = []
    for zz, tp in zip(EBZ, EBT):
        pred.append(bg.comoving_radial_distance(zz)/rd if tp == 'DM/rd' else (C/bg.hubble_parameter(zz))/rd)
    d = EBM - np.array(pred); return float(d @ EBICOV @ d)

def bg_of(Om, H0, ombh2, lcdm):
    h = H0/100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2-ombh2, mnu=0.06)
    if not lcdm:
        a, w = w_of_a(Om, 1.4964, 0.5); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; return camb.get_background(pa)

def desi_pred(bg, rd, idx):
    out = []
    for i in idx:
        zz, tp = Z[i], T[i]
        if tp == 'DM/rd':   out.append(bg.comoving_radial_distance(zz)/rd)
        elif tp == 'DH/rd': out.append((C/bg.hubble_parameter(zz))/rd)
        else:               out.append(((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd)
    return np.array(out)

def chi2_desi_sub(bg, rd, idx):
    """Held-out DESI χ² on a subset, using the MARGINAL sub-covariance."""
    d = M[idx] - desi_pred(bg, rd, idx)
    icov_sub = np.linalg.inv(COV[np.ix_(idx, idx)])
    return float(d @ icov_sub @ d)

def train_chi2(Om, H0, ombh2, MB, lcdm, desi_idx=None, use_cmb=False, use_eboss=False):
    bg = bg_of(Om, H0, ombh2, lcdm); dd = bg.get_derived_params(); rd = dd['rdrag']; c = 0.0
    if desi_idx is not None and len(desi_idx):
        c += chi2_desi_sub(bg, rd, desi_idx)
    if use_eboss:
        c += eboss_chi2(bg, rd)
    zp, mu, chol = data['pan']; dmu = mu - (5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB)
    c += float(dmu @ cho_solve(chol, dmu))
    zc, H, icc = data['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    if use_cmb:
        zs = dd['zstar']; rs = dd['rstar']; DMz = bg.comoving_radial_distance(zs)
        R = np.sqrt(Om)*H0*DMz/C; lA = np.pi*DMz/rs; v = np.array([R-R_PL, lA-LA_PL]); c += float(v @ CMB_CINV @ v)
    c += ((ombh2-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2 + ((H0-SHOES[0])/SHOES[1])**2
    return c

def fit_train(lcdm, desi_idx, use_cmb, use_eboss=False):
    def obj(v):
        Om, H0, ob, MB = v
        if not (0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and -20.5<MB<-18.5): return 1e9
        try: return train_chi2(Om, H0, ob, MB, lcdm, desi_idx, use_cmb, use_eboss)
        except Exception: return 1e9
    return minimize(obj, [0.30, 68.5, 0.02237, -19.40], method='Nelder-Mead',
                    options=dict(xatol=1e-4, fatol=1e-3, maxiter=5000)).x

def heldout(train_idx, test_idx, use_cmb, label, use_eboss=False):
    res = {}
    for lcdm, name in [(False, 'SEDE'), (True, 'ΛCDM')]:
        x = fit_train(lcdm, train_idx, use_cmb, use_eboss)
        bg = bg_of(x[0], x[1], x[2], lcdm); rd = bg.get_derived_params()['rdrag']
        res[name] = chi2_desi_sub(bg, rd, test_idx)
    d = res['SEDE'] - res['ΛCDM']
    print(f"  {label:42s}: held-out χ² SEDE={res['SEDE']:.2f}  ΛCDM={res['ΛCDM']:.2f}  "
          f"Δχ²={d:+.2f}  {'SEDE predicts better' if d<0 else 'ΛCDM better'}")
    return d

if __name__ == "__main__":
    print("="*72); print("OUT-OF-SAMPLE PREDICTION — does SEDE predict data it never saw better than ΛCDM?"); print("="*72)
    print(f"\n[CV2 redshift-split]  (train one half of DESI BAO + SN+CC+priors, predict the other)")
    d_fwd = heldout(LOW, HIGH, False, "low-z (z<0.8) -> predict high-z")
    d_rev = heldout(HIGH, LOW, False, "high-z (z>=0.8) -> predict low-z")
    print(f"\n[CV1 blind-DESI]  (train pre-DESI, NO DESI; predict ALL DESI BAO)")
    d_proxy = heldout(np.array([], dtype=int), np.arange(len(Z)), True, "proxy: SN+CC+CMB -> all DESI BAO")
    d_blind = None
    if EZ is not None:
        d_blind = heldout(np.array([], dtype=int), np.arange(len(Z)), False,
                          "REAL: eBOSS DR16 FS+SN+CC -> all DESI BAO", use_eboss=True)
    print("\n" + "="*72)
    allneg = d_fwd < 0
    blind_str = f" | blind-DESI(eBOSS) {d_blind:+.2f}" if d_blind is not None else ""
    print(f"VERDICT: low->high Δχ²={d_fwd:+.2f} (the DESI evolving-DE direction) | high->low {d_rev:+.2f} | "
          f"blind-DESI(proxy) {d_proxy:+.2f}{blind_str}")
    print(f"  Out-of-sample {'FAVOURS SEDE where it matters (low->high)' if allneg else 'does not favour SEDE'} "
          f"-> {'preference is signal, not just flexibility' if allneg else 'consistent with flexibility'}.")
