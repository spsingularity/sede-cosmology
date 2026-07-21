#!/usr/bin/env python3
"""
TRACER-LEVEL leave-one-out (§5) — the post-DESI-DR2 robustness test the literature demands.

Recent analyses (e.g. tracer-focused DESI DR2 studies) argue the evolving-dark-energy preference
is driven mainly by the LRG1 (z≈0.51) and LRG2 (z≈0.71) bins, and may be unstable at tracer level.
Our headline leave-one-out (run_xval_loo.py) drops whole *datasets*; this drops individual DESI DR2
*tracer bins* and re-fits the parameter-free SEDE-H (λ=0.5, γ=theory) vs ΛCDM in the same CAMB-pinned
joint. If the preference is real it should survive removing the very bins said to drive it.

Drops: none, −LRG1, −LRG2, −(LRG1+LRG2), and −QSO as a control.  Run:  python run_tracer_loo.py  (~10 min)
"""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve
from sede import friedmann as fr
from sede.data_loader import load_desi_dr2
from run_lambda_verify import w_of_a, data, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES
C = 299792.458

Z_D, T_D, M_D, COV_D = load_desi_dr2()                 # same rows/order as data['desi']
TRACER_Z = {'BGS':0.295, 'LRG1':0.510, 'LRG2':0.706, 'LRG+ELG':0.934,
            'ELG':1.321, 'QSO':1.484, 'Lya':2.330}


def desi_subset(drop_names):
    keep = np.ones(len(Z_D), bool)
    for nm in drop_names:
        keep &= ~np.isclose(Z_D, TRACER_Z[nm], atol=2e-3)
    zs = Z_D[keep]; ts = [T_D[i] for i in range(len(keep)) if keep[i]]; ms = M_D[keep]
    icov = np.linalg.inv(COV_D[np.ix_(keep, keep)])     # re-invert the SLICED covariance (not the inverse)
    return zs, ts, ms, icov


def chi2(Om, H0, ombh2, MB, s8, desi, lcdm=False):
    h = H0/100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2-ombh2, mnu=0.06)
    if not lcdm:
        a, w = w_of_a(Om, 1.4964, 0.5); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; bg = camb.get_background(pa)
    dd = bg.get_derived_params(); rd = dd['rdrag']; zs = dd['zstar']; rs = dd['rstar']; c = 0.0
    z, t, m, icov = desi
    pred = np.array([bg.comoving_radial_distance(zz)/rd if tp == 'DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp == 'DH/rd'
        else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz, tp in zip(z, t)])
    c += float((m-pred) @ icov @ (m-pred))
    zp, mu, chol = data['pan']; dmu = mu - (5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB)
    c += float(dmu @ cho_solve(chol, dmu))
    DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om)*H0*DMz/C; lA = np.pi*DMz/rs
    v = np.array([R-R_PL, lA-LA_PL]); c += float(v @ CMB_CINV @ v)
    c += ((ombh2-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2 + ((H0-SHOES[0])/SHOES[1])**2
    zc, H, icc = data['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    zf, fo, fe = data['fs8']; Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz))/H0)
    c += float(np.sum(((fo - fd*s8*Dd)/fe)**2))
    return c


def fit(desi, lcdm=False):
    def obj(v):
        Om, H0, ob, MB, s8 = v
        if not (0.2 < Om < 0.45 and 60 < H0 < 78 and 0.019 < ob < 0.025 and -20.5 < MB < -18.5 and 0.62 < s8 < 0.90):
            return 1e9
        try: return chi2(Om, H0, ob, MB, s8, desi, lcdm)
        except Exception: return 1e9
    return minimize(obj, [0.30, 68.5, 0.02237, -19.40, 0.78], method='Nelder-Mead',
                    options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000))


if __name__ == "__main__":
    print("=" * 78)
    print("TRACER-LEVEL LEAVE-ONE-OUT — does SEDE survive dropping the LRG bins said to drive DESI?")
    print("=" * 78)
    scenarios = [("none (all 7)", []), ("− LRG1 (z0.51)", ['LRG1']), ("− LRG2 (z0.71)", ['LRG2']),
                 ("− LRG1+LRG2", ['LRG1', 'LRG2']), ("− QSO (control)", ['QSO'])]
    print(f"\n  {'drop':>16s} {'Δχ²(SEDE−ΛCDM)':>16s} {'Ωm':>7s} {'H0':>7s} {'σ8':>7s}  (SEDE best-fit)")
    rows = []
    for tag, drop in scenarios:
        desi = desi_subset(drop)
        rS = fit(desi, lcdm=False); rL = fit(desi, lcdm=True)
        Om, H0, ob, MB, s8 = rS.x; d = rS.fun - rL.fun
        rows.append((tag, d))
        print(f"  {tag:>16s} {d:>16.2f} {Om:>7.3f} {H0:>7.2f} {s8:>7.3f}", flush=True)
    ds = np.array([d for _, d in rows])
    print(f"\n  Δχ² range: [{ds.min():+.2f}, {ds.max():+.2f}]  mean={ds.mean():+.2f}")
    allneg = bool(np.all(ds < 0))
    drop_both = dict(rows)["− LRG1+LRG2"]
    print(f"  SEDE preferred (Δχ²<0) with LRG1+LRG2 BOTH removed: {drop_both < 0}  (Δχ²={drop_both:+.2f})")
    print(f"  SEDE preferred in EVERY tracer holdout: {allneg}  ->",
          "robust at tracer level — NOT an LRG1-2 artefact." if allneg else
          "NOT tracer-robust — the preference rides the LRG bins (report honestly).")
    import sys
    sys.exit(0)
