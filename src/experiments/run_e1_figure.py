#!/usr/bin/env python3
"""
Fig 9 — E1: the structure-sourcing mechanism, tested on CURRENT data (the honest open frontier)
================================================================================================
Companion plot to the §6 *forecast* (Fig 8). This shows the ACTUAL present-day E1 reconstruction
(`run_e1_mechanism.py`): SEDE's founding claim is that f_sat read off the EXPANSION (geometry, DESI
D_H/r_d → H(z)) equals f_sat read off STRUCTURE GROWTH (RSD fσ8 → σ8(z) → D(z)). The two legs TRACK
at z<1 but DIVERGE at z>1 (structure flattens, geometry declines) — so the distinctive mechanism is
not yet confirmed; the test is RSD-limited and inconclusive (χ²/dof = 0.40 RSD-budget, 1.40
geometry-budget). This is the figure for the paper's most honest sentence: SEDE's EOS is favoured, its
*mechanism* is the open question Euclid/LSST will settle (Fig 8).

Run: python run_e1_figure.py
"""
import numpy as np, camb
from camb import dark_energy
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from sede import friedmann as fr
from sede import data_loader as dl
from run_lambda_verify import w_of_a, data

C = 299792.458
Om, H0, s8, GAM, LAM = 0.298, 68.81, 0.760, 1.4964, 0.5
Or = 9e-5; ODE0 = 1 - Om - Or
Es = lambda z: fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, LAM)
fsat_model = lambda D: (1 - np.exp(-GAM * D ** 2)) / (1 - np.exp(-GAM))


def main():
    h = H0 / 100.; pa = camb.CAMBparams()
    pa.set_cosmology(H0=H0, ombh2=0.02237, omch2=Om * h ** 2 - 0.02237, mnu=0.06)
    a, w = w_of_a(Om, GAM, LAM); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False
    bg = camb.get_background(pa); rd = bg.get_derived_params()['rdrag']

    # geometry leg — DESI D_H/r_d -> H(z) -> f_sat_geom
    z, t, m, icov = data['desi']; cov = np.linalg.inv(icov)
    zg, fg, fg_e = [], [], []
    for i, (zz, tp, val) in enumerate(zip(z, t, m)):
        if tp != 'DH/rd':
            continue
        H = C / (val * rd); E = H / H0
        rhoDE = E ** 2 - Om * (1 + zz) ** 3 - Or * (1 + zz) ** 4
        fsg = rhoDE / (ODE0 * E ** (2 * LAM))
        sig = np.sqrt(cov[i, i]); dE = -(C / (val ** 2 * rd)) / H0 * sig
        drho = 2 * E * dE; dfg = abs(drho / (ODE0 * E ** LAM) - rhoDE * LAM * E ** (LAM - 1) * dE / (ODE0 * E ** (2 * LAM)))
        zg.append(zz); fg.append(fsg); fg_e.append(abs(dfg))
    zg, fg, fg_e = np.array(zg), np.array(fg), np.array(fg_e)

    # structure leg — RSD fσ8 -> σ8(z) -> D(z) -> f_sat_struct
    zf, fo, fe = dl.load_fss8()
    D, frate = fr.compute_growth_model(zf, Om, lambda z: Es(z))
    s8z = fo / frate; Dstr = s8z / s8; fss = fsat_model(Dstr)
    dfdD = (2 * GAM * Dstr * np.exp(-GAM * Dstr ** 2)) / (1 - np.exp(-GAM))
    fss_e = np.abs(dfdD * (fe / frate) / s8)

    # model f_sat curve (the SEDE prediction both legs should follow)
    zc = np.linspace(0.01, 2.0, 200)
    Dm, _ = fr.compute_growth_model(zc, Om, lambda z: Es(z))
    fmod = fsat_model(Dm)

    # significance (as in run_e1_mechanism)
    Dr, _ = fr.compute_growth_model(zf, Om, lambda z: Es(z)); fmod_r = fsat_model(Dr)
    chi2_rsd = float(np.sum(((fss - fmod_r) / fss_e) ** 2)) / len(zf)
    o = np.argsort(zf); fss_g = np.interp(zg, zf[o], fss[o])
    chi2_geo = float(np.sum(((fg - fss_g) / fg_e) ** 2)) / len(zg)
    print(f"χ²/dof: RSD-budget={chi2_rsd:.2f}  geometry-budget={chi2_geo:.2f}")

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ax.axvspan(1.0, 2.0, color='orange', alpha=0.10)
    ax.text(1.5, 0.94, 'z > 1: the legs diverge\n(structure flattens, geometry declines)',
            ha='center', fontsize=8, color='darkorange')
    ax.text(0.45, 0.94, 'z < 1: the legs track', ha='center', fontsize=8, color='green')
    ax.plot(zc, fmod, 'k-', lw=1.5, label='SEDE prediction  f$_{sat}$ = (1−e$^{−γD²}$)/(1−e$^{−γ}$)')
    ax.errorbar(zg, fg, yerr=fg_e, fmt='o', color='C0', capsize=3, ms=6,
                label='GEOMETRY leg  (DESI D$_H$/r$_d$ → H(z))')
    ax.errorbar(zf, fss, yerr=fss_e, fmt='s', color='C3', capsize=2, ms=5, alpha=0.8,
                label='STRUCTURE leg  (RSD fσ8 → σ8(z))')
    ax.set_xlabel('redshift z'); ax.set_ylabel('activation fraction  f$_{sat}$(z)')
    ax.set_xlim(0, 2.0); ax.set_ylim(0, 1.05)
    ax.set_title('E1 — the structure-sourcing mechanism on CURRENT data\n'
                 'SEDE\'s founding claim: the two legs should coincide')
    ax.legend(fontsize=8, loc='lower left')
    ax.text(0.62, 0.30,
            f'χ²/dof = {chi2_rsd:.2f} (RSD budget: consistent)\n'
            f'χ²/dof = {chi2_geo:.2f} (geometry budget: marginal)\n'
            '⟹ inconclusive, RSD-limited\n'
            'Euclid/LSST decides ~4σ (Fig 8; PREREG §F)',
            transform=ax.transAxes, ha='center', va='center', fontsize=7.5,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.92))
    fig.savefig('results/fig_e1_mechanism.png', dpi=130, bbox_inches='tight', pad_inches=0.02)
    fig.savefig('results/fig_e1_mechanism.pdf', bbox_inches='tight', pad_inches=0.02)
    print("  saved -> results/fig_e1_mechanism.png / .pdf")


if __name__ == '__main__':
    import os; os.makedirs('results', exist_ok=True)
    main()
