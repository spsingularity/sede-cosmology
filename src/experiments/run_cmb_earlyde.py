#!/usr/bin/env python3
"""
Fig 9 — Why SEDE survives the CMB where standard holographic DE dies
====================================================================================================
The full-CMB result (full non-lite plik Δχ²=−0.33, plik_lite −0.43, ACT DR6 lensing −0.31, all
|Δχ²|<1) is the single highest-leverage result in the paper, and it had no figure. The physical
reason it holds: SEDE carries *less* early dark energy than ΛCDM and its deformation is confined to
the late-time expansion, so the sound horizon, the distance to last scattering, and the acoustic peaks
are essentially unchanged — there is no early–late wall of the kind that disfavours standard HDE.

This figure shows that argument directly:
  Panel A — the early-DE fraction Ω_DE(z): SEDE dips BELOW ΛCDM at high z (gate f_sat→0), while
            un-gated standard Barrow-HDE carries MORE early DE; at z=3, SEDE 0.029 < ΛCDM 0.035.
  Panel B — the fractional expansion deviation H/H_ΛCDM − 1: SEDE's deviation is a few % at z≲1 and
            decays to ~0 by z≈3, flat through BBN/recombination ⟹ the early universe is standard;
            standard Barrow-HDE deviates over a much wider range ⟹ it perturbs the CMB.

Run: python run_cmb_earlyde.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from sede.friedmann import E_SEDE_volume, E_LCDM
from run_delta_orthogonality import E_barrowHDE          # un-gated standard Barrow-HDE solver

Om, GAM, Or = 0.30, 1.4964, 9.0e-5


def OmDE_of(Efun, z):
    E = Efun(z); M = Om * (1 + z) ** 3 + Or * (1 + z) ** 4
    return (E ** 2 - M) / E ** 2


def main():
    z = np.logspace(-1.3, 3.18, 400)                      # z ≈ 0.05 → ~1500 (past recombination)
    E_sede = E_SEDE_volume(z, Om, GAM)
    E_lcdm = E_LCDM(z, Om)
    E_hde = E_barrowHDE(z, 1.0)

    ode_sede = OmDE_of(E_SEDE_volume.__call__ if False else (lambda zz: E_SEDE_volume(zz, Om, GAM)), z)
    ode_lcdm = OmDE_of(lambda zz: E_LCDM(zz, Om), z)
    ode_hde = OmDE_of(lambda zz: E_barrowHDE(zz, 1.0), z)

    # key values
    z3 = 3.0
    s3 = OmDE_of(lambda zz: E_SEDE_volume(zz, Om, GAM), np.array([z3]))[0]
    l3 = OmDE_of(lambda zz: E_LCDM(zz, Om), np.array([z3]))[0]
    print(f"Ω_DE(z=3): SEDE={s3:.3f}  ΛCDM={l3:.3f}  (SEDE has LESS early DE)")

    # Early-DE suppression is MONOTONE in the Barrow deformation Δ (§8.1 third
    # argument for Δ=1): exact identity Ω_DE(z)=Ω_DE0·f_sat·E(z)^(−Δ), E>1 in the
    # past, so larger Δ suppresses early DE more. With the geometric ceiling Δ≤1,
    # the data push Δ to the volume endpoint. (A diagnostic, not a derivation.)
    from sede.friedmann import E_SEDE_barrow
    print("Ω_DE(z=3) vs Barrow Δ (monotone ⇒ Δ=1 maximally CMB-safe):")
    for D in (0.0, 0.5, 1.0):
        oD = OmDE_of(lambda zz, _D=D: E_SEDE_barrow(zz, Om, GAM, Delta=_D), np.array([z3]))[0]
        print(f"   Δ={D:.1f} (λ={1-D/2:.2f}):  Ω_DE(z=3)={oD:.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)

    # Panel A — early-DE fraction
    ax[0].plot(z, ode_hde, 'C3-', lw=2, label='standard Barrow-HDE (un-gated)')
    ax[0].plot(z, ode_lcdm, 'k--', lw=1.6, label='ΛCDM')
    ax[0].plot(z, ode_sede, 'C0-', lw=2, label='SEDE (gated, Δ=1)')
    ax[0].axvspan(1000, 1500, color='gray', alpha=0.15)
    ax[0].text(1090, 0.45, 'recombination', rotation=90, va='center', fontsize=7, color='gray')
    ax[0].plot([z3], [s3], 'C0o', ms=7); ax[0].plot([z3], [l3], 'ks', ms=6, mfc='none')
    ax[0].annotate(f'z=3:  SEDE {s3:.3f} < ΛCDM {l3:.3f}\n(less early DE ⟹ peaks unshifted)',
                   xy=(z3, s3), xytext=(5, 0.30), fontsize=7.5, color='C0',
                   arrowprops=dict(arrowstyle='->', color='C0', lw=1))
    ax[0].set_xscale('log'); ax[0].set_xlabel('redshift z')
    ax[0].set_ylabel('early-DE fraction Ω$_{DE}$(z)')
    ax[0].set_title('SEDE carries less early dark energy than ΛCDM\n(the gate suppresses it; HDE does not)')
    ax[0].legend(fontsize=8, loc='upper right'); ax[0].set_ylim(-0.02, 0.78)

    # Panel B — fractional expansion deviation (the flat SEDE line IS the CMB-safety message)
    dev_sede = (E_sede / E_lcdm - 1) * 100
    dev_hde = (E_hde / E_lcdm - 1) * 100
    sede_max = np.max(np.abs(dev_sede))
    ax[1].axhline(0, color='k', lw=0.8)
    ax[1].plot(z, dev_hde, 'C3-', lw=2, label='standard Barrow-HDE (deviates widely)')
    ax[1].plot(z, dev_sede, 'C0-', lw=2.4, label=f'SEDE (Δ=1): |dev| ≤ {sede_max:.1f}%')
    ax[1].axvspan(1000, 1500, color='gray', alpha=0.15)
    ax[1].text(1090, 8.5, 'recombination', rotation=90, va='center', fontsize=7, color='gray')
    ax[1].annotate('SEDE ≈ ΛCDM at every z\n(sub-% deviation, → 0 by z≈3)\n⟹ sound horizon, peaks &\n'
                   'last-scattering distance\nunchanged: CMB-safe',
                   xy=(2.5, np.interp(2.5, z, dev_sede)), xytext=(0.43, 0.50), textcoords='axes fraction',
                   ha='left', fontsize=7.5, color='C0', arrowprops=dict(arrowstyle='->', color='C0', lw=1))
    ax[1].set_xscale('log'); ax[1].set_xlim(0.05, 1500); ax[1].set_ylim(-1.5, 11.5)
    ax[1].set_xlabel('redshift z')
    ax[1].set_ylabel('expansion deviation  H/H$_{ΛCDM}$ − 1  [%]')
    ax[1].set_title('SEDE deforms only the late expansion → 0 by z≈3\n(BBN/recombination standard)')
    ax[1].legend(fontsize=8, loc='upper right')
    ax[1].text(0.03, 0.97,
               'full plik −0.33 · lite −0.43\nACT DR6 lensing −0.31\nall |Δχ²|<1 ⟹ no early–late wall',
               transform=ax[1].transAxes, fontsize=7, va='top', ha='left',
               bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))

    fig.savefig('results/fig_cmb_earlyde.png', dpi=130, bbox_inches='tight', pad_inches=0.02)
    fig.savefig('results/fig_cmb_earlyde.pdf', bbox_inches='tight', pad_inches=0.02)
    print("  saved -> results/fig_cmb_earlyde.png / .pdf")


if __name__ == '__main__':
    import os; os.makedirs('results', exist_ok=True)
    main()
