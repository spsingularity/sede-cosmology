#!/usr/bin/env python3
"""
Δ-orthogonality: why SEDE's Δ=1 is NOT excluded by standard Barrow-HDE fits (Δ<0.5)  [R1]
==========================================================================================
Both models share the H-scaling ρ_DE ∝ H^{2−Δ}; the ONLY difference is SEDE's structure gate f_sat(z),
which suppresses dark energy at high z. That gate decouples Δ from the **early-dark-energy fraction** —
the CMB's main lever on Δ:

  standard Barrow-HDE:  ρ_DE ∝ H^{2−Δ}            ⟹ Ω_DE(z=3) is large and Δ-sensitive (0.70→0.15
                                                     over Δ=0→1) ⟹ the CMB/EOS data pin Δ
                                                     (Luciano+ 2025: Δ<0.47–0.54).
  SEDE:                 ρ_DE ∝ H^{2−Δ} f_sat(z)    ⟹ f_sat→0 at high z keeps Ω_DE(z=3) SMALL for all
                                                     Δ (0.13→0.03); at Δ=1, Ω_DE(z=3)=0.029, CMB-safe.

So at the SAME Δ the two models have different observables (at Δ=1: SEDE early-DE 0.029, w₀=−1.0 crossing
−1; Barrow-HDE early-DE 0.147, w₀=−0.77, no crossing). A Δ-constraint derived from Barrow-HDE's
observables therefore does NOT transfer to SEDE: SEDE's Δ=1 is viable precisely where Barrow-HDE's is not.
This script computes Ω_DE(z=3) and w₀ vs Δ for both, and makes results/fig_delta_orthogonality.png.

Run: python run_delta_orthogonality.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from sede.friedmann import E_SEDE_barrow

Om, GAM, Or = 0.30, 1.4964, 9.0e-5
OmDE = 1 - Om - Or


def w0_of_E(Efun):
    z = np.linspace(0, 0.3, 60); E = Efun(z); M = Om*(1+z)**3 + Or*(1+z)**4
    rho = np.maximum(E**2 - M, 1e-9)
    w = -1 + (1/3.)*(1+z)*np.gradient(np.log(rho), z)
    return w[0]


def E_barrowHDE(z, Delta, n=400):
    """Standard Barrow-HDE (Hubble cutoff, NO structure gate): E²=M+Ω_DE0·E^{2−Δ}, solved self-consistently."""
    z = np.atleast_1d(z).astype(float); M = Om*(1+z)**3 + Or*(1+z)**4
    E = np.sqrt(M + OmDE)                              # init
    for _ in range(200):
        En = np.sqrt(M + OmDE * E**(2 - Delta))
        if np.max(np.abs(En - E)) < 1e-10: E = En; break
        E = En
    return E


def OmDE_z(Efun, zq=3.0):
    E = Efun(np.array([zq]))[0]; M = Om*(1+zq)**3 + Or*(1+zq)**4
    return (E**2 - M) / E**2


def main():
    Deltas = np.linspace(0.0, 1.2, 25)
    w0_sede = np.array([w0_of_E(lambda z, D=D: E_SEDE_barrow(z, Om, GAM, Delta=D)) for D in Deltas])
    w0_bhde = np.array([w0_of_E(lambda z, D=D: E_barrowHDE(z, D)) for D in Deltas])
    ede_sede = np.array([OmDE_z(lambda z, D=D: E_SEDE_barrow(z, Om, GAM, Delta=D)) for D in Deltas])
    ede_bhde = np.array([OmDE_z(lambda z, D=D: E_barrowHDE(z, D)) for D in Deltas])

    print("Δ     Ω_DE(z=3) SEDE   Ω_DE(z=3) Barrow-HDE     w0 SEDE   w0 Barrow-HDE")
    for D, es, eb, ws, wb in zip(Deltas, ede_sede, ede_bhde, w0_sede, w0_bhde):
        if abs(D*5 - round(D*5)) < 0.05:
            print(f"  {D:.2f}      {es:.4f}            {eb:.4f}            {ws:+.3f}      {wb:+.3f}")
    print(f"\n  Ω_DE(z=3): at Δ=1, SEDE = {np.interp(1,Deltas,ede_sede):.3f} (CMB-safe) vs "
          f"Barrow-HDE = {np.interp(1,Deltas,ede_bhde):.3f} (5× larger) — the gate.")

    Dh_sede = 1.0     # SEDE: amplitude/early-DE channels prefer the volume endpoint (§5.6, Δ̂≈0.98–1.03)
    Dh_bhde = 0.30    # standard Barrow-HDE: EOS/early-DE fits prefer Δ<0.5 (Luciano+ 2025), Δ=0 within 1σ

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    # Panel A: Ω_DE(z=3) vs Δ — the CMB's main lever on Δ
    ax[0].plot(Deltas, ede_bhde, 'C3-', lw=2, label='standard Barrow-HDE (ρ∝H$^{2-Δ}$, no gate)')
    ax[0].plot(Deltas, ede_sede, 'C0-', lw=2, label='SEDE (ρ∝H$^{2-Δ}$·f$_{sat}$, gated)')
    ax[0].axhline(0.04, color='gray', ls='--', lw=1, label='≈ ΛCDM early-DE (CMB-safe)')
    s1 = np.interp(1, Deltas, ede_sede); b1 = np.interp(1, Deltas, ede_bhde)
    ax[0].plot([1.0], [s1], 'C0o', ms=8)
    ax[0].plot([1.0], [b1], 'C3s', ms=8)
    # annotate that SEDE's Δ=1 is INDEPENDENTLY CMB-safe (≲ ΛCDM early-DE), not merely < Barrow-HDE.
    # Text boxes are placed in the open lower-centre region (below both curves) so neither
    # label overlaps the descending red Barrow-HDE curve or the flat blue SEDE curve.
    ax[0].annotate(f'Barrow-HDE Δ=1: {b1:.3f}\n(5×, CMB-disfavoured)',
                   xy=(1.0, b1), xytext=(0.16, 0.30), textcoords='axes fraction',
                   fontsize=7.5, color='C3', ha='left',
                   bbox=dict(boxstyle='round', fc='white', ec='C3', alpha=0.85),
                   arrowprops=dict(arrowstyle='->', color='C3', lw=1))
    ax[0].annotate(f'SEDE Δ=1: {s1:.3f}\n≲ ΛCDM early-DE\n(independently CMB-safe)',
                   xy=(1.0, s1), xytext=(0.40, 0.11), textcoords='axes fraction',
                   fontsize=7.5, color='C0', ha='left',
                   bbox=dict(boxstyle='round', fc='white', ec='C0', alpha=0.85),
                   arrowprops=dict(arrowstyle='->', color='C0', lw=1))
    ax[0].set_xlabel('Barrow deformation Δ'); ax[0].set_ylabel('early-DE fraction Ω$_{DE}$(z=3)')
    ax[0].set_title('Early-DE fraction vs Δ: the gate holds SEDE low\nwhere un-gated Barrow-HDE is not')
    ax[0].legend(fontsize=8, loc='upper right')
    # Panel B: the data-preferred Δ in each model's own parametrisation
    ax[1].axvspan(0.0, 0.54, color='C3', alpha=0.18)
    ax[1].errorbar([Dh_bhde], [1], xerr=[[0.30], [0.24]], fmt='s', color='C3', capsize=4,
                   label='Barrow-HDE Δ̂ < 0.5 (EOS fit; Luciano+ 2025)')
    ax[1].errorbar([Dh_sede], [0], xerr=[[0.11], [0.11]], fmt='o', color='C0', capsize=4,
                   label='SEDE Δ̂ ≈ 1.0 (amplitude/early-DE; §5.6)')
    ax[1].axvline(1.0, color='C0', ls='--', lw=0.8)
    ax[1].set_xlim(-0.1, 1.3); ax[1].set_ylim(-0.6, 1.6); ax[1].set_yticks([0, 1])
    ax[1].set_yticklabels(['SEDE', 'Barrow-HDE'])
    ax[1].set_xlabel('preferred Δ in each model'); ax[1].set_title('The inferred Δ in each parametrisation\n(the gate is the difference)')
    ax[1].legend(fontsize=8, loc='upper center')
    fig.savefig('results/fig_delta_orthogonality.png', dpi=130, bbox_inches='tight', pad_inches=0.02)
    fig.savefig('results/fig_delta_orthogonality.pdf', bbox_inches='tight', pad_inches=0.02)
    print("\n  ⟹ Same H^{2−Δ} scaling, but the gate decouples Δ from the EOS: the Barrow-HDE EOS")
    print("     constraint (Δ<0.5) does NOT transfer to SEDE, whose Δ is set by the amplitude (→1).")
    print("  saved -> results/fig_delta_orthogonality.png")


if __name__ == '__main__':
    import os; os.makedirs('results', exist_ok=True)
    main()
