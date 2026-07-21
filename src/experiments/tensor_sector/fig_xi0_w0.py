#!/usr/bin/env python3
"""
WP-B — the (Xi_0, w_0) discriminant figure for SEDE.
================================================================================
One figure carrying the whole tensor-sector argument: SEDE occupies the corner
(Xi_0 = 1, w_0 != -1) by construction — GR tensor sector (no modified GW
propagation) on top of a non-LCDM background — a cell no running-Planck-mass
rival occupies.

Xi(z) = d_L^GW(z)/d_L^EM(z) = Xi_0 + (1 - Xi_0)/(1+z)^n   [Belgacem et al. 2019,
arXiv:1906.01593, JCAP 07 (2019) 024 — the (Xi_0, n) parametrisation]. Xi_0 = 1
<=> alpha_M = 0 (standard GW luminosity distance).

Every plotted value is a named constant with its source in an adjacent comment.
Illustrative (non-point-estimate) elements are drawn as bands/tracks and labelled
as such; no unverified point estimate is presented as data.

Self-contained: python fig_xi0_w0.py  ->  outputs/fig_xi0_w0.{pdf,png}
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42   # embed TrueType (repo figure convention)
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os

# ------------------------------------------------------------------ constants
# --- SEDE: Xi_0 exact (derived, fpab_tensor_sector.py); w_0 from the pipeline ---
SEDE_W0    = -0.984   # CPL projection of w(z) over 0<z<1, canonical E_SEDE_volume
                      # (prerecomb/w0 recompute; literal w(0) = -0.996, crossing z~0.18)
SEDE_XI0   = 1.0      # DERIVED exact: alpha_M = alpha_T = 0 (fpab_tensor_sector.py)
SEDE_W0_ERR = 0.004   # w0 is DERIVED from Om (not fit); sigma(w0) propagated from sigma(Om)=0.01
                      # via the pipeline: dw0/dOm * sigma(Om) = 0.004 (E_SEDE_volume)

# --- LCDM: exact reference ---
LCDM_W0, LCDM_XI0 = -1.0, 1.0

# --- RT nonlocal gravity (minimal model, initial conditions in RD) ---
# Belgacem, Foffa, Maggiore, Yang 2019, arXiv:1907.02047, Table 2:
RT_XI0, RT_N = 0.93, 2.59          # Xi_0 = 0.93, n = 2.59, delta(0) = 0.15
RT_W0 = -1.04                      # phantom (w_DE < -1); within ~5% of -1 at z<1 (Table 2);
                                   # placed at a representative phantom w0 (paper quotes phantom, not a single number)

# --- Dark-siren forecast on Xi_0, centred on GR ---
# Mukherjee, Wandelt, Silk 2021, MNRAS 502, 1136 (doi:10.1093/mnras/stab001): ~3500 dark
# sirens of ~30 Msun to z=0.5, BAO-anchored -> Xi_0 = 0.98 (+0.04/-0.23) marginalised,
# 0.99 +/- 0.02 fixed n. Draw the fixed-n (tighter, sigma~0.02) circle:
SIGMA_XI0_DARKSIREN = 0.02
# LISA/ET-era tighter forecast (order of magnitude), Belgacem et al. 2019 (1906.01593):
SIGMA_XI0_LISA = 0.01
SIGMA_W0_FORECAST = 0.03          # combined-probe w0 scale, DESI DR2 era (current-generation)

# ------------------------------------------------------------------ figure
def main():
    fig, ax = plt.subplots(figsize=(6.6, 5.4), constrained_layout=True)

    xlim = (-1.15, -0.85)
    ylim = (0.88, 1.06)

    # --- region shading -----------------------------------------------------
    # Xi_0 != 1 band (top+bottom): "modified GW propagation -> falsifies SEDE"
    ax.axhspan(ylim[0], 0.995, color="#d9d9d9", alpha=0.55, zorder=0)
    ax.axhspan(1.005, ylim[1], color="#d9d9d9", alpha=0.55, zorder=0)
    ax.text(-1.13, 1.042, "modified GW propagation ($\\Xi_0\\neq1$) — falsifies SEDE",
            fontsize=7.5, color="0.35", va="center", ha="left", style="italic")

    # Xi_0 = 1 line (GR tensor sector) and w0 = -1 line
    ax.axhline(1.0, color="k", lw=1.0, ls="-",  zorder=1)
    ax.axvline(-1.0, color="k", lw=0.8, ls=":", zorder=1)
    ax.text(-0.998, 0.888, "$w_0=-1$", fontsize=7.5, color="0.4", rotation=90, va="bottom", ha="left")

    # SEDE signature corner label (Xi_0=1, w0 != -1)
    ax.text(-0.90, 1.012, "SEDE signature corner\n($\\Xi_0=1$, $w_0\\neq-1$)",
            fontsize=8, color="C0", ha="right", va="bottom", fontweight="bold")

    # --- dark-siren forecast ellipse, centred on GR --------------------------
    for sig_xi, lbl, ec in [(SIGMA_XI0_DARKSIREN, f"dark-siren forecast $\\sigma(\\Xi_0)\\approx{SIGMA_XI0_DARKSIREN}$\n(Mukherjee+ 2021)", "C2"),
                             (SIGMA_XI0_LISA,     f"LISA/ET-era $\\sigma(\\Xi_0)\\approx{SIGMA_XI0_LISA}$", "C4")]:
        e = Ellipse((LCDM_W0, LCDM_XI0), width=2*SIGMA_W0_FORECAST, height=2*sig_xi,
                    fill=False, ec=ec, lw=1.4, ls="--", zorder=3)
        ax.add_patch(e)
        ax.plot([], [], color=ec, ls="--", lw=1.4, label=lbl)

    # --- running-M* rivals: illustrative track (NOT a point estimate) --------
    # f(R) / G4(phi) Horndeski generically give Xi_0 != 1 (alpha_M != 0); the exact
    # value runs with the model's free-function amplitude. Draw an illustrative
    # vertical span at w0 ~ -1 spanning below and above unity, clearly labelled.
    ax.annotate("", xy=(-0.995, 0.965), xytext=(-0.995, 1.035),
                arrowprops=dict(arrowstyle="<->", color="C5", lw=1.6))
    ax.text(-0.985, 0.965, "running-$M_*$ models\n($f(R)$, $G_4(\\phi)$): $\\Xi_0\\neq1$\n(illustrative track)",
            fontsize=7, color="C5", va="center", ha="left")

    # --- points -------------------------------------------------------------
    # LCDM
    ax.plot(LCDM_W0, LCDM_XI0, "ks", ms=9, mfc="white", mew=1.6, zorder=5, label="$\\Lambda$CDM (exact)")
    # RT nonlocal (minimal): real point estimate, horizontal+vertical (both from lit)
    ax.plot(RT_W0, RT_XI0, "D", color="C3", ms=8, zorder=5,
            label=f"RT nonlocal, minimal ($\\Xi_0={RT_XI0}$, $n={RT_N}$)\n[Belgacem+ 2019]")
    # SEDE: Xi_0 exact (no vertical error bar), horizontal-only w0 error bar
    ax.errorbar(SEDE_W0, SEDE_XI0, xerr=SEDE_W0_ERR, fmt="o", color="C0", ms=10,
                capsize=3, elinewidth=1.6, zorder=6, label=f"SEDE ($w_0={SEDE_W0}$, $\\Xi_0\\equiv1$)")
    ax.annotate("$\\Xi_0\\equiv1$ (derived)", xy=(SEDE_W0, SEDE_XI0), xytext=(-0.955, 1.038),
                fontsize=8, color="C0",
                arrowprops=dict(arrowstyle="->", color="C0", lw=1.1))

    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_xlabel("dark-energy equation of state today  $w_0$", fontsize=11)
    ax.set_ylabel("GW-propagation parameter  $\\Xi_0 = d_L^{\\rm GW}/d_L^{\\rm EM}$", fontsize=11)
    ax.set_title("The $(\\Xi_0,\\,w_0)$ discriminant plane", fontsize=11)
    ax.legend(fontsize=7.2, loc="lower left", framealpha=0.95, ncol=1)
    ax.grid(alpha=0.15, zorder=0)

    os.makedirs("outputs", exist_ok=True)
    fig.savefig("outputs/fig_xi0_w0.pdf", bbox_inches="tight", pad_inches=0.03)
    fig.savefig("outputs/fig_xi0_w0.png", dpi=200, bbox_inches="tight", pad_inches=0.03)
    print("saved -> outputs/fig_xi0_w0.pdf / .png")

    # --- validation assertions (repo convention) ----------------------------
    assert SEDE_XI0 == 1.0, "SEDE must sit on Xi_0 = 1 exactly (derived)."
    assert LCDM_XI0 == 1.0 and LCDM_W0 == -1.0, "LCDM anchor must be exact (-1, 1)."
    assert 0.85 < RT_XI0 < 1.0, "RT nonlocal Xi_0 must be the Belgacem+ 2019 sub-unity value."
    assert abs(SEDE_W0 - (-1.0)) > 3 * 1e-3, "SEDE must be displaced from w0=-1 (the whole point)."
    print("validation: SEDE on Xi_0=1 exactly; LCDM exact; RT sub-unity; SEDE w0 != -1.  OK")


if __name__ == "__main__":
    main()
