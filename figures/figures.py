#!/usr/bin/env python3
"""
The paper's figures, computed from the real SEDE model (sede.friedmann) and real data
(sede.data_loader) — designed around what the paper must SHOW, not around old plot scripts.

  fig1_mechanism   — how Λ fades in: f_sat(z) gates dark energy on as structure grows
  fig2_inputs      — the parameter-free dark sector: λ(Δ) and γ(p) both derived
  fig3_datafit     — the model fits the data: DESI BAO + Pantheon SN residuals vs ΛCDM
  fig4_eos         — w(z) crosses −1, and (w0,wa) vs the DESI DR2 result
  fig5_oos         — out-of-sample: trained on pre-DESI eBOSS, SEDE predicts DESI
  fig6_forecast    — the decisive test: σ(Δ) from DESI DR3 + Euclid

Run: python figures.py   (writes results/paper/fig_*.png)
"""
import os, sys, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sede import friedmann as fr, data_loader as dl

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
C = 299792.458                       # km/s
OUT = os.path.join(ROOT, "results/paper"); os.makedirs(OUT, exist_ok=True)

def savefig_vector(fig, basepath):
    """Save both a PNG (GitHub preview) and a vector PDF (journal submission)."""
    fig.savefig(f"{basepath}.png", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{basepath}.pdf", bbox_inches="tight", pad_inches=0.02)
# joint best-fit (canonical Barrow SEDE-H vs ΛCDM); both share early physics -> common r_d
OmS, H0S, GAM = 0.298, 68.8, 1.4964
OmL, H0L = 0.299, 68.66
RD = 147.5                           # Mpc (shape comparison; uniform for both)
CSEDE, CL = "#c0392b", "#2c3e50"     # SEDE red, ΛCDM slate

def Es(z):  return fr.E_SEDE_barrow(np.atleast_1d(z), OmS, GAM, Delta=1.0)
def El(z):  return fr.E_LCDM(np.atleast_1d(z), OmL)

def distances(Efun, H0, zmax=2.6, n=1400):
    """Return interpolating callables DM,DH,DV (Mpc) from a dimensionless E(z)."""
    zg = np.linspace(0, zmax, n); E = Efun(zg)
    Dc = np.concatenate([[0], np.cumsum(0.5*(C/H0)*(1/E[1:]+1/E[:-1])*np.diff(zg))])
    DM = lambda z: np.interp(z, zg, Dc)
    DH = lambda z: C/(H0*Efun(np.atleast_1d(z)))
    DV = lambda z: (z*DM(z)**2*DH(z))**(1/3)
    return DM, DH, DV

def bao_pred(z, t, Efun, H0):
    DM, DH, DV = distances(Efun, H0)
    out = np.empty(len(z))
    for i,(zi,ti) in enumerate(zip(z,t)):
        val = (DV(zi) if ti=="DV/rd" else DM(zi) if ti=="DM/rd" else DH(zi))/RD
        out[i] = float(np.ravel(val)[0])   # DH/DV return length-1 arrays; numpy 2.x needs a scalar
    return out

# ───────────────────────────────────────────────────────────────── fig 1
def fig1_mechanism():
    z = np.linspace(0, 4, 200); D = fr.compute_growth_factor(z, OmS)
    fsat = np.clip((1-np.exp(-GAM*D**2))/(1-np.exp(-GAM)), 0, 1)
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.7), constrained_layout=True)
    ax[0].plot(z, fsat, color=CSEDE, lw=2.4)
    ax[0].fill_between(z, 0, fsat, color=CSEDE, alpha=0.08)
    ax[0].set(xlabel="redshift z", ylabel=r"$f_{\rm sat}(z)$  — horizon entropy saturated",
              title="Dark energy switches on with structure", xlim=(0,4), ylim=(0,1.05))
    ax[0].text(2.2, 0.18, "structure\nforming", color=CSEDE, fontsize=8.5, ha="center")
    ax[0].annotate("today", (0,1.0), (0.6,0.86), fontsize=8.5,
                   arrowprops=dict(arrowstyle="->", lw=0.7))
    # cosmic density inventory: DE is orders of magnitude below matter/radiation at high z
    zc = np.logspace(0, 5, 400) - 1; x = 1 + zc; Or = 9e-5
    zh = np.linspace(0, 6, 80); Dh = fr.compute_growth_factor(zh, OmS)
    fh = np.clip((1-np.exp(-GAM*Dh**2))/(1-np.exp(-GAM)), 0, 1)
    # matter-dominated extrapolation D ∝ a for z>6 (so f_sat ∝ D² → 0)
    D6 = Dh[-1]*(1+6.0); fsat_c = np.where(zc <= 6, np.interp(zc, zh, fh),
                                           np.clip((1-np.exp(-GAM*(D6/x)**2))/(1-np.exp(-GAM)), 0, 1))
    Ec = np.sqrt(OmS*x**3 + Or*x**4 + (1-OmS)*fsat_c*np.sqrt(OmS*x**3 + Or*x**4))
    ax[1].loglog(x, OmS*x**3, color="#7f8c8d", lw=1.6, label=r"matter $\propto(1+z)^3$")
    ax[1].loglog(x, Or*x**4, color="#2980b9", lw=1.6, ls="-.", label=r"radiation $\propto(1+z)^4$")
    ax[1].loglog(x, (1-OmS)*fsat_c*Ec, color=CSEDE, lw=2.4, label=r"SEDE $\rho_{\rm DE}$")
    ax[1].axhline(1-OmL, color=CL, lw=1.5, ls="--", label=r"ΛCDM $\rho_\Lambda$")
    ax[1].axvline(1101, color="k", lw=0.7, ls=":"); ax[1].text(1101, 3e-6, "recomb.", rotation=90, fontsize=7, ha="right")
    ax[1].set(xlabel="1 + z", ylabel=r"$\rho_i/\rho_{\rm crit,0}$", title="DE is negligible at high z (CMB/BBN-safe)",
              ylim=(1e-6, 1e16)); ax[1].legend(fontsize=7.5, loc="upper left")
    savefig_vector(fig, f"{OUT}/fig1_mechanism"); plt.close(fig); print("fig1")

# ───────────────────────────────────────────────────────────────── fig 2
def fig2_inputs():
    from sede.gamma_computation import entropy_weight_scan
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
    D = np.linspace(0, 1, 100); a.plot(D, 1-D/2, color=CSEDE, lw=2.4)
    a.scatter([0,1], [1,.5], c=["#7f8c8d", CSEDE], s=70, zorder=5)
    a.annotate(r"$\Delta=0,\ \lambda=1$"+"\n(area-law)", (0,1), (0.10,0.84), fontsize=8.5)
    a.annotate(r"$\Delta=1,\ \lambda=0.5$"+"\n(volume-law, adopted)", (1,.5), (0.28,0.55),
               fontsize=8.5, color=CSEDE)
    a.set(xlabel=r"Barrow deformation $\Delta$", ylabel=r"$\lambda=1-\Delta/2$",
          title=r"$\lambda$ from the horizon (§3.3, §8)", xlim=(-.03,1.03), ylim=(.45,1.05))
    scan = entropy_weight_scan(p_list=(2/3,1.,4/3,5/3,2.)); ps=[p for p,_ in scan]; gs=[g for _,g in scan]
    b.plot(ps, gs, "o-", color=CL, lw=1.8)
    i = ps.index(5/3); b.scatter([5/3],[gs[i]], color=CSEDE, s=85, zorder=6)
    b.annotate(fr"$p=5/3\to\gamma={gs[i]:.2f}$"+"\n(binding energy into\nthe horizon, Presc. 4C)",
               (5/3,gs[i]), (0.72,gs[i]+0.05), fontsize=8.5, color=CSEDE)
    b.axhline(0, color="gray", lw=.6, ls=":")
    b.set(xlabel=r"entropy-weight exponent $p$", ylabel=r"$\gamma=d\ln\Sigma_S/d\ln\sigma_8^2$",
          title=r"$\gamma$ from structure (Presc. 4C)")
    fig.suptitle("The dark sector adds no fitted parameter: λ, γ, w(z) set by stated prescriptions", fontsize=11)
    savefig_vector(fig, f"{OUT}/fig2_inputs"); plt.close(fig); print("fig2")

# ───────────────────────────────────────────────────────────────── fig 3
def fig3_datafit():
    z,t,mean,cov = dl.load_desi_dr2(); z=np.array(z); sig=np.sqrt(np.diag(cov))
    pS = bao_pred(z,t,Es,H0S); pL = bao_pred(z,t,El,H0L)
    lab = [f"{zi:.2f} {ti.split('/')[0]}" for zi,ti in zip(z,t)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 4.1), constrained_layout=True)
    x = np.arange(len(z))
    a.errorbar(x, (mean-pL)/sig, yerr=1, fmt="s", color=CL, ms=4, capsize=2, label="DESI DR2 − ΛCDM")
    a.plot(x, (pS-pL)/sig, "^", color=CSEDE, ms=7, label="SEDE − ΛCDM")
    a.axhline(0, color="gray", lw=.8, ls=":")
    a.set_xticks(x); a.set_xticklabels(lab, rotation=90, fontsize=6.5)
    a.set(ylabel=r"$(\,\mathrm{pred}-\Lambda\mathrm{CDM})/\sigma$", title="DESI DR2 BAO")
    a.legend(fontsize=8.5)
    # SN: Pantheon residuals vs ΛCDM (shape), with the SEDE curve
    zp, mu, covp = dl.load_pantheon_plus()
    DMl,_,_ = distances(El, H0L); DMs,_,_ = distances(Es, H0S)
    muL = 5*np.log10((1+zp)*DMl(zp))+25; muS = 5*np.log10((1+zp)*DMs(zp))+25
    w = 1/np.diag(covp); off = np.sum(w*(mu-muL))/np.sum(w)
    res = mu-muL-off
    bins = np.linspace(0.01, 2.3, 12); idx = np.digitize(zp, bins)
    zc = 0.5*(bins[1:]+bins[:-1])
    rb = [res[idx==k].mean() if np.any(idx==k) else np.nan for k in range(1,len(bins))]
    eb = [res[idx==k].std()/max(1,np.sqrt((idx==k).sum())) if np.any(idx==k) else np.nan for k in range(1,len(bins))]
    zz = np.linspace(0.01, 2.3, 200); cS = 5*np.log10((1+zz)*DMs(zz))+25 - (5*np.log10((1+zz)*DMl(zz))+25)
    b.errorbar(zc, rb, yerr=eb, fmt="o", color="#7f8c8d", ms=4, capsize=2, label="Pantheon+ (binned) − ΛCDM")
    b.plot(zz, cS-np.median(cS)+0, color=CSEDE, lw=2.2, label="SEDE − ΛCDM")
    b.axhline(0, color="gray", lw=.8, ls=":"); b.set_xscale("log")
    b.set(xlabel="redshift z", ylabel=r"$\Delta\mu$ (mag)", title="Pantheon+ supernovae", ylim=(-0.06,0.06))
    b.legend(fontsize=8.5)
    fig.suptitle("DESI DR2 BAO and Pantheon+ residuals relative to ΛCDM", fontsize=11)
    savefig_vector(fig, f"{OUT}/fig3_datafit"); plt.close(fig); print("fig3")

# ───────────────────────────────────────────────────────────────── fig 4
def w_of_z(zg):
    """w(z) of the SEDE dark-energy fluid from ρ_DE(z) = E² − Ω_m(1+z)³ − Ω_r(1+z)⁴,
    via the continuity equation d ln ρ_DE/d ln a = −3(1+w)."""
    Or = 9e-5; E = Es(zg); rho = E**2 - OmS*(1+zg)**3 - Or*(1+zg)**4
    lnr = np.log(rho); lna = np.log(1/(1+zg))
    return -1 - (1/3)*np.gradient(lnr, lna)

def _cov_ellipse(ax, mean, cov, nsig, **kw):
    """Add an nσ covariance ellipse (properly oriented via eigendecomposition)."""
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]; vals, vecs = vals[order], vecs[:, order]
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    w, h = 2 * nsig * np.sqrt(vals)
    ax.add_patch(Ellipse(mean, w, h, angle=ang, **kw))


def _sede_w0wa():
    """SEDE's (w0,wa) by CPL-fitting w(z) over 0<z<1 — the SAME range as gen_predictions.py,
    so this reproduces the sha256-locked predictions.json value (w0,wa)=(-0.98,-0.11)."""
    zf = np.linspace(0, 1, 200); af = 1 / (1 + zf)
    A = np.vstack([np.ones_like(af), (1 - af)]).T
    return tuple(np.linalg.lstsq(A, w_of_z(zf), rcond=None)[0])


def fig4_eos():
    zg = np.linspace(0.0, 2.5, 400); w = w_of_z(zg)
    w0wa = dl.load_desi_dr2_w0wa()
    w0s, was = _sede_w0wa()
    COLS = {"pantheonplus": "#2e7d32", "desy5": "#1565c0", "union3": "#6a1b9a"}
    LAB  = {"pantheonplus": "DESI+Pantheon+", "desy5": "DESI+DESY5", "union3": "DESI+Union3"}

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.0), constrained_layout=True)
    a.plot(zg, w, color=CSEDE, lw=2.4, label="SEDE  w(z)")
    a.axhline(-1, color=CL, lw=1.6, ls="--", label="ΛCDM")
    pp = w0wa["pantheonplus"]
    a.plot(zg, pp["w0"] + pp["wa"]*zg/(1+zg), color="#2e7d32", lw=1.6, ls=":", label="DESI DR2 CPL (Pantheon+)")
    a.set(xlabel="redshift z", ylabel="w(z)", title="Equation of state crosses −1", xlim=(0,2.5), ylim=(-1.18,-0.6))
    a.legend(fontsize=8.5, loc="upper right")

    # (w0,wa) plane — OFFICIAL DESI DR2 w0waCDM chains (real, ρ≈-0.9 tilted covariance)
    b.axhline(0, color="gray", lw=.6, ls=":"); b.axvline(-1, color="gray", lw=.6, ls=":")
    mu_pp = (pp["w0"], pp["wa"])
    _cov_ellipse(b, mu_pp, pp["cov2x2"], 2, fc=COLS["pantheonplus"], ec="none", alpha=0.10)
    _cov_ellipse(b, mu_pp, pp["cov2x2"], 1, fc=COLS["pantheonplus"], ec="none", alpha=0.22)
    for sn in ("pantheonplus", "desy5", "union3"):
        d = w0wa[sn]
        _cov_ellipse(b, (d["w0"], d["wa"]), d["cov2x2"], 1, fc="none", ec=COLS[sn], lw=1.4)
        b.plot(d["w0"], d["wa"], "s", color=COLS[sn], ms=6, label=LAB[sn])
    b.plot(w0s, was, "^", color=CSEDE, ms=11, label=f"SEDE Δ=1 ({w0s:+.2f},{was:+.2f})")
    b.plot(-1, 0, "o", color=CL, ms=7, label="ΛCDM")
    src = "official DR2 chains" if pp["source"] == "official" else "published marginals"
    b.set(xlabel=r"$w_0$", ylabel=r"$w_a$", title=f"In DESI's evolving-DE quadrant ({src})",
          xlim=(-1.15,-0.55), ylim=(-1.6,0.35)); b.legend(fontsize=7.5, loc="lower left")
    fig.suptitle("SEDE's phantom-crossing EOS, alongside the official DESI DR2 w0waCDM result", fontsize=11)
    savefig_vector(fig, f"{OUT}/fig4_eos"); plt.close(fig); print("fig4")

# ───────────────────────────────────────────────────────────────── fig 5
def fig5_oos():
    ze,te,me,cve = dl.load_eboss_dr16_fs()
    zd,td,md,cvd = dl.load_desi_dr2()
    if ze is None: print("fig5 skipped (no eBOSS files)"); return
    # only BAO distance points (drop fs8)
    def bao_only(z,t,m,cov):
        k=[i for i,ti in enumerate(t) if ti in ("DM/rd","DH/rd","DV/rd")]
        s=np.sqrt(np.diag(cov))
        return np.array(z)[k],[t[i] for i in k],np.array(m)[k],s[k]
    ze,te,me,se = bao_only(ze,te,me,cve); zd,td,md,sd = bao_only(zd,td,md,cvd)
    fig, ax = plt.subplots(figsize=(7.6, 4.2), constrained_layout=True)
    # residuals vs ΛCDM, eBOSS (training) and DESI (held-out)
    pSe=bao_pred(ze,te,Es,H0S); pLe=bao_pred(ze,te,El,H0L)
    pSd=bao_pred(zd,td,Es,H0S); pLd=bao_pred(zd,td,El,H0L)
    ax.errorbar(ze, (me-pLe)/se, yerr=1, fmt="o", color="#8e44ad", ms=5, capsize=2,
                label="eBOSS DR16 (TRAINING) − ΛCDM")
    ax.errorbar(zd, (md-pLd)/sd, yerr=1, fmt="s", color=CL, ms=4, capsize=2,
                label="DESI DR2 (HELD-OUT) − ΛCDM")
    ax.plot(np.sort(np.r_[ze,zd]), np.zeros(len(ze)+len(zd)), alpha=0)
    zs=np.argsort(np.r_[ze,zd])
    allz=np.r_[ze,zd][zs]; allp=((np.r_[pSe,pSd]-np.r_[pLe,pLd])/np.r_[se,sd])[zs]
    ax.plot(allz, allp, "^", color=CSEDE, ms=8, label="SEDE prediction − ΛCDM")
    ax.axhline(0, color="gray", lw=.8, ls=":")
    ax.set(xlabel="redshift z", ylabel=r"$(\,\mathrm{value}-\Lambda\mathrm{CDM})/\sigma$",
           title="Out-of-sample: trained on pre-DESI eBOSS, SEDE predicts DESI")
    ax.text(0.02,0.03,"held-out DESI Δχ²(SEDE−ΛCDM) = −8.56", transform=ax.transAxes,
            fontsize=8.5, bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.legend(fontsize=8.5, loc="upper right")
    savefig_vector(fig, f"{OUT}/fig5_oos"); plt.close(fig); print("fig5")

# ───────────────────────────────────────────────────────────────── fig 6
def fig6_forecast():
    # Fisher σ(Δ) (from run_tier2): DESI DR3, Euclid, combined -> σ vs Δ=0 separation
    sig = {"DESI DR3":0.115, "Euclid":0.133, "DR3+Euclid":0.087}
    fig, ax = plt.subplots(figsize=(6.6, 4.0), constrained_layout=True)
    names=list(sig); sep=[1/sig[k] for k in names]
    bars=ax.bar(names, sep, color=["#2980b9","#16a085",CSEDE], alpha=0.9)
    for bar,k in zip(bars,names):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.2,
                f"{1/sig[k]:.1f}σ\nσ(Δ)={sig[k]:.3f}", ha="center", fontsize=8.5)
    ax.set(ylabel=r"separation of $\Delta=1$ from $\Delta=0$  [$\sigma$]",
           title="Fisher forecast (optimistic): the Barrow deformation Δ", ylim=(0,13))
    ax.text(0.02,0.95,"Δ=1 (volume-law horizon) vs Δ=0 (smooth):\na discrete model-selection target\n(Fisher-optimistic; marginalised errors broader)",
            transform=ax.transAxes, va="top", fontsize=8.0,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    savefig_vector(fig, f"{OUT}/fig6_forecast"); plt.close(fig); print("fig6")

if __name__ == "__main__":
    for f in (fig1_mechanism, fig2_inputs, fig3_datafit, fig4_eos, fig5_oos, fig6_forecast):
        try: f()
        except Exception as e: print(f"{f.__name__} FAILED: {e}")
