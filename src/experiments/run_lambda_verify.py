#!/usr/bin/env python3
"""Independently verify SEDE_V2's claim: parameter-free lambda=0.5 + gamma=theory
SEDE-H beats LCDM under the gold-standard CAMB CMB. Uses MY E_SEDE_lambda (same
H^{2 lambda}-coupling fixed-point form as their cosmology.py) in MY CAMB joint."""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve
from sede import friedmann as fr
from run_camb_joint import load, R_PL, LA_PL, CMB_CINV, SR, SLA, SHOES, OMBH2_PRIOR
C=299792.458; data=load()

def w_of_a(Om, gamma, lam, n=250, Omega_r=9.0e-5):
    # rho_DE must exclude BOTH matter AND radiation (E_SEDE_lambda carries Omega_r in E^2);
    # omitting the radiation term contaminates w(z) at high z and, since CAMB's PPF re-adds
    # radiation, double-counts it. (checkpoint-review bug fix; conservative for the headline.)
    z=np.linspace(0,8,n); E=fr.E_SEDE_lambda(z,Om,gamma,lam,Omega_r=Omega_r)
    rho=np.maximum(E**2-Om*(1+z)**3-Omega_r*(1+z)**4,1e-8); w=-1+(1/3.)*(1+z)*np.gradient(np.log(rho),z)
    a=1/(1+z); i=np.argsort(a); return a[i],w[i]

def chi2(Om,H0,ombh2,MB,s8,gamma,lam,use_shoes=True):
    h=H0/100.; pa=camb.CAMBparams(); pa.set_cosmology(H0=H0,ombh2=ombh2,omch2=Om*h**2-ombh2,mnu=0.06)
    a,w=w_of_a(Om,gamma,lam); de=dark_energy.DarkEnergyPPF(); de.set_w_a_table(a,w); pa.DarkEnergy=de
    pa.WantCls=False; pa.WantTransfer=False; bg=camb.get_background(pa)
    d=bg.get_derived_params(); rd=d['rdrag']; zs=d['zstar']; rs=d['rstar']; c=0.0
    z,t,m,icov=data['desi']
    pred=np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
        else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz,tp in zip(z,t)])
    dd=m-pred; cb=float(dd@icov@dd); c+=cb
    zp,mu,chol=data['pan']; dmu=mu-(5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB); c+=float(dmu@cho_solve(chol,dmu))
    DMz=bg.comoving_radial_distance(zs); R=np.sqrt(Om)*H0*DMz/C; lA=np.pi*DMz/rs
    v=np.array([R-R_PL,lA-LA_PL]); c+=float(v@CMB_CINV@v)
    c+=((ombh2-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2
    if use_shoes: c+=((H0-SHOES[0])/SHOES[1])**2
    zc,H,icc=data['cc']; dH=H-bg.hubble_parameter(zc); c+=float(dH@icc@dH)
    zf,fo,fe=data['fs8']; Dd,fd=fr.compute_growth_model(zf,Om,lambda zz:bg.hubble_parameter(np.atleast_1d(zz))/H0)
    c+=float(np.sum(((fo-fd*s8*Dd)/fe)**2)); return c,R,lA,cb

def fit(gamma_fixed,lam,gfree=False,use_shoes=True):
    def obj(v):
        if gfree: Om,H0,ob,MB,s8,g=v
        else: Om,H0,ob,MB,s8=v; g=gamma_fixed
        if not(0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and -20.5<MB<-18.5 and 0.62<s8<0.90 and 0.3<g<6): return 1e9
        try: return chi2(Om,H0,ob,MB,s8,g,lam,use_shoes)[0]
        except Exception: return 1e9
    best=None
    seeds=[[0.305,68.4,0.02237,-19.40,0.78]+([1.5] if gfree else [])]
    for s in seeds:
        r=minimize(obj,s,method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-3,maxiter=6000))
        if best is None or r.fun<best.fun: best=r
    return best

if __name__=="__main__":
    import sys
    # LCDM baselines (full + no-SHOES) from run_camb_joint marginalised best fits
    LC_full, LC_core = 1440.79, None
    # compute LCDM core (no SHOES) once
    def lcdm(use_shoes):
        def obj(v):
            Om,H0,ob,MB,s8=v
            if not(0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and -20.5<MB<-18.5 and 0.62<s8<0.90): return 1e9
            h=H0/100.; pa=camb.CAMBparams(); pa.set_cosmology(H0=H0,ombh2=ob,omch2=Om*h**2-ob,mnu=0.06); pa.WantCls=False; pa.WantTransfer=False
            bg=camb.get_background(pa); d=bg.get_derived_params(); rd=d['rdrag']; zs=d['zstar']; rs=d['rstar']; c=0.0
            z,t,m,icov=data['desi']
            pred=np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
                else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz,tp in zip(z,t)])
            dd=m-pred; c+=float(dd@icov@dd)
            zp,mu,chol=data['pan']; dmu=mu-(5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB); c+=float(dmu@cho_solve(chol,dmu))
            DMz=bg.comoving_radial_distance(zs); R=np.sqrt(Om)*H0*DMz/C; lA=np.pi*DMz/rs; v=np.array([R-R_PL,lA-LA_PL]); c+=float(v@CMB_CINV@v)
            c+=((ob-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2
            if use_shoes: c+=((H0-SHOES[0])/SHOES[1])**2
            zc,H,icc=data['cc']; dH=H-bg.hubble_parameter(zc); c+=float(dH@icc@dH)
            zf,fo,fe=data['fs8']; Dd,fd=fr.compute_growth_model(zf,Om,lambda zz:bg.hubble_parameter(np.atleast_1d(zz))/H0); c+=float(np.sum(((fo-fd*s8*Dd)/fe)**2))
            return c
        return minimize(obj,[0.30,68.5,0.02237,-19.40,0.78],method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-3,maxiter=6000)).fun
    LCf=lcdm(True); LCc=lcdm(False)
    print("LCDM: full=%.2f  core(no SH0ES)=%.2f"%(LCf,LCc))
    for lam in [0.5, 1.0]:
        # gamma fixed to theory 1.5, parameter-free
        bf=fit(1.4964,lam,gfree=False,use_shoes=True); bc=fit(1.4964,lam,gfree=False,use_shoes=False)
        # gamma free
        bg_=fit(None,lam,gfree=True,use_shoes=True); gval=bg_.x[5]
        print("λ=%.2f: γ=theory PARAM-FREE  Δχ²(full)=%+.2f  Δχ²(core)=%+.2f  | γ-free best γ=%.2f Δχ²(full)=%+.2f"
              %(lam,bf.fun-LCf,bc.fun-LCc,gval,bg_.fun-LCf), flush=True)
