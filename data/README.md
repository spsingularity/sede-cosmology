# Data vectors

Small, standard public data vectors used by `reproduce.py` are vendored here:

- `desi_dr2_desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt` — DESI DR2 BAO (D_M/r_d, D_H/r_d)
- `Pantheon+SH0ES.dat` — Pantheon+ distance moduli
- `sdss_DR16_BAOplus_{LRG,QSO}_FSBAO_DMDHfs8*.{dat,txt}` — eBOSS DR16 full-shape (out-of-sample)
- `desi_dr2_w0wa_{pantheonplus,desy5,union3}.{covmat,margestats}` — DESI DR2 (w0,wa) chains

## Not vendored (fetched separately — large or external)

- **`Pantheon+SH0ES_STAT+SYS.cov` (~32 MB)** — the full Pantheon+ stat+sys covariance. Without
  it, `figures.py` falls back to the stat-only errors for `fig3_datafit` (a WARN is printed).
  Download from the Pantheon+ release (Brout et al. 2022) and place it in this directory for the
  full covariance.
- **Moresco (2022) cosmic-chronometer covariance** — cache as `data/moresco_cov.npy` if needed.

The full reproduction pipeline (marginalised MCMC, the false-preference mock calibration behind
`fig_F5_calibration`, and the CAMB-in-the-loop numbers) lives in the code repository named in the
paper's Data Availability statement; this manuscript repo carries the figure-level reproduction.
