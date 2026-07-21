# Tensor sector & standard-siren predictions (WP-A / WP-B)

Deliverables integrating the FPAB tensor-sector result (Ξ₀ = 1) into the SEDE papers.
Written into the **canonical** `research/sede-dev` tree (source of the public `papers/`).

## Label crosswalk (handoff package → this repo)

The originating handoff package used a two-paper naming that predates the current split:

| Package label | This repo |
|---|---|
| "Paper I (`SEDE_cosmology.tex`)" | **Paper III** — `paper/SEDE_cosmology.md` (`papers/3-sede-cosmology`) |
| "Paper II" (driven non-equilibrium) | **Paper I** — `paper/SEDE_foundations.md` (`papers/1-sede-foundations`) |
| the θ₊ = 0 elliptic result it cites | **Paper II** — `SEDE_count.md` eq. 2.9 (`papers/2-sede-count`) |
| "§7.2 roughening (corrected)" | retracted (v2.8, 2026-07-09); now the count paper's eq. 2.9 result — P3–P5 already withdrawn in §6 |
| "§5.2 consistency item" | not in foundations; do not use that pointer |
| action `G₂ + G₃□φ` | paper convention is `K(X,φ) − G₃□φ` (reconciled in the inserted text) |

## Files

- `fpab_tensor_sector.py` — two-tier SymPy proof (first-principles O(ε²) expansion + Bellini–Sawicki cross-check). Reproduces `verification_output.txt` line-for-line. Needs `sympy ≥ 1.12`.
- `verification_output.txt` — frozen reference output; diff any rerun against it.
- `prerecomb_earlyde.py` — WP-A2. Recomputes f_sat and ρ_DE/ρ_tot at recombination from the pipeline (`sede.friedmann.E_SEDE_volume`, canonical Δ=1/λ=0.5). Run: `PYTHONPATH=<sede-dev root> python prerecomb_earlyde.py`.
- `fig_xi0_w0.py` — WP-B. The (Ξ₀, w₀) discriminant figure → `outputs/fig_xi0_w0.{pdf,png}`. Every plotted value is a named constant with a citation comment; ends with validation asserts.

## Key numbers (pipeline / script-traced)

- f_sat(z=1089.9) = **3.4×10⁻⁶**, ρ_DE/ρ_tot = **1.0×10⁻¹⁰** (note asserted ~10⁻⁶, ~10⁻¹⁰). z=3 cross-check Ω_DE = 0.0286 reproduces Fig. 5.
- SEDE w₀: CPL(0<z<1) **−0.984** (literal w(0)=−0.996, crossing z≈0.18) — matches the paper's (−0.98,−0.11).
- Tensor action `L₂ = (M_pl²a³/8)[ḣ² − (k²/a²)h²]`: α_M=α_T=0, c_T=1, zero mass, zero residual — for arbitrary K, G₃.
- RT-nonlocal minimal (rival): Ξ₀=0.93, n=2.59 (Belgacem+ 2019, arXiv:1907.02047, Table 2).

## Insertions made (additive only; no existing equation/claim altered)

- `SEDE_cosmology.md` §3.4 — c_T=1 assertion upgraded to "derived" with a pointer (no sentence altered).
- `SEDE_cosmology.md` §6 — new subsection "Tensor sector and standard-siren predictions", incl. WP-A2 pre-recombination corollary.
- `SEDE_foundations.md` §8.1 — new clarification (iv): driven-NESS tensor cross-check (four channels), the GW-silence honesty item, and the Δ_BH=0 ringdown/no-echo line.
- `refs.bib` — added BelliniSawicki:2014, Belgacem:2019lisa, Belgacem:2019nonlocal, Mukherjee:2021siren.

## Notes on the rival / forecast values (all resolved in-text)

1. **Screening at the source** — argued from the paper's own near-horizon scale-dependence (μ turns on only across k₅₀=0.70 aH … k₉₀=2.1 aH, i.e. Gpc scales), so the enhancement is absent at sub-parsec binary scales; Vainshtein screening reinforces it. No dependence on an unpinned r_V.
2. **σ(Ξ₀) forecast** — verified against Mukherjee+ 2021: Ξ₀=0.98⁺⁰·⁰⁴₋₀·²³ marginalised, 0.99±0.02 fixed-n, ~3500 events of ~30 M⊙ to z=0.5. Text quotes 0.02 (fixed) → 0.04 (marginalised).
3. **SEDE w₀ error bar** = ±0.004 in the figure, propagated from σ(Ω_m)=0.01 via the pipeline (w₀ is derived from Ω_m, not fit).
4. **f(R)/G₄(φ)** — drawn as an illustrative running-M* band (its Ξ₀ runs with the model's free-function amplitude); a Barrow-modified-gravity point is omitted rather than fabricated. RT-nonlocal is a real, cited point (Ξ₀=0.93, n=2.59).
