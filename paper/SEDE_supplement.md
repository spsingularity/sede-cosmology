# Supplementary Material for "Structural Entropy Dark Energy: a fixed-parameter, growth-gated holographic dark-energy model without a cosmological constant"

## Abstract

This supplement collects the supporting motivations and numerical probes for the volume-law horizon
postulate, demoted from the main article to keep its main line focused: three further consistency
motivations for the volume endpoint, the Majorana–SYK diagonalisation that fixes the horizon state,
numerical probes of the state and counting halves of the postulate, the degrees-of-freedom accounting
that locates the postulate as a counting (not thermalisation) claim — with the SYK and causal-set
comparisons — and the summary of how far the postulate is derivable. None is a derivation; each is a
consistency check. Cross-references of the form §N refer to sections of the main article; §S-N to sections of this supplement.

## 1. Supporting motivations and numerical probes for the volume-law postulate

This supplement collects the material demoted from §8 of the main article to keep the main line
focused: three further consistency *motivations* for the volume endpoint (§8.1) and the *numerical
probes* of the state and counting halves of the postulate (§S2). None is a derivation; each is a
consistency check that makes the postulate less *ad hoc* without selecting it. (The CHR near-critical
layer, the cross-horizon black-hole reading, and the Verlinde dark-matter connection are developed in a
separate foundations companion paper [@Pandev:2026foundations].)

### 1.1 Three further consistency motivations for the volume endpoint
*(supplementing the geometric ceiling, the CKN scale/form split, and the direct area-law disfavouring
of §8.1)*

- **A driven steady state.** The horizon need not self-thermalise — its scrambling time
 (∼ ln S / H ≈ 280/H) far exceeds its age. Instead the structure-growth gate (the rise of f_sat)
 *drives* it off the ground state; with a scrambling-limited relaxation rate the volume-law state is
 *maintained* as a non-equilibrium steady state throughout structure formation (the long scrambling
 time becomes an asset). The volume law is then a fixed point of SEDE's own dynamics, not an external
 input.
- *A roughening universality class* (**withdrawn**). This motivation modelled the horizon as a
 stochastic growing surface, ∂h/∂t = ν∇²h + (λ_K/2)(∇h)² + η, and argued the deformation was fixed by
 the roughening universality class. It is retracted: the count companion paper [@Pandev:2026count] shows the
 classical apparent horizon is constraint-slaved — θ₊ = 0 is elliptic, with no growth equation at any
 order (its eq. 2.9) — and has no minimum-cut surface, so it does not roughen and there is no height
 field or roughening universality class. The ceiling Δ ≤ 1 and the value Δ = 1 rest instead on the count paper's
 leg-budget capacity theorem and on birth-selection (foundations §5), not on any surface-growth picture.
- An independent emergent-gravity proposal [@Verlinde:2016toy]. In Verlinde's programme the de Sitter
 vacuum carries a volume-law contribution to its entanglement entropy, S ∝ V_H, in addition to the
 horizon area law — the same super-area scaling SEDE adopts at Δ = 1. A second, independently-motivated
 route arriving at a volume-law entropy makes the postulate less *ad hoc*; it remains a motivation, not
 a proof, since Verlinde's volume-law derivation is itself heuristic.

### 1.2 The state is maximally entangled (Majorana-SYK diagonalisation)
The *state* half of the postulate (§S2) — that the horizon dof are maximally entangled (thermal),
not in a low-entanglement ground state — is settled and not peculiar to SEDE. Exact diagonalisation of
the Majorana Sachdev–Ye–Kitaev model [@Sachdev:1992fk; @Kitaev:2015], the canonical chaotic,
black-hole-dual system saturating the Maldacena–Shenker–Stanford bound λ_L ≤ 2πT [@Maldacena:2015waa]
(V60), confirms: the gap ratio reproduces Wigner–Dyson statistics with the
correct N-mod-8 symmetry class [@GarciaGarcia:2016mno; @Cotler:2016fpe; @Atas:2012np],
eigenstates saturate the Page value [@Page:1993df; @Vidmar:2017pak], and the OTOC scrambles completely
— while an integrable control does none of these. But this **does not fix Δ**: a black hole is a
*maximal* scrambler yet area-law (S = A/4, reproduced by string microstate counts), and SYK, being
all-to-all (geometry-free), cannot even pose the counting question.

*(The SYK scrambling figure — the gap ratio climbing through the N-mod-8 RMT classes, the Page-value
saturation of the eigenstate entanglement, and the complete OTOC scrambling of the q=4 model against
the q=2 integrable control — is omitted from this version; the underlying diagnostics are the
verification-suite check V60 and regenerate from the accompanying code.)*

### 1.3 The counting question — numerical probes (allowed, realisable, not selected)
The *counting* half (§S2) — does the horizon dof number grow as area (N ∝ R^{d-2}, Δ=0) or bulk
volume (N ∝ R^{d-1}, Δ=1)? — is where SEDE's postulate lives. Three numerical probes
(V62, V61, V63):

| question | probe | result |
|---|---|---|
| Allowed? by the established entropy bounds | volume entropy vs the smooth-area bounds | S ∝ R³ overshoots the holographic/Bousso/Bekenstein bound A/4 by ∝ R (∼10⁶¹ at the cosmic horizon); consistent only as a genuine d_H=3 fractal horizon (true area ∝ R³ ⟹ S = A_fractal/4 *exactly* saturates the bound) |
| Selected? by quantum gravity | causal-set sprinkle, Minkowski and de Sitter | both scalings present; the canonical horizon object — causal *links* — is area (R^{d-2}, the Bekenstein–Hawking recovery), in flat space and, because dS is conformally flat, identically in de Sitter. **Volume is not the default.** (measured: bulk count 1.92≈2 ⟹ Δ=1; link count 0.98≈1 ⟹ Δ=0, in 2+1D) |
| Realisable? in a holographic code | Ryu–Takayanagi min-cut, local vs nonlocal | local (lattice) connectivity ⟹ area law (exponent 1.0); nonlocal (expander) connectivity ⟹ volume law (exponent 1.9). Δ=1 is realisable; its required ingredient is *nonlocal* horizon connectivity (the all-to-all structure of SYK). |

The pattern: volume-counting is allowed (not by evading the bounds but by the horizon being
genuinely space-filling — the postulate itself), realisable (a nonlocal holographic code produces
S ∝ R^{d-1} exactly), but not selected (the frameworks we can compute in default to area for the
canonical horizon object, and a black hole is area-law). No consistency argument available to us
*derives* the volume count; equally, none *excludes* it. The two deciders are stated in §S2: the dS
static-patch Hilbert-space dimension (theoretical) and the deformation Δ (empirical, §5.6/§6).

## 2. The postulate is a dof-*counting* claim, not a thermalisation claim

The phrase "volume-law" conflates two claims of very different status. **(i) The state:** are the
horizon degrees of freedom maximally entangled (thermal) or in a low-entanglement ground state?
**(ii) The counting:** does the *number* of horizon dof grow as the area, N ∝ R^{d-2}, or the spatial
volume, N ∝ R^{d-1}? Since S ∼ (state factor) × N, the deformation Δ — through S ∝ A^{1+Δ/2} ∝
R^{(d-2)(1+Δ/2)} — is fixed by the counting (area→Δ=0, volume→Δ=1 in d=4); the state only sets
whether the prefactor is maximal. SEDE's postulate is, precisely, the counting claim.

*The state half is settled, and is not peculiar to SEDE.* Maximal scrambling has a precise meaning in
the Sachdev–Ye–Kitaev model [@Sachdev:1992fk; @Kitaev:2015], the canonical chaotic, black-hole-dual
system. But maximal scrambling **does not fix Δ**: a black hole is a *maximal* scrambler yet
area-law (S = A/4, Δ=0, reproduced by string microstate counts). Scrambling is therefore shared
with the Δ=0 black hole; and SYK, being all-to-all (geometry-free), cannot even pose the counting
question. (An explicit Majorana-SYK diagonalisation confirming the state is maximally entangled —
Wigner–Dyson statistics, Page-value saturation, complete OTOC scrambling — is given in §S1.2;
it settles the *state*, not the counting.)

*The counting half is where SEDE lives, and area is the default.* In causal-set quantum gravity — the
setting where the counting question is natural — both scalings appear, but the *canonical*
horizon-entropy object, the causal-link count that recovers the Bekenstein–Hawking law (Sorkin and
collaborators), is the area one (R^{d-2}), in flat space and identically in de Sitter. So even in
the framework most hospitable to volume-counting, area-law is the default, and SEDE's Δ=1 is the
specific, non-default identification of horizon entropy with the bulk count. (The causal-set
sprinkle measurements, the entropy-bound check, and a Ryu–Takayanagi holographic-code realisation —
volume-law requires *nonlocal* horizon connectivity — are in §S1.3; their verdict is that the
volume count is **allowed and realisable but not selected**: no consistency argument available to us
derives it, and none excludes it.)

The sharpened postulate, properly located. The one open item is therefore *not* about
thermalisation or chaos (settled, and shared with area-law black holes) but is purely a
dof-counting statement:

> the cosmic horizon's entropy counts its bulk (volume-scaling, N ∝ R³ ⟺ d_H = 3) degrees of
> freedom rather than its boundary (area-scaling, links) ones — a count that exceeds the
> holographic area bound by R/ℓ_P (the seam of §8.3).

It reduces to two deciders. *Theoretically* it is precisely *"is the Hilbert-space dimension of the de
Sitter static patch e^{Area/4} or e^{Volume}?"* — a frontier problem a complete dS holography would
settle. *Empirically* it is the deformation Δ: Δ=0 (area) is already disfavoured by present data
(§5.6), and DESI DR3 + Euclid will measure Δ to ∼0.09 (§6). SEDE's postulate is thus not an arbitrary
assumption but a *consistent, realisable bet on a well-posed open question of quantum gravity*, with a
decisive test imminent. Tellingly, the one horizon we can check today — the black hole — is area-law; on
the evidence here this asymmetry is a genuine conceptual cost. The count companion paper [@Pandev:2026count] *offers* a resolution
— a state-dependent count in which a black hole is area-law because *undriven* and only the
structure-driven cosmic horizon activates the volume law — which, if its (conjectural) mechanism holds,
would turn the cost into a predicted feature, with the black-hole side derivable and a sharp falsifier (an
isolated volume-law black hole).

## 3. Is the postulate derivable? Summary of the reduction (full treatment in the foundations companion paper)

The volume-law postulate is not monolithic: it splits into a *state* (the horizon dof are maximally
entangled/thermal), a *form* (the entropy scales as the volume), a *scale* (the magnitude is ρ_crit, not
ρ_Planck), and a *count* (the number of dof grows as the bulk, N ∝ R³, not the boundary, N ∝ R²). The
deformation Δ is fixed by the count alone. A systematic study of how derivable each piece is — five
routes (de Sitter holography/SYK, entanglement first law in a thermal state, Verlinde emergent gravity,
gravitational non-additivity, and a Bekenstein no-go), the reduction of the count to a *driven
non-equilibrium steady state*, and three closure routes — is developed in the foundations companion paper [@Pandev:2026foundations]; we state the
conclusions here, as they bear on how *ad hoc* the postulate is.

The *state* is first-principles (a maximal scrambler is volume-law-entangled; §S1.2, §S2); the
*form* reduces to thermalisation (volume-law entanglement is the generic behaviour of a thermal reference
state for regions larger than the thermal length); the *scale* is CKN (§8.2). What remains irreducible is
purely the count, and a Bekenstein no-go makes precise *why*: maximum-entropy/energy bounds give the
*area* law (§8.3), so the volume count cannot be obtained from any energy argument — it is a genuine
state/counting statement. The count is, in turn, downstream of *connectivity* (locally-connected horizon
dof give area-law entropy, nonlocal give volume-law) and connectivity of *interaction range*: gravity is
strongly long-range (1/r with α = 1 ≤ d), hence non-additive — the volume class is the *natural* one for a
self-gravitating horizon, with area-law arising only when the holographic bound intervenes for an
*isolated, equilibrium* horizon (a black hole). The cosmic horizon differs in being continuously *driven*
by structure formation, which converts the black-hole/cosmic asymmetry (§7) into a discriminator
(equilibrium → area, driven → volume) and reduces the residue to a single, falsifiable **driven-NESS
statement**: a strongly-long-range horizon, driven by structure formation, settles in the volume-law
steady state rather than relaxing to the area-law equilibrium. The foundations companion paper supplies
the two mechanisms this requires — *selection at the horizon's birth* (the N-dependent cooperative
coupling J/J_c = 2πGm²N carries the horizon, grown from Planckian size, through the ordering
bifurcation, with the Barrow tilt selecting the volume branch and long-range hysteresis locking it; no
barrier is ever crossed) and the structure *drive*, whose envelope is the gate f_sat, which *populates*
the selected volume capacity rather than selecting the branch — so the static *counting* postulate is
replaced by the weaker, more physical statement of a birth-selected branch plus derived occupancy
kinetics.

Closure — deriving even that bistable free energy from first principles — would reduce the dark sector
to Ω_m plus established physics. It takes three precise forms: a de Sitter-holography computation of the
bulk count (open, and provably not shortcut-able through energy bounds); exhibiting the volume branch as a
real saddle/state (the bistable landscape follows from gravitational cooperativity given that branch); or
*measuring* Δ. The last is decisive and already in hand: DESI DR3 + Euclid pin Δ to ∼0.09 (§6), so Δ = 1
establishes the volume count as fact and Δ = 0 refutes it. The postulate is thus a single, sharply-posed,
falsifiable statement with an imminent test — not an unconstrained assumption.

## Code and data availability

All code and analysis pipelines are publicly available at <https://github.com/spsingularity/sede-cosmology>; the
numerical diagnostics above are regenerated by the accompanying code. A tagged release is archived at Zenodo, DOI
[10.5281/zenodo.21050314](https://doi.org/10.5281/zenodo.21050314).

## References

::: {#refs}
:::
