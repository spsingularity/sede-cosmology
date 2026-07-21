"""
FPAB (cubic kinetic-gravity-braiding) tensor sector:  alpha_M, alpha_T, Xi_0
===========================================================================
Part B: first-principles second-order tensor expansion of
        L = sqrt(-g) [ (Mpl^2/2) R + G2(phi, X) + G3(phi, X) Box(phi) ]
        around flat FRW, using the exponential TT parametrization
        g_ij = a^2 (exp[eps * h_ij])_ij ,  h_xx = -h_yy = h(t) cos(k z).
        In this parametrization sqrt(-g) = a^3, X = phidot^2/2 and
        Box(phi) = -(phiddot + 3 H phidot) are EXACTLY eps-independent,
        so G2 and G3 contribute nothing at O(eps^2): the entire quadratic
        tensor action comes from the Einstein-Hilbert term.
Part A: cross-check with the general Horndeski tensor-sector formulas
        (Bellini-Sawicki alpha-parametrization).
"""
import sympy as sp

# ---------- Part B: direct second-order expansion ----------
t, x, y, z = sp.symbols('t x y z')
k, eps, Mpl = sp.symbols('k epsilon M_pl', positive=True)
a   = sp.Function('a', positive=True)(t)
h   = sp.Function('h')(t)
phi = sp.Function('phi')(t)

psi = eps * h * sp.cos(k * z)                      # single + polarization, wave along z
g   = sp.diag(-1, a**2 * sp.exp(psi), a**2 * sp.exp(-psi), a**2)
coords = [t, x, y, z]
ginv = g.inv()
detg = sp.simplify(g.det())
print("det g  =", detg, "   (exactly -a^6: sqrt(-g)=a^3, eps-independent)")

# Christoffel symbols
def christoffel(g, ginv, coords):
    n = len(coords)
    Gamma = [[[0]*n for _ in range(n)] for _ in range(n)]
    for mu in range(n):
        for nu in range(n):
            for rho in range(n):
                s = 0
                for sig in range(n):
                    s += ginv[mu, sig] * (sp.diff(g[sig, nu], coords[rho])
                                          + sp.diff(g[sig, rho], coords[nu])
                                          - sp.diff(g[nu, rho], coords[sig]))
                Gamma[mu][nu][rho] = sp.simplify(s / 2)
    return Gamma

Gamma = christoffel(g, ginv, coords)

# Ricci tensor and scalar
def ricci_scalar(g, ginv, Gamma, coords):
    n = len(coords)
    R = 0
    for nu in range(n):
        for rho in range(n):
            Rnr = 0
            for mu in range(n):
                Rnr += sp.diff(Gamma[mu][nu][rho], coords[mu])
                Rnr -= sp.diff(Gamma[mu][nu][mu], coords[rho])
                for sig in range(n):
                    Rnr += Gamma[mu][mu][sig] * Gamma[sig][nu][rho]
                    Rnr -= Gamma[mu][rho][sig] * Gamma[sig][nu][mu]
            R += ginv[nu, rho] * Rnr
    return R

R = ricci_scalar(g, ginv, Gamma, coords)

# Verify the matter-sector scalars are eps-independent (exact, no expansion)
X_exact   = sp.Rational(1, 2) * phi.diff(t)**2          # X = -g^{00} phidot^2 /2, g^{00}=-1
sqrtg     = a**3
box_phi   = sp.simplify(sum(sp.diff(sqrtg * ginv[mu, 0] * phi.diff(t), coords[mu])
                            for mu in range(4)) / sqrtg)
print("Box(phi) =", box_phi, "   (eps-independent -> G2, G3 give ZERO at O(eps^2))")

# Einstein-Hilbert Lagrangian density, expanded to O(eps^2)
L_EH  = sp.Rational(1, 2) * Mpl**2 * sqrtg * R
L2    = sp.series(L_EH, eps, 0, 3).removeO().expand()
L2_e2 = L2.coeff(eps, 2)

# average over one wavelength in z  (cos^2 -> 1/2, sin^2 -> 1/2, cross terms -> 0)
period = 2 * sp.pi / k
L2_avg = sp.integrate(L2_e2, (z, 0, period)) / period
L2_avg = sp.simplify(sp.expand(L2_avg))
print("\nRaw O(eps^2) Lagrangian (z-averaged):")
sp.pprint(L2_avg)

# extract canonical coefficients (fully expanded polynomial in h, hdot, hddot)
hd, hdd = h.diff(t), h.diff(t, 2)
Lc = sp.expand(L2_avg)
K      = sp.simplify(Lc.coeff(hd, 2).coeff(h, 0))      # coeff of hdot^2
c_h2   = sp.simplify(Lc.coeff(h, 2).coeff(hd, 0))      # coeff of h^2
c_hddh = sp.simplify(Lc.coeff(hdd, 1))                 # coeff of hddot*h  (should be 0)
c_hdh  = sp.simplify(Lc.coeff(hd, 1).coeff(h, 1))      # coeff of hdot*h   (should be 0)
Gr     = sp.simplify(-c_h2 * a**2 / k**2)              # gradient coeff:  -Gr k^2/a^2 h^2
mass0  = sp.simplify(c_h2.subs(k, 0))                  # k-independent mass term (must vanish)
resid  = sp.simplify(Lc - (K * hd**2 + c_h2 * h**2))   # anything else left over

print("\nCanonical quadratic action  L2 = K hdot^2 - Gr (k^2/a^2) h^2")
print("K  (kinetic)        =", K)
print("Gr (gradient)       =", Gr)
print("c_T^2 = Gr/K        =", sp.simplify(Gr / K))
print("h''h and h'h coeffs =", c_hddh, ",", c_hdh, "  (no IBP even needed)")
print("mass term at k=0    =", mass0, "  (no graviton mass, no background-EOM input required)")
print("residual terms      =", resid)

Mstar2 = sp.simplify(8 * K / a**3)   # matches (Mstar^2/8) a^3 hdot_ij hdot_ij with <h_ij h_ij> = h^2 here
print("M*^2 inferred = 8K/a^3 =", Mstar2, " -> constant in a  ->  alpha_M = 0")

# ---------- Part A: Horndeski alpha-parametrization cross-check ----------
print("\n--- Part A: general Horndeski tensor formulas (Bellini-Sawicki) ---")
ph, Xs, H = sp.symbols('phi X H')
pd = sp.symbols('phidot')
G4 = sp.Function('G4')(ph, Xs); G5 = sp.Function('G5')(ph, Xs)
Mstar2_gen = 2 * (G4 - 2 * Xs * G4.diff(Xs) + Xs * G5.diff(ph) - pd * H * Xs * G5.diff(Xs))
alphaT_gen = 2 * Xs * (2 * G4.diff(Xs) - 2 * G5.diff(ph)
                       - (sp.Symbol('phiddot') - pd * H) * G5.diff(Xs)) / Mstar2_gen
subs_FPAB = {G4: Mpl**2 / 2, G5: 0}
Mstar2_FPAB = Mstar2_gen.subs(subs_FPAB).doit().simplify()
alphaT_FPAB = alphaT_gen.subs(subs_FPAB).doit().simplify()
print("FPAB:  G2, G3 arbitrary;  G4 = Mpl^2/2;  G5 = 0")
print("M*^2   =", Mstar2_FPAB, "  (constant -> alpha_M = 0)")
print("alpha_T =", alphaT_FPAB, "  (c_T = 1)")
print("\ndelta(z) = -alpha_M/2 = 0  for all z   =>   d_L^GW(z) = d_L^EM(z)   =>   Xi_0 = 1 exactly, n unconstrained.")
