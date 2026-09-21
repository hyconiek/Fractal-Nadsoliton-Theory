# MP7-035 — local finite-size rounding and equal-cap shift

Scientific state: **PROVED_GLOBAL_LEADING_ASYMPTOTIC; PROVED_INTERVAL_ASSISTED_LOCAL_FINITE_N_ERROR; GLOBAL_QUANTITATIVE_REMAINDER_OPEN**.

This result is conditional on the finite-N statistical extension of MP7-031/032.
The leading equal-mass shift is global after MP7-017: the uniform state and the
D12 orbit of twelve localized minima exhaust the global minimizers at the
coexistence point.  The explicit finite-N error constants derived below are
currently for fixed local caps at the certified equal-energy gain only; a
uniform quantitative error for the full measure in `g=g_eq+c/N` remains open.

## 1. Local Laplace asymptotic

For a fixed sufficiently small cap around a nondegenerate stable minimum,
Taylor expansion gives

`Phi(theta_*+x)=Phi_* + 1/2 x^T H x + O(||x||^3)`.

The explicit positive Hessian and cap-boundary gap from MP7-034 permit the
standard scaling `x=u/sqrt(N)`.  On every fixed bounded u-set the cubic
remainder tends uniformly to zero, while positive quadratic growth supplies a
dominating Gaussian.  The cap boundary contributes exponentially less than the
core.  Therefore, for each named cap,

`I_N = exp(-N Phi_*) (2 pi/N)^(7/2) det(H)^(-1/2) [1+o(1)]`.

Because the compared minima are at the same g, `G7=g H7` and the common g
factor cancels from their determinant ratio.  D12 invariance makes all twelve
localized prefactors identical.

The initial Krawczyk root boxes were far too small to give useful explicit
finite-N constants.  The later interval Hessian/T3 calculation below enlarges
the certified convex neighborhoods and supplies a quantitative local error.

## 2. Ratio near the local equal-energy event

Let

`Delta V(g) = V_localized(g)-V_uniform(g)`.

At the certified local event `g_eq`, `Delta V=0`, and R7P-028 gives

`Delta V'(g_eq) in
 [-0.452413732176580, -0.452413730675052]`.

MP7-034 gives

`A = log 12 + 0.5 log(det G_uniform/det G_localized)
   in [-0.609295503686632, -0.609295351426053]`.

Hence the local-family cap ratio has the asymptotic form

`log R_N(g) = A - N Delta V(g) + o(1)`.

For the finite-size window

`g = g_eq + c/N`,

smooth continuation of the nondegenerate localized branch yields

`log R_N -> A - c Delta V'(g_eq)`.

Thus the two-family conditional share approaches the logistic expression

`P_loc,conditional -> R/(1+R),  R=exp(A-c Delta V')`.

## 3. Equal-cap-mass displacement

The leading equal-cap condition `R=1` gives

`c_* = A / Delta V'(g_eq)`.

Outward interval division yields

`c_* in [1.34676582095488, 1.34676616197633]`.

Therefore

`g_N = g_eq + c_*/N + o(1/N)`

for equality of the named localized-family and uniform cap masses.

The sign is **positive**.  This is the opposite of what one would infer from
multiplicity alone: the localized family must move slightly to the higher-g
side, where its energy is lower, to compensate for the much larger Gaussian
width of the uniform minimum.

## 4. Explicit local finite-N error at exact coexistence

This subsection quantifies the two local Laplace integrals at the certified
`g_eq`.  It does not yet quantify the global complement or a moving-minimum
`tube g=g_eq+c/N`.

Let `d=7`, let `H0` be the Hessian at the relevant minimum, and take a Euclidean
cap of radius `r` around that minimum.  A direct interval evaluation of the
third derivative tensor over a coordinate cube containing that ball gives an
upper Frobenius bound `B3`; hence it is also an operator-norm bound.  Combining
this with the certified central Hessian gives

`lambda_min(H(theta)) >= m := lambda_min(H0)-B3*r`

throughout the ball.  The optimized calculation uses `r=0.005` for both wells:

- localized: `B3 <= 0.217933205897420`, `m >= 0.104735758532760`;
- uniform: `B3 <= 0.376250164153595`, `m >= 0.0718737719957623`.

Choose the Gaussian core radius

`rho_N = c sqrt(log N/N) <= r`.

Taylor's theorem on that core gives the exponent error

`delta_N <= B3 c^3 (log N)^(3/2)/(6 sqrt(N))`.

For the reference Gaussian with quadratic form `H0`, the union bound and
`H0 >= m I` give

`q_G <= 2d exp[-N m rho_N^2/(2d)]`.

For the actual cap integral, strong convexity bounds the part outside the core
relative to the full Gaussian prefactor by

`q_A <= [sqrt(det H0)/m^(d/2)] 2d exp[-N m rho_N^2/(2d)]`.

Therefore, if `L_N` denotes the local cap integral divided by its leading
Gaussian Laplace term,

`exp(-delta_N)(1-q_G) <= L_N <= exp(delta_N)+q_A`.

The optimization over the tail exponent/c-scale is deterministic.  Sufficient
common thresholds for **both** the uniform and localized cap integrals are:

| requested relative error | sufficient N |
|---:|---:|
| 25% | `N >= 6.084117415e8` |
| 10% | `N >= 4.413892270e9` |
| 5% | `N >= 2.058778730e10` |
| 1% | `N >= 7.552785569e11` |

These are deliberately conservative theorem bounds, not estimates of the N at
which the asymptotic becomes numerically accurate in practice.  They use only
interval derivative/curvature bounds and Gaussian tail inequalities; no Monte
Carlo or fitted asymptotic is used.

The executable records are
`results/MP7-035_cap_lipschitz_scan.json` and
`results/MP7-035_optimized_explicit_error.json`.

## 5. Upgrade after MP7-017: global finite-N asymptotic

MP7-017 changes the interpretation of the preceding local calculation.  At
`g_eq` the full seven-dimensional dual has exactly thirteen global minima: the
uniform state and twelve localized D12 images.  All are nondegenerate stable
minima.  The potential is coercive.

Choose disjoint fixed caps around these thirteen minima.  On the closed
complement, coercivity permits restriction to a sufficiently large compact
ball, and continuity plus the fact that no other point is globally minimizing
gives a strictly positive energy gap `delta>0`.  The outer tail is also
separated by coercivity.  Therefore the full auxiliary-field partition is, to
leading order, the sum of these thirteen Laplace contributions.

The same statement is uniform for `g=g_eq+c/N` with c in any fixed bounded set
once N is large: the nondegenerate minima continue smoothly, while the compact
complement gap persists by continuity.  Consequently

`R_N(c) = Z_localized-family / Z_uniform
        -> exp(A-c DeltaV'(g_eq))`

is now the **global** equilibrium weight ratio of the two competing phase
families in the finite-N extension, not merely a conditional ratio of selected
caps.  The resulting global limiting localized fraction is

`R/(1+R)`.

The equal-mass displacement remains

`g_N = g_eq + c_*/N + o(1/N)`,

`c_* in [1.34676582095488,1.34676616197633]`.

The leading global asymptotic statement is therefore proved.  The local
Laplace factors now also have explicit fixed-`g_eq` error/N0 bounds from
Section 4.  What remains open is a **global quantitative** remainder: an
explicit complement-gap contribution and derivative/continuation constants
uniform for the moving minima in `g=g_eq+c/N`.
