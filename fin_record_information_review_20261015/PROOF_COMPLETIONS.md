# RI proof audit: bounded prediction, filter and likelihood scopes

Date: 2026-10-15. Analytical model arguments, not mechanized or peer reviewed.
Inherited exact spectral and OBS/MT/UF premises are bound in verification.json.

## 1. Finite-resolution prediction genuinely follows

The admitted joint weak path law has a continuous limit. Its two observation
increments, conditional on the entire signal path, have independent Gaussian
errors with variance 1/2. Thus the quantizer thresholds have zero mass. Z1
has the stationary continuous quartic density, so its ZERO-target probability
is zero. The continuous mapping theorem pays each of the 27 joint cells.

On this fixed finite alphabet the cell convergence is joint TV convergence.
The risk functional `inf_delta E[1{delta(R)!=C}]` is TV-Lipschitz, since every
decision loss lies in [0,1]. Consequently the optimal risk for this declared
finite record converges. This is not full-path TV, posterior convergence or
statistical equivalence of all microscopic histories.

## 2. Complete the qualitative strict-advantage proof

Let A=1{Z1>0}, H=1{V1>1/2}, and condition on the initial Z0=z. The cubic
additive-noise scalar diffusion is order preserving in its initial state.
One must not confuse monotonicity in Brownian **increments** with pointwise
ordering of cumulative Brownian paths.

Here is a direct route to the source's association claim. Replace the drift
temporarily by a globally Lipschitz smooth truncation agreeing on [-L,L].
In an Euler scheme with step times its Lipschitz constant <=1, each update
`x -> x+dt*b_L(x)` is nondecreasing. Terminal state and the accumulated
signal integral are then nondecreasing functions of all preceding independent
Gaussian increments. The independent observation error is one further
coordinate. For independent coordinates, increasing bounded functions have
nonnegative covariance: in one coordinate use

    Cov(f(X),g(X))=E[(f(X)-f(X'))(g(X)-g(X'))]/2>=0,

and extend by conditioning/induction over coordinates. Therefore the
conditional covariance Cov(A,H|Z0=z) is nonnegative on these schemes.
Pass first to the Lipschitz diffusion, then L->infinity by nonexplosion.
The event boundaries have zero probability, so bounded indicators pass
through these limits. No unproved finite-N information theorem is used.

The functions a(z)=P_z(Z1>0) and h(z)=P_z(V1>1/2) are strictly increasing.
For h, a common-noise ordered pair of initial states stays strictly ordered,
so its integral is strictly ordered; independent Gaussian observation noise
gives a strictly increasing tail probability. For a, the uniformly elliptic
scalar transition is strictly stochastically ordered. Explicitly, the
derivative of `P_t 1_(0,infinity)` solves the differentiated backward equation
with potential b'=-z^2 and initial delta at zero. Its Feynman–Kac kernel is
strictly positive for t>0; localization makes this statement valid for the
nonexplosive cubic drift. Hence a is strictly increasing, not merely weakly
monotone. Smooth step approximations supply the distributional initial datum.

The nondegenerate quartic prior then gives

    Cov(A,H)=E Cov(A,H|Z0)+Cov(a(Z0),h(Z0))>0.

The second term is strictly positive by the independent-copy formula above.
Since the nine-outcome record contains H, it is not independent of the sign.
Writing d_r=p(PLUS,r)-p(MINUS,r), the balanced binary limit has

    r_9=1/2-(1/2) sum_r |d_r| < 1/2.

This supplies a rigorous **qualitative** advantage, not a decimal lower
bound. Finite-N ZERO mass vanishes and the no-record risk tends to 1/2;
therefore for sufficiently large admitted N the finite-record risk improves
its own no-record baseline as well. No numerical threshold N0 is obtained.

## 3. Pay conditional projection, not only Novikov

The source's limiting filter is correct in its stated conditional-law class.
Under the reference law Q let Z have the stationary cubic diffusion and
Y be independent Brownian with variance 2/k. For any r>1 its likelihood obeys

    E_Q[L_T^r | Z_path]
      =exp[(r^2-r) k/4 integral_0^T Z_s^2 ds].

Quartic tails and stationary Jensen bounds make this finite for every fixed
r,T. In particular E_Q L_T=1, and all these likelihood powers have finite
expectations. Conditional integration over Y shows that weighting by L_T
leaves the signal's marginal path law unchanged. The resulting measure is
the mathematical model law P, not new physical evidence.

These bounds plus Hölder and polynomial signal moments make the stochastic
integrals in the product/conditional-projection argument true martingales
on finite horizons. Doob's inequality controls likelihood maxima; compact
localization and these integrable bounds justify the limits. Thus for compact
C2 tests (and the stated integrable extensions), the Zakai equation and its
Itô quotient give

    d pi(f)=pi(Gf)ds+(k/2)[pi(zf)-pi(z)pi(f)] [dY-pi(z)ds].

The normalized innovation is Brownian in the observation filtration. It is
not the hidden W, and it is not independent of the signal by that fact alone.

The primary [Campillo–Corso paper](https://arxiv.org/html/2606.09272v1)
was checked at sections 5.3–5.5, 6.1–6.2 and 8.1–8.2. Its general uniqueness
result in section 7 assumes bounded first/second coefficient derivatives;
the cubic drift does not satisfy that premise. We do **not** invoke it for
arbitrary formal Zakai solutions. Uniqueness here is of the normalized
Kallianpur–Striebel conditional law for the frozen joint model, up to usual
versions. No finite-N microscopic filter limit is thereby proved.

## 4. Fix the initial normalization of the exact likelihood

For full hidden law mu at the observation origin and observed initial parity
u0, the unconditional-density starting row is

    alpha0(x)=mu(x) 1{o(x)=u0}.

Its mass is P(o0=u0), not one. If instead a normalized conditional prior is
used, multiply the history likelihood by that initial-parity probability.
The source's phrase “conditioned ... retaining the probability” is accepted
only with this explicit convention. Do not multiply or divide that factor
twice. Hidden correlations in mu are retained.

Within each no-visible-switch block, retain the **full original diagonal**.
Outgoing visible transitions are killing, not reflected rates. The source
ordered matrix product is then the exact density relative to initial-parity
counting measure and ordered switch-time Lebesgue measure. Future target
probabilities require the full hidden semigroup after the observation ends.
The numerical checker compares this construction with a distinct Poisson
uniformization/counted-switch calculation in a small actual q12 parent.

## 5. Tighten what the open atom actually requires

“Uniformly over every microscopic history” is stronger than necessary for
the bounded Bayes-risk target. An integrated L1 likelihood/decision comparison
with exceptional histories paid in probability is sufficient. It must be
uniform over the **decision rules** in the risk comparison; this is not a
pointwise supremum over every rare path.

Conversely, joint reference convergence and UI of the joint likelihood alone
do not automatically imply convergence of conditional likelihoods/risks.
A finite-cylinder approximation or equivalent deficiency argument must
show that the fine parity history supplies no extra target information.
This conditional-projection payment is part of the remaining gate.

## 6. Numerical anomaly retained, not repaired by symmetrizing

The supplied COARSE table has a sign-reversal cell residual about .01535246,
despite balanced total sign marginals. FINE and REFINED residuals are below
5e-15. Grid atoms at the inclusive/exclusive +/-1/2 quantizer boundaries can
break the continuum null-boundary symmetry; the coarse y grid has such atoms.
The asymmetry is not a violation of the exact diffusion theorem.

The intake does not average mirrored cells or retune thresholds to make a
diagnostic pass. All three risk decimals remain uncertified. The 27 printed
cell approximations are exactly the REFINED table, not independent bounds.
Mass renormalization, discarded y-grid tails and grid differences do not pay
spatial/time/advection/quantizer error. No decimal risk is promoted.
