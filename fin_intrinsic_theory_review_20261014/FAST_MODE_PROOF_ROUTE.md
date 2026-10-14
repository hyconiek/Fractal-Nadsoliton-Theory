# A derived route through the TH-005 convergence gap

Date: 2026-10-14.
Status: **ARCHITECT_PROOF_STRATEGY**, with exact deterministic identities.
The stochastic bounds below are proposed lemmas, not intake-certified results.

## 1. Use the full occupation tangent, not an invented fast environment

At the accepted fold, set

    Sigma_f=diag(p_f)-p_f p_f^T,
    M_f=X^T Sigma_f X,
    H7=I/g_f-M_f,
    J=g_f Sigma_f A-I, restricted to 1^T dp=0.

Let v be the unit null vector of H7, in the source's fixed sign convention.
Define a right and left occupation null pair

    r_f=Sigma_f Xv,
    ell_f^T=g_f v^T X^T.

Then `Jr_f=0`, `ell_f^T J=0` and `ell_f^T r_f=1`, directly from
`M_f v=v/g_f`. Thus

    Q_f=I-r_f ell_f^T

is an oblique projection on the simplex tangent, with ten-dimensional range.
It is not a covariance/spacetime metric. Normalization has already removed
the constant occupation direction.

In coordinates `(theta,r_hidden)`, with R the orthogonal projection onto the
complement of `span(1,X)`, the deterministic equations are exactly

    theta_dot=-g grad Phi_g(theta),
    (r_hidden)_dot=R softmax(Xtheta)-r_hidden.

The Jacobian has diagonal blocks `-g_f H7` and `-I4`. Consequently its tangent
spectrum is one zero, six negative retained eigenvalues and four `-1`s.
The source fold inertia licenses stability; numerical eigenvalues are a
diagnostic, not how their signs are proved. Do not copy MP7's nonorthogonal
three-coordinate principal-block lower bound directly into an orthogonal
ten-coordinate kinetic gap. Establish the actual restriction/bound you use.

Let B be J restricted to `range(Q_f)`, in a fixed chosen basis. Its Hurwitz
property permits a positive quadratic Lyapunov form obtained from

    B^T P+PB=-I.

This is an analysis construction from a supplied drift, not a new stochastic
law, physical metric, extra state or fitted damping matrix.

## 2. Make the prepared/stopped contract quantitative

Fix a compact lambda range, finite scaled horizon S and finite soft radius R.
Use the exact fold, not its decimal midpoint:

    g_N=g_f+lambda N^(-2/3),
    theta_N=g_N X^T p_N,
    Z_N(s)=N^(1/3) v^T(theta_N(N^(1/3)s)-theta_f),
    w_N=Q_f(p_N-p_f).

One sufficient **target** initial class is

    Z_N(0) => Z_0, with initial support strictly inside (-R,R),
    E||w_N(0)||^2 <= C/N,
    initial states in the interior local domain with probability tending to 1.

Freeze a support/moment or bad-initial-event condition strong enough for all
estimates; "prepared" alone is not a theorem hypothesis. This is restriction
to an initial-law class, not an added microscopic law or physical preparation.

The class is nonempty: for bounded deterministic z0, round to actual counts
the positive vector `p_f+N^(-1/3) z0 r_f`. It sums to one because `1^T r_f=0`;
balanced rounding changes occupations by O(1/N). Then `w_N(0)=O(1/N)` and
the actual g_N-coordinate has `Z_N(0)->z0`. This construction should be proved
and exported, not used as an unspecified ideal preparation.

Stop on `|Z_N|>=R` or `||w_N||>=delta_N`, with `delta_N=N^(-1/4)` as one
candidate auxiliary tube; also stop at S. Norms and threshold constants must
be fixed. Show the auxiliary tube exit is negligible before the soft exit;
only then remove it from the observable limiting statement.

Why not the source's fixed original neighbourhood alone? Its soft boundary
scales to infinity. That would require a separate explosion/reinjection/global
matching theorem. A compact-z stopped law is the minimal defensible first bridge.

## 3. A sufficient generator-energy bound

Write x=`ell_f^T(p-p_f)`. On the stopped domain, x=O(N^(-1/3)); its difference
from the actual soft theta/g_N convention is O(N^(-2/3)). Keep that conversion
explicit in the final Taylor estimate.

The smooth finite-alphabet drift has, on the fast subspace,

    Q_f b_N(p,g_N)
      =B w+O(x^2+|g_N-g_f|+|x| ||w||+||w||^2+1/N).

There is no linear x forcing because `Jr_f=0`. The finite-N difference is
the genuine leave-one-out correction, not a new empirical-refresh convention.
For `V=w^T Pw`, bounded jumps and smooth rates should give, after absorbing
local cross terms into the stable drift,

    L_N V <= -kappa V+C(x^4+|g_N-g_f|^2+1/N).

All constants must be independent of N, within the frozen R/S/lambda scope.
Prove this inequality with actual generator jumps and an explicitly chosen
local domain. Do not justify it by saying fluctuations "should be N^-1/2".
If source interval gaps are too conservative, existence of a positive gap
can suffice for an asymptotic theorem; do not manufacture certified numerics.

For physical model horizon `T_N=N^(1/3) S`, stopped Dynkin estimates would
then yield

    E integral_0^(T_N wedge tau_N) ||w_N(t)||^2 dt = O(N^(-2/3)),
    E integral_0^(S wedge tau_N/N^(1/3)) ||w_N(N^(1/3)s)||^2 ds = O(1/N).

Use stopped/killed occupation estimates, not an unjustified exponential bound
on `E V(t wedge tau)` after a nonzero absorbing exit value. The nonnegative
stopped energy estimate also targets

    P(fast tube exits before soft exit/horizon)
      <= C N^(-2/3)/delta_N^2 = O(N^(-1/6)).

Handle initial tail events, jump overshoot and norm-equivalence constants.
These displayed rates are sufficient proof targets, not computed finite-N
confidence bounds.

## 4. Use the bound directly in the soft martingale problem

For `|Z|<=R`, the projected drift remainder contains `|x|||w||` and `||w||^2`.
The energy bound and Cauchy–Schwarz target

    E N^(2/3) integral |x| ||w|| ds = O(N^(-1/6)),
    E N^(2/3) integral ||w||^2 ds = O(N^(-1/3)).

Terms `|x|^3`, `|g_N-g_f||x|`, `|g_N-g_f|||w||` and the LOO O(1/N) error also
vanish in the required integrated scaling. Compute them explicitly rather
than writing an unbounded o(1) that depends on a path still to be controlled.

Because `theta_dot` closes exactly at leading order, hidden occupation modes
do **not** create an extra leading theta drift. They do affect quadratic
variation. The vanishing auxiliary-tube event and the full occupation control
give `p_N->p_f` uniformly on stopped intervals, sufficient to replace the
actual soft covariance by `2g_f`. Use conditional predictable brackets;
stationary covariance or marginal histograms are insufficient.

The rescaled jumps have size O(N^(-2/3)); their absolute third characteristic
is O(N^(-2/3)). For every smooth compactly supported test function, Taylor
remainder then vanishes. Tightness and identification should give the killed
generator

    G phi(z)=-g_f(a lambda+b z^2/2) phi'(z)+g_f phi''(z)

on `(-R,R)`. The nonzero diffusion and smooth bounded-domain drift supply
uniqueness of the stopped one-dimensional law, once its precise absorption
convention and boundary mapping are defined and verified. Prove hitting-map
continuity under the limit; do not infer exit-time convergence solely from
finite-dimensional convergence of coordinates.

This route may avoid a general memory-pole or perturbed-test-function campaign.
If direct bounds fail, export the offending term/rate and only then decide
whether a corrector is needed. The route itself remains **unexecuted**.

## 5. Conditions that would refute this candidate scope

- A second nonnegative eigenvalue in the matching full occupation chart.
- A nonvanishing integrated drift/noise correction under the paid preparation.
- Auxiliary tube exit with nonnegligible limiting probability.
- A clock or feature normalization needed to fit, rather than derive, coefficients.
- Uncontrolled initial layers or explosion hidden by a stopped/global scope swap.

Retain the accepted path-LDP theory if this convergence scope fails. A failed
reduction need not invalidate the microscopic parent or justify adding a new law.
