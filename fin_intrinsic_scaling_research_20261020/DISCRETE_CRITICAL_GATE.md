# IS-004: exact critical defect interface for the actual discrete source

Date: 2026-10-20. Exact interface plus a sufficient conditional criterion.
The criterion's actual A7 hypotheses below are NOT proved in this packet.

## 1. Do not insert CTMC Gamma into a discrete martingale

Reflect the already selected global source on D_N and start exactly pi_D.
Let d=d_N, h_N=1/(d sqrt(N)), Z=N^(1/4)m and v_n be its true global
scaling. On internal count edges, P_D(n,n')=Q_N(n,n')v_n v_n'/d.
Define the actual per-critical-time drift and raw quadratic coefficient

    beta_N(n)=N^(3/4)sum_(internal a,b) p_a R_ab^N v_n v_n' (sigma_b-sigma_a),
    G_N(n)=sum_(internal a,b) p_a R_ab^N v_n v_n' (sigma_b-sigma_a)^2.

For the next single intrinsic step,

    E[Delta Z|n]=h_N beta_N,
    E[(Delta Z)^2|n]=h_N G_N,
    h_N^(-1)Var(Delta Z|n)=G_N-h_N beta_N^2.         (1)

This last subtractive term belongs to the actual DISCRETE martingale.
Q=d(P_D-I) may be used as an analysis operator but not as a replacement
path law with an automatically identical predictable bracket.

Since a source step changes at most one label,

    |Delta Z|<=2N^(-3/4),
    |beta_N|<=2d N^(-1/4),
    0<=G_N<=4d/N<=4q.                               (2)

## 2. Two exact unpaid stationary payments

Freeze f(z)=q(chi lambda z-z^3/3), with the inherited d-normalization;
do not fit the coefficient after inspecting the answer. A sufficient next
pair of estimates is

    D_N=E_piD[(1+|Z|)|beta_N-f(Z)|] ->0,
    C_N=E_piD|G_N-2q| ->0.                          (3)

The Gibbs moments are already supplied by MT. These are actual v_n v_n'
defects of the unaugmented source, not consequences of its new polynomial
lower bound or near-uniform local fields alone. No pointwise v=1, native
statewise parity-drift equality or arbitrary kinetic rigidity is assumed.

If the first estimate holds, E|beta_N|=O(1), and (2) pays the discrete
centering correction without demanding a new L2 drift estimate:

    h_N E beta_N^2<=2N^(-3/4)E|beta_N| ->0.          (4)

Thus (3) really gives the predictable bracket limit 2q s, not merely the
raw conditional second moment. The bound is in the stationary comparison;
unbounded iid moments are not transferred through preparation TV.

## 3. Why these payments would suffice

Let h=h_N. On the critical grid Z_j=Z(Y_j), the exact martingale has
increments Delta M_j=Delta Z_j-h beta_N(Y_j), bounded by 4N^(-3/4), and
predictable variance h[G_N-h beta_N^2]. Its bracket is uniformly bounded
on fixed horizons by (2). Fourth-power martingale expansion and the
max-jump-squared control of the fourth compensator give uniform fourth
moments, as in the paid SQ proof, now with centered discrete increments.

For compact containment use the stopped discrete identity
E[Delta(Z^2)|n]=h[2Z beta_N+G_N]. The target 2Z f(Z) is bounded above;
stationarity and (3) control the integral of |Z(beta_N-f)| even with the
stop indicator. The nonnegative stopped terminal square then bounds the
probability of reaching a large radius, including the jump overshoot.
No unjustified stationary law at the stopping time is used.

On compact stops the target drift is bounded and total residual tends to
zero in L1. Bounded brackets and vanishing jumps pay tightness. The drift
Riemann sum differs from its piecewise-constant integral by at most
h max|beta_N|<=2N^(-3/4); bracket endpoint error also vanishes. Equations
(3),(4), fourth-martingale UI and cylinder tests pass both martingale
identities in the joint natural past filtration. The same unique
nonexplosive SDE as SQ follows, with Brownian increments independent of
the initial quartic variable.

This is a conditional application of the already used martingale method,
not a new generic diffusion theorem and NOT an executed SDE theorem for
P_N. If (3) fails, retain the actual defect/outcome; do not change f, the
clock, Gibbs measure, source postulate or observation to force agreement.
The new local window would transfer a PAID bounded-path result to late iid
origins, but cannot supply the unpaid estimates (3) itself.
