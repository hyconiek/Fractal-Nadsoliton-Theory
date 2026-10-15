# DKO audit: boundary payments, joint path law and Gaussian resource scope

Date:2026-10-15. Analytical model audit under inherited exact spectral,
MT/KO local-regime and labelled-coupling premises. Not mechanized or peer reviewed.

## 1. Do not use a full-generator martingale for reflected paths unchanged

The source correctly derives characteristics from the original heat bath.
Its stationary reflected proof comparison needs this additional payment.
Let B_N be the actual departure layer and pi_D the lifted local Gibbs law.
The number of suppressed attempts has expected rate at most C N pi_D(B_N),
with pi_D(B_N) bounded by a polynomial times exp(-Nb).

At fixed alpha>0, an F/U jump has size at most C_alpha/sqrt(N). Hence the
omitted reflected drift is bounded by C_alpha sqrt(N)1_B and the omitted
bracket by C_alpha 1_B. Their integrated stationary expectations vanish on
every fixed fast horizon. The soft coordinate jump is O(N^(-3/4)), whose
omitted drift is O(N^(1/4)1_B), also negligible. Alternatively stop the full
parent at its first departure and couple; do not silently identify L with L_D.

The actual iid-start initial-to-local labelled TV error and boundary exits
remain errors of the path comparison. No actual trajectory is postselected
or physically reflected.

## 2. Native occupation moments pay the fixed-horizon residuals

MT gives E m^2=O(N^(-1/2)) and E||w||^2=O(1/N) under pi_D. Thus

    E[|m|sqrt(N)||w||] <= sqrt(N)*sqrt(E m^2 E||w||^2)=O(N^(-1/4)),
    E[sqrt(N)||w||^2]=O(N^(-1/2)).

Stationarity turns these into L1 bounds on integrated global drift errors
over fixed T. Private errors are uniformly bounded for fixed alpha. Sampling
errors are bounded conditionally on complete counts using exact block
hypergeometric moments. These moment calculations are done in the stationary
comparison; only probabilities/path laws pass to the actual parent via TV.

For general Fourier features the uniform old-label jump product is
`delta_rs+f_r(old)f_s(old)`, not the constant2 of a parity sign. Conditional
sampling pays the global/private cancellation in expectation and integrated
probability, not pointwise in each block. The intake's finite fixtures
explicitly reject that shortcut.

Freeze notation: `Zhat_N(t)=N^(1/4)m_N(t_N+t)` on **fast** model time.
It is not the UF critically accelerated process evaluated at t_N+t. The
parity-line drift cancellation bounds its drift expectation by O(N^(-1/2))
and its bracket density by O(N^(-1/2)) on a fixed fast horizon, so its change
vanishes while its initial quartic amplitude remains random.

## 3. Complete tightness and the joint martingale identification

Use the semimartingale equations for F and U. On a fixed horizon, Gronwall
or variation of constants bounds their suprema in probability by their tight
initial states, bounded-bracket martingale suprema and vanishing integrated
errors. A uniform moment bound at deterministic times alone would not control
maxima at stopping times; this equation supplies that missing step.

The linear drift increments then satisfy Aldous' criterion after localization;
the global integrated residual controls every subinterval. Martingale brackets
have uniformly bounded densities at fixed alpha and jumps vanish. The
fourth-jump compensator is at most maximal-jump-squared times the quadratic
compensator. The fourth-power martingale formula/BDG supplies uniform fourth
moments, hence UI of the martingale and square-minus-bracket. Mere bounded
expected quadratic variation would not suffice for every passage claimed.

Pass the identities jointly with the initial vector in the natural past
filtration. The six-component continuous limiting martingale has deterministic
diagonal bracket diag(2I,2c_alpha I), so its Brownian coordinates and future
increments are independent of the initial law. The unique linear SDE law
and the paid product initial Gaussian/quartic law identify the entire limit.

Therefore the fixed-fast-horizon joint global/private OU theorem is accepted
pointwise for each fixed alpha in(0,1). The observed V is a stationary Gaussian
sum of two rates, generally non-Markov. No full microscopic observation-law
TV equivalence or uniform growing-horizon theorem follows.

## 4. Bounded signs and the exact gauge scope

The source's angular proof of the Gaussian sign identity is correct. The
two-time pair is nonsingular for h>0 and finite positive gamma, so finite-N
zero outcomes vanish in the limit without deletion. Bounded sign expectations
and the frozen agreement event receive the actual-parent weak law cleanly.

With q from the null mode, the derivative of
`rho(gamma)=q+(q^gamma-q)/(1+c gamma)` is strictly negative on gamma in(0,1).
Its endpoints approach1 and q, proving unique inversion when rho lies in(q,1).
Only ensemble probabilities are identified, not gamma from one observed bit.

Fixed nonzero per-channel amplitude multipliers cancel exactly in sign products.
The source amplitude-verification flag repeats its own formula; it is not an
independent gain test. The intake tests the actual pointwise sign identity
and the failure under a time-changing sign gain. Offsets, distortion, sign
errors and channel-dependent clocks remain outside the theorem. A common
clock-unit gauge is not a physical clock calibration.

The [sign-event certificate](SIGN_EVENT_CERTIFICATE.json) adds independent
outward probability bounds from the inherited rational spectra. Its positive
arcsin series uses110 terms and a geometric tail, with exact raw-dyadic
outward decimal rounding. It assumes mpmath.iv correctness, not physical
evidence or a finite-N accuracy certificate.

## 5. Gaussian information and sign achievability are correctly separated

For R_q(i,j)=q^|i-j|, its finite-section quadratic form is the Fourier
integral with symbol `(1-q^2)/(1-2q cos(theta)+q^2)`. This proves the source's
uniform eigenvalue bounds for every finite record length, not just a
numerically sampled Toeplitz matrix. Absolute exponential summability pays
D1,D2; d0!=0 pays the strict Frobenius lower bound.

The relative covariance E is symmetric and has norm<=1/2 on the stated
alpha interval. Taylor's theorem for x-log(1+x), followed by eigenvalue
comparison, gives the source's explicit c M alpha^2 and C M alpha^2 KL bounds.
The constants depend on the two fixed gammas and fixed spacing h.

For sign achievability use the standardized latent pair (F,U0), with U0
variance1 and observation W=sqrt(alpha)F+sqrt(1-alpha)U0. The latent stationary
OU law and L2 gap are independent of alpha. For two bounded lag-block
functions, condition on the intervening Markov endpoints and apply Hermite
L2 contraction across the gap between their intervals. This pays the claimed
summable segment covariances and O(1/M) empirical variance uniformly in alpha.
The nonzero linear-in-alpha sign-mean difference and midpoint test then give
vanishing error when M alpha^2 grows. Large KL alone is not used to claim it.

Bookkeeping correction: M amplitudes yield M-1 adjacent sign products. The
source's estimator with M products uses M+1 amplitudes. Either convention
has the same asymptotic scaling, but data contracts must use one consistently.

This theorem belongs entirely to the limiting Gaussian family. It is not a
simultaneous alpha_N->0, M_N->infinity theorem of the finite-N FIN parent.
It also assumes the frozen common spacing/clock model, not arbitrary
composite physical nuisance classes.

## 6. Exact terminal admission and genuine remaining question

The previously deferred fixed-fast-path bridge is now paid in its local,
fixed-alpha/fixed-horizon scope. The dynamic sign law and ideal gauge
identifiability are paid, and the Gaussian resource theorem is paid separately.
No finite-N triangular resource theorem, physical source, apparatus,
locality/conservation, quantum/gravity or theory-of-everything result follows.

A direct finite-parent bounded-statistic proof can be a more focused next
target than full growing-horizon path equivalence. The next plan proposes
a **different bounded amplitude statistic**, not a silent extension of the
sign-only theorem. Its corrector and native error bounds remain research
targets; none is admitted by this review's arithmetic checks.
