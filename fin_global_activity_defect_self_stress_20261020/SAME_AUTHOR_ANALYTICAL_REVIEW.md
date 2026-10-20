# Same-author full-input analytical recheck of AD

Date: 2026-10-20. **SELF_REVIEW_SUPPORTED_CONDITIONAL_ARGUMENT; NOT INDEPENDENT.**
The user requested continued work without delegation. Actual AD source
files were read; original hashes and statements remain unchanged. This
note records the scope checked and one explicit grid-time payment, not a
new source/physics theorem or a proof-assistant certificate.

## 1. Global trial and faces

Rewriting Q_N=n_a exp(E_ab+ell_ab/N) and d_N=Nr+1 gives exactly the
displayed F_N row, including its single self contribution. Positive
tau=sum p_a exp(E_ab) is bounded below by r exp(-H) on the whole simplex.
It therefore has a positive fixed-neighbourhood extension uniformly over
the compact gain interval. This supports a0,a1 smoothness even on faces;
there is no hidden 1/p_a factor. Terms with unavailable departures vanish.

The Taylor numerator coefficients are tau a0^2=r and
2tau a0a1+(1+kappa)a0^2+a0 b dot grad a0=1. Thus the denominator r+1/N
cancels through first order. The O(N^-2) remainder needs bounded second
derivatives of a0 and first derivatives of a1 on that fixed neighbourhood,
which their explicit smooth formulas supply. A numeric all-simplex C or
N0 is not exported. The four native interval tests are supporting controls,
not the source of that analytical compactness argument.

## 2. Local comparison and weights

Ordered internal pairs and self edges give exactly the factor-one-half
log-flow identity. Margins of the auxiliary table differ by pi*delta;
they are not silently replaced by pi. IS's global v bounds give the
needed logarithmic bound before localization. Both exact and auxiliary
outward flow sums are controlled by actual departure-state mass. Hence
the proof obtains its conditional Hellinger bound directly, rather than
dividing a global error by exponentially small pi(D).

Weighted Cauchy--Schwarz uses the ACTUAL row mass bounds to absorb
(1+|Z|)^2 into MT moments. The drift prefactor is 2d_N N^-1/4, of order
N^3/4, whereas raw noise has prefactor 4d_N/N, of order one. These factors
and the O(N^-2) row residual pay the displayed drift/noise powers. The
source's pointwise v need not converge uniformly. The comparison with
SQRT uses actual internal off-diagonal rates, not parity lumpability.

## 3. Discrete bracket and grid clarification

Write h=1/(d_N sqrt(N)), M_N(t)=M_floor(t/h), in the piecewise-constant
natural filtration. Every predictable bracket step is at most 4q h.
Thus, pathwise for ANY real times s<=t,

    <M_N>_t-<M_N>_s<=4q(t-s+h).                     (1)

This is the correct grid version of a bounded bracket density. For bounded
stopping times tau, optional isometry consequently gives

    E|M_N(tau+theta)-M_N(tau)|^2<=4q(theta+h),       (2)

after the standard finite-horizon localization. The explicitly embedded
filtration is used throughout; a continuous bracket density is not asserted.
The extra h tends to zero. Fourth-moment UI is paid by
bounded conditional variance and |Delta M|<=4N^-3/4. Compact-stopped drift
increments use bounded f, the vanishing total residual and a single-grid
endpoint error <=2N^-3/4. Together (1),(2) make the tightness step explicit.

The centered variance is G-h beta^2, with h E beta^2 controlled by
2N^-3/4 E|beta|. Weighted drift convergence gives this L1 drift moment;
no unproved L2 drift bound is required. The stopped square identity keeps
its positive raw quadratic term and terminal square. Martingale and
square-minus-bracket identities pass with the random initial state in the
joint natural past filtration, not just via marginal covariance limits.

## 4. Preparation and disposition

The late iid coupling uses the already paid IS window of THIS source,
then the same source transitions. Its initial error is paid once; the
stationary reflected comparison's suppressed-event probability is at most
the number of future steps times pi_D(B_D). These are bounded-path transfers,
not automatic unbounded moment/entropy or global-equilibrium statements.

No concrete analytical contradiction was identified in these checked
steps. The argument remains conditional on admitted spectral/MT/source
premises, self-audited and not mechanized. The vertex interval-positivity
noncertificate remains retained; it neither supplies a numeric global N0
nor refutes the compactness/existence proof. Independent review has not
occurred. All natural F01 and physical/ontology gates stay open.
