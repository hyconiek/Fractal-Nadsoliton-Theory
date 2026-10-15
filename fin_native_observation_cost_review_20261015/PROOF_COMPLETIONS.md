# NOC audit: a native classifier, not a transferred Gaussian resource claim

Date:2026-10-15. Analytical model audit under inherited exact spectral and
MT/KO local-regime premises. Not mechanized or externally peer reviewed.

## 1. Native generator comparison is genuinely uniform in block fraction

The source normalized chart is correct. For R the exact coefficient sum is
`sum_a |1_B-alpha|=2k(1-alpha)`, so the drift remainder is O(sqrt(k)/N),
not an inverse-alpha blow-up. The global remainder uses native MT moments:
`sqrt(N) E|m| ||w||=O(N^(-1/4))` and `sqrt(N)E||w||^2=O(N^(-1/2))`.

The exact R conditional variance is bounded uniformly over1<=k<=N/2.
Mixed brackets are handled by the actual finite hypergeometric conditional
law. FF error uses p-u; RR error additionally uses block-composition noise;
FR has zero conditional mean and a concentration bound. No finite-N driver
independence is assumed. The source's conservative FR bound suffices even
though a sharper conditional calculation is available.

Splitting inside and outside jumps pays the absolute third compensator:
inside contributes O(k^(-1/2)), outside/global contributes O(N^(-1/2)).
Reflection deletes genuine outward rates without renormalization. A test
increment is at most C||Dv||/sqrt(k), so its stationary L1 boundary payment
is `(CN/sqrt(k))pi_D(B_D)`. The inherited positive boundary energy makes it
polynomial-times-exponentially small. All tests are defined also at outside
chart points; reflection is not silently ignored in a full-generator formula.

Thus the stated stationary L1 comparison holds for bounded derivatives1–3,
with constants uniform in N,k for the two frozen models. It is not a
statewise uniform approximation or a path-law Gaussian equivalence.

## 2. Corrector and semigroup domain are paid

The explicit OU semigroup formula is correct. Couple its deterministic-start
transition to its stationary Gaussian using the same noise: Lipschitzness
of the cosine and stable contraction bound each centered term by
`C exp(-gamma_min ell)(1+|f|+|r|)`, with a summable variance correction.
This proves pointwise absolute/local uniform convergence and linear growth.

Derivatives1–3 commute with the sum because each differentiated term is
uniformly bounded by a summable exponential. The same domination justifies
`P_s psi=sum_ell(P_(s+ell)g-mu)` and derivative bounds for0<=s<=1. Gaussian
transition moments and linear growth justify semigroup expectation/sum
interchange. Localized Itô or the explicit formula gives
`partial_s P_s psi=G P_s psi`; no arbitrary formal PDE uniqueness is needed.

The class used in Duhamel therefore has exactly the bounded derivatives and
growth needed by the native estimate. A bound on psi itself in sup norm is
neither available nor required.

## 3. Duhamel does not assume a Markov chart

On the full finite native state space differentiate
`P_D,1-s[(P_s psi)oPhi]`. Its derivative is
`P_D,1-s[(G P_s psi)oPhi-L_D((P_s psi)oPhi)]`. Integration yields precisely
the source sign convention for e_N. Stationarity/positivity contracts L1,
so the one-step defect costs C epsilon_N,k.

This proof is valid even though Phi=(F,R) is not generally an autonomous
Markov state. A separate intake diagnostic exhibits N3 states(0,1,1) and
(0,1,3) with the same chart but different finite observed drift. That floating
witness is not an interval no-go certificate; it demonstrates why the proof
uses the full native semigroup. A full q12 N2 matrix/OU quadrature independently
checks the nonzero Duhamel defect and its integral identity numerically.

## 4. Native concentration and threshold conventions

The exact discrete Poisson identity, stationary second moment of psi_N and
orthogonal martingale differences pay the source probability bound. No
constant full-chain gap, iid record assumption or growing-time OU limit is
used. Stationarity bounds the **average** defect without a factor M.

The wrong-decision event uses a weak inequality at the midpoint. Apply the
strict-tail bound at, for example, one quarter of the ideal gap to cover
ties as well; constants change but no exponent or conclusion changes.
There are M observed snapshots0,...,M-1. The terminal state Y_M appears
only in the proof and is not extra classifier data.

## 5. Transfer only the final bounded event to the original parent

An unambiguous sufficient transfer is: lift the actual late count law to
the full exchangeable labelled law, at TV error epsilon_prep,N. Compare
its future path with the stationary local reflected law by same-law initial
coupling and suppression of departures. The latter probability on M native
units is at most `CN M pi_D(B_D)` by the stationary outgoing compensator.

Therefore under each unchanged parent,

    P_actual(wrong) <= P_stationary_D(wrong)
                       +epsilon_prep,N+CN M pi_D(B_D).

This is enough for the frozen polynomial M and common late polynomial origin.
The source's more expanded outside/mixing/exit ledger is also a conservative
o(1) route. Initial exceptions and exit outcomes are retained; no survivor
postselection is used. No division of epsilon_prep by alpha and no multiplication
by M is licensed or necessary when transferring the final decision event.

## 6. Exact floors, ideal mean gap and arithmetic ceiling

The floor orders in the source are correct. The intake checks its finite
examples by integer16th/fourth roots, not floating exponentiation. Small-N
generator fixtures are not executions in the tiny MT trapping domain or
instances of the frozen asymptotic block when k>N/2.

[MEAN_GAP_CERTIFICATE.json](MEAN_GAP_CERTIFICATE.json) adds an independent
outward interval for the ideal first-order coefficient d, under the inherited
rational spectral and mpmath.iv assumptions. The mean-value theorem gives

    mu_A7(alpha)-mu_control(alpha)
       >=d alpha exp(-alpha beta_control),
    beta_control=(1/gamma_control-1)/2.

The certificate proves exp(-beta_control/4)>1/2, so a d/2 lower coefficient
works for0<alpha<=1/4. This is an ideal mean estimate, not native constants
or a useful N0. It does not endorse all70 printed producer digits.

The source point-tail check estimates an infinite tail using one final
floating difference; cancellation can make that difference zero. It is not
a rigorous bound. The actual proof's geometric estimates are valid. The
intake separately sums explicit bounds for variance-prefactor and phase tails,
certifying the stated four finite240-term fixtures without using cancellation.

## 7. Terminal theorem and limits

Combining these payments proves the source's asymptotic consistency for
actual finite-N parents along the frozen k_N,M_N. The polynomial upper terms
are N^(-1/8), N^(-3/16), N^(-13/32), plus native boundary and preparation errors.
Constants/N0 remain existential and inherited-premise dependent.

This is a one-sided sufficient result for a specific amplitude-sensitive
cosine test. It closes that native achievability bridge, not native necessity,
minimax optimality, sign-only cost or all Gaussian resource transfer.

## 8. Why the next task is a readout-scale bridge

A sharper native lower bound is interesting but not a blocker for using this
effective model. Existing incidence/locality obstructions in reports148,
154,200,202,203 are reused facts; repeating them is not a new physics path.

The selected next task removes one explicit operational restriction: the
cosine test assumes a known amplitude normalization. A simultaneous two-mode
angular event is invariant to a common nonzero multiplier and is already a
function of FIN's labels. Its native proof needs new nonsmooth/anti-concentration
payments. It is a different experiment and new plan, not a result of NOC.
No physical common-gain instrument, A7 source, locality or units are supplied.
