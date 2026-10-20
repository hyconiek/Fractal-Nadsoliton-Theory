# AH proof audit: filtration, common bounds and arithmetic

Date: 2026-10-20. Conditional analytical mathematics, not formal proof-
assistant verification or external peer review. Original source is preserved.

## 1. Frozen class and the common-bound quantifier

Admit the source's one-label, twelve-state, symmetric radial true-jump chain,
uniform initial distribution and ideal equality/first-appearance record.
For AH-002, make the quantifier explicit:

    0<a<= min{J_ij,K_ij: i!=j}, a<=1/11,
    0<delta<=d_iso(J,K), 0<delta<=1, 0<eta<1.

The same a must apply to **both** supplied hypotheses for a single L to
bound each model's error. If separate lower bounds are available, take
their minimum (or the larger of the resulting valid lengths). A lower
bound for J alone does not pay the written concentration proof under K.
This is the explicit admitted scope of the source's generic class-bound
wording, not a bound inferred from observed frequencies.

The estimator and classifier remain analysis objects on an assumed equality
record, not an independently derived observer or physical apparatus.

## 2. Pay the martingale inequality without an unread theorem import

For a martingale difference D with conditional mean zero and |D|<=1,
convexity of exp gives, for real theta,

    exp(theta D)<=((1+D)/2)exp(theta)+((1-D)/2)exp(-theta).

Conditional expectation is at most cosh(theta). Since
`log cosh(theta)<=theta^2/2` (integrate tanh(theta)<=theta for theta>=0
and use evenness), iterating conditional expectations gives

    E exp(theta sum_(k=1)^n D_k)<=exp(n theta^2/2).

Markov's inequality and theta=u/n imply the one-sided bound
`Pr(sum D_k>=u)<=exp(-u^2/(2n))`; replace D by -D for the other side.
Thus the two-sided bound has factor 2, with **no independence assumption**.
For a fixed bound |D_k|<=c, rescale by c. For u=n epsilon the exponential
tail is summable. Subadditivity of the union over n>=M, followed by M->infty,
shows almost-sure convergence of the empirical martingale average to zero
(first for epsilon=1/r, then all epsilon>0). These elementary arguments pay
the source's Azuma/Borel--Cantelli use, including its constants.

## 3. AH-001: coverage, consistency and random names

As long as a particular vertex has not appeared, the next-state probability
of that vertex is at least min off-diagonal J. Therefore the geometric
coverage tail in the source is valid; a union over twelve states pays finite
cover time almost surely. Uniform initial data do not hurt this estimate.

For every i,j, there are ten two-step intermediate vertices if i!=j and
eleven if i=j. Hence `J^2_ij>=b=10a^2`. With `beta=12b`, one has
`0<beta<=120/121<1`. The residual
`R=(J^2-beta U)/(1-beta)` is nonnegative stochastic, and is doubly stochastic
because J is symmetric and U is uniform. Thus uniform pi is invariant, and
the source's geometric bound on centered state functions holds (even a
slightly smaller constant would suffice). This pays uniform convergence
of `g=sum_(r>=0) J^r f` and the Poisson identity `g-Jg=f`.

The displayed telescoping decomposition has the correct sign. Its
martingale differences are bounded by `2||g||_infty`, so section 2 pays
state-frequency convergence. The edge increments
`1_{Y_t=i,Y_(t+1)=j}-J_ij 1_{Y_t=i}` have absolute value <=1 and conditional
mean zero in the original-state filtration. Their convergence gives
`N_ij/N_i -> J_ij`, since `N_i/n ->1/12`.

Use a **simultaneous** event of convergence for all 144 entries, and only
then apply the path's first-visit bijection. No conditioning on that random,
path-dependent naming is needed or licensed. Once formed, the bijection
never changes. Define all finite estimators by padding unvisited rows with
zeros, making them measurable. After coverage their limit is in the finite
permutation orbit of the true matrix. The equality of all cylinder laws and
uniqueness on the generated sequence sigma-field then gives the source's
orbit-intersection argument. Conversely conjugacy with uniform initial law
preserves every equality record. All steps are valid in the frozen class.

## 4. AH-002: every time index, denominator and constant

Put `n=2m+1`, `L=n+1=2m+2`. Observe states Y0,...,Yn. Count **departures**
at t=0,...,n-1; do not count the terminal state as an extra departure.
All even states Y2,...,Y2m are departures and their following transition
is inside the observed prefix. With `F_(2k)=sigma(Y0,...,Y_(2k))`,

    E[1_{Y_(2k)=i}|F_(2k-2)]=J^2(Y_(2k-2),i)>=b.

Subtract these conditional means. A deficit of the even-time count below
bm/2 forces a martingale deficit below -bm/2. Section 2 pays
`exp(-b^2 m/8)` and the twelve-row union costs `12 exp(-b^2 m/8)`.
Since m>=2/b, all twelve departure counts are positive on the good event.
Merely seeing twelve names, with one only at the final state, would **not**
pay those denominators; the proof uses the stronger departure event.

For each ordered off-diagonal original pair, the edge-martingale sum M_ij
has n differences bounded by 1. At `u=delta b m/8`, and `n<=3m`,

    u^2/(2n)>=delta^2 b^2 m/384.

There are exactly 132 such pairs and the bound is two-sided, giving 264
times that exponential. On the simultaneous good event,

    |N_ij/N_i-J_ij| <= (delta b m/8)/(bm/2)=delta/4.

The diagonals are exactly zero. Union bounds are taken over the original
fixed labels **before** applying the random first-visit permutation; hence
the same max-entry bound in a random naming is valid without a false
renamed-filtration martingale assertion.

The lengths in the source make the two bad-event contributions at most
eta/2 each. The permutation action is an isometry for the max-entry norm,
so its finite-orbit quotient distance satisfies the triangle inequality
even when the empirical matrix is not radial. It is within delta/4 of
the true orbit and at least 3delta/4 from the other. Missing-row fallback
is used only outside the paid good event. The nearest-orbit rule therefore
has each model's error <=eta under the **common** a bound. Its existence
does not provide a fast algorithm for the permutation minimization.

## 5. Independent finite argument for completeness

There is a second way to check the main uniqueness implication. If J and K
are not conjugate, the minimum over their finitely many orbit distances is
some delta>0. Strict positivity provides a common a>0. Apply AH-002 with
eta=1/4. If the length-L anonymous laws were identical, every deterministic
decision would have the same probability of selecting J under both models.
The guarantee would require this same probability to be >=3/4 under J and
<=1/4 under K, a contradiction. Thus some finite prefix separates every
nonconjugate pair. This checks completeness through the finite probability
estimate rather than relying only on the infinite estimator/cylinder proof.
It is not an optimal uniform prefix-length theorem.

## 6. Exact arithmetic for the displayed illustrative budget

The source's `math.log`/floating `ceil` is a numerical preview, not an
outward certificate by itself. The new certificate uses rational inputs
`a=1/20`, `delta=1/10`, `eta=1/20`, so `b=1/40`.

For x>=1 write `x=2^e y`, `1<=y<2`. With `z=(v-1)/(v+1)`,

    log(v)=2sum_(k=0)^(N-1) z^(2k+1)/(2k+1)+remainder,
    0<=remainder<=2 z^(2N+1)/[(2N+1)(1-z^2)]

for v=y and v=2. Integrating the positive geometric series for 1/(1-z^2)
proves both bounds; z<=1/3 makes 32 terms amply sufficient. Exact Fraction
intervals enclose log(480) and log(10560), then every threshold in m. The
lower and upper maximum thresholds have the same ceiling. Therefore

    m=569231067, L=1138462136

are rigorously certified for this formula. The checker rejects m-1.
This strengthens arithmetic only, not theoretical sharpness or feasibility.

These parameters do **not** apply to the source's separate weights
`(1,2,3,4,5,6)` path fixture, whose true minimum entry is 1/36<1/20.
They can, for example, be instantiated by the flat radial jump matrix
and weights `(1,1,1,1,1,10)`: their common minimum is 1/20, and their
isomorphism distance is exactly 9/22>1/10 because the flat matrix is
permutation invariant. This is an exact analysis control, not a selected
physical comparator or a simulation of a billion-state record.

## 7. Hierarchy and nonconclusions

Clarify AH-003's compressed table: **all-lag** return traces determine the
finite spectral multiset in this class; a single return probability at one
lag is just one scalar, not a spectral determination. Labelled modal data
at a known lag can identify rates in that supplied comparison clock; an
uncalibrated lag or event-index data leave activity scale free. The earlier
mean-bracket inverse remains a distinct labelled one-hot theorem.

Do not promote the sufficient a-dependent length into necessity/optimality,
an efficient orbit algorithm, native N2 observability, physical seconds,
an equality detector, a nontrivial history source or emergent geometry.
F01/F03/F06 and the separate RP004 audit obligation are unchanged.
