# SRN audit: native angular consistency and its precise boundary

Date:2026-10-15. Analytical audit under the inherited exact spectral and
MT/KO local-regime assumptions, not mechanized or external peer review.
The original submitted source files remain unchanged.

## 1. A complete maximum-atom route without an unread asymptotic theorem

The source's conditional sampling law is correct. Given total counts and
all nonselected block counts, X has mass

    p(x)=C(A,x)C(B,t-x)/C(H,t), H=A+B.

The good event is t>=c k and t<=3H/4. Both selected populations have at
least cN entries, so A/H,B/H>=c. Chebyshev gives bad probability C_c/k:
E t=kH/N>=2ck and E t<=H/2, and Var(t)<=k/4. These statements precede
conditioning on the remaining block counts; no conditional good-event
probability is assumed afterwards.

For completeness a direct mass-ratio argument pays the source's abbreviated
Stirling paragraph. On the support,

    q(x)=p(x+1)/p(x)=(A-x)(t-x)/((x+1)(B-t+x+1)).

This decreases in x. An interior mode m is within a fixed constant of
tA/H, by solving q(x)>=1; the usual expression is
floor((t+1)(A+1)/(H+2)), with either adjacent mode allowed at equality.
For sufficiently large t, all four cells m,t-m,A-m,B-t+m are >=c' t.
The last two use H-t>=H/4, not an infinite-population approximation.

In the integer band |x-m|<=sqrt(t), reduced by a fixed constant factor if
necessary, the four cells remain >=c''t. The continuous logarithm of q has
derivative

    -1/(A-x)-1/(t-x)-1/(x+1)-1/(B-t+x+1),

whose absolute value is <=C_c/t there. Since q(m-1)>=1>=q(m),
|log q(m)|<=C_c/t. Summing log ratios gives

    p(m+ell)>=exp[-C_c(ell+ell^2)/t]p(m)

for 0<=ell<=c'''sqrt(t); the negative side is treated identically.
There are at least c''''sqrt(t) such support points. Normalization implies
p(m)<=C_c/sqrt(t)<=C_c/sqrt(k). Bounded t is absorbed by increasing C_c.
This is a uniform analytical proof; the finite fixtures do not prove it.

For h_plus=f3+f1, labels0/6 differ by4sqrt(2); for h_minus=f3-f1,
labels2/8 differ by3sqrt(2). Conditional on the other counts, an interval
of radius eta in the normalized sum contains at most C(1+eta sqrt(k))
values of X. Multiply by the maximum atom and add C_c/k. Thus

    P(|W3+/-W1-z|<=eta | n)<=C_c(eta+k^(-1/2))

uniformly in shift, permitted count population and 1<=k<=N/2. Ties are
included. The positive MT cylinder is essential; do not extend this theorem
to populations missing a selected type.

## 2. Four-coordinate generator: reuse the method, pay joint terms

Use the full two-block count chain, not an assumed Markov process in
Phi=(F1,R1,F3,R3). Conditional block assignments are hypergeometric under
the lifted local Gibbs measure. Every F/R moment needed here is bounded
uniformly in k by MT fast moments and exact sampling variance.

For each feature pair r,s the exact old/new jump product is
(f_r(j)-f_r(i))(f_s(j)-f_s(i)). At uniform reference it averages over j
to delta_rs+f_r(i)f_s(i), not a pointwise constant. Averaging old labels
and paying conditional block fluctuation bounds gives the displayed
diagonal OU diffusion and the vanishing cross-mode/global-private errors.
The intake explicitly checks nonzero finite cross brackets; independence
is not a native premise.

The global drift uses cancellation on the parity line plus
sqrt(N)E|m|||w||=O(N^(-1/4)) and sqrt(N)E||w||^2=O(N^(-1/2)).
The null feature has A f1=0, not the original Laplacian's mode1 eigenvalue.
Private forcing uses sum_a|1_B-alpha|=2k(1-alpha); it costs O(sqrt(k)/N).
Splitting inside/outside jumps pays the absolute third compensator.
Reflection deletes actual outward increments; its stationary payment is
(CN/sqrt(k))pi_D(B_D), not zero by definition.

Taylor expansion therefore yields the submitted stationary L1 bound with

    epsilon_N,k=N^(-1/4)+k^(-1/2)+(N/sqrt(k))pi_D(B_D).

The comparison is uniform over k<=N/2 for these two frozen matrices and
bounded derivatives1–3. It is neither statewise nor path-law equivalence.

## 3. Corrector, bandwidth and boundary replacement

H_eta is an analysis-only compactly supported convolution. Its difference
from raw H is confined to the two boundary strips. Native anti-concentration
pays E_piD|H-H_eta|<=C(eta+k^(-1/2)); independent Gaussian strip densities
pay E_G|H-H_eta|<=C eta, uniformly in alpha<=1/2.

The reference OU's contraction rate is bounded below by the two fixed
positive gammas. Lipschitz coupling to its stationary initial law proves
local uniform summability of the centered corrector terms and linear growth
C eta^(-1)(1+||x||). Differentiating the linear OU transition gives the
summable bounds C eta^(-j)exp(-j gamma_min ell) for j=1,2,3.
These bounds and Gaussian moments also justify expectation/sum interchange,
P_s psi and its generator equation for0<=s<=1; there is no unpaid global
zero-order sup bound.

Full-state Duhamel, with the sign convention checked in the NOC audit,
costs C eta^(-3)epsilon_N,k. Native stationary E psi(Phi)^2<=C eta^(-2)
pays the skeleton martingale and endpoint terms. The averaged defect is
bounded in L1 without multiplying by M. This gives the source's raw tail
bound

    C[eta^(-2)/(M alpha^2)+eta^(-3)epsilon_N,k/alpha
                                +(eta+k^(-1/2))/alpha].

If eta is not small compared with alpha, use the trivial probability bound
1 and enlarge the constant. In the frozen sequence eta/alpha->0, so the
nontrivial argument applies. Raw values and all ties remain observed;
the mollifier never becomes noise, readout data or a new microscopic law.

## 4. Midpoint ties, floors and one final-event transfer

Ideal p_control>p_A7. The raw rule selects control above the midpoint and
A7 otherwise. A wrong decision may satisfy a weak rather than strict
deviation inequality. Apply the tail result at, for example, one quarter
of the ideal gap to cover equality; only the constant changes.

Exact k=floor(N^(63/64)), M=floor(N^(1/4)), eta=N^(-1/32) give the
five powers -5/32,-9/64,-49/128,-1/64,-61/128. Boundary terms are polynomial
times exponentially small under the inherited MT departure bound. k<=N/2
is eventually true, for example at N>=2^64; the exact N in12N restriction
is retained. This is not practical resource guidance.

For the actual iid parent use a common deterministic late origin whose
count-law TV error is epsilon_prep,N=o(1) for both models. Exchangeability
lifts that to labelled TV. Couple future paths to the stationary reflected
chain; its departure probability on M native units is <=CN M pi_D(B_D).
Thus transfer the final bounded decision event by adding

    epsilon_prep,N+CN M pi_D(B_D),

once, without alpha/eta division or a factor M on the preparation error.
This avoids interpreting an undefined reflected semigroup on initial states
outside D; any expanded source initial-exception ledger must include those
states explicitly. There is no survivor postselection or unbounded-moment
transfer through TV. The hidden endpoint at M is not classifier data.

## 5. Independent arithmetic and what the checks do not prove

[ANGULAR_GAP_CERTIFICATE.json](ANGULAR_GAP_CERTIFICATE.json) encloses the
positive ideal slope from inherited rational spectra with mpmath.iv and
exact raw-dyadic outward decimal rounding. For v=1+alpha beta,

    dp/dv=1/(pi sqrt(v)(1+v)).

The mean-value theorem and v<=1+beta_control/16 pay gap>=d_ang alpha/2
on0<alpha<=1/16. This removes the need for a merely numerical slope check;
it supplies no numerical native constants, practical N0 or physical gain.

The source's nominal high-precision spectral and angular values are not
outward certificates. Several source flags called exact gain checks first
convert to float. The intake instead compares actual Fourier sums exactly
in Q(sqrt(3)), after removing their common sqrt(2), and keeps exact ties.
Its labelled finite jump sums use a separately versioned A7 construction.
These small tests validate conventions, not the all-N analytical estimates.
