# One-bit critical information: conditional entropy, native transfer and source floor

Date: 2026-10-21. Conditional analytical theorem, not independent review.

## BI-001. The bit is not assumed Markov

Let X_j be a stationary Markov chain (here the full reflected count skeleton)
and Y_j=f(X_j) a binary readout. Conditioning on the full previous state
can only reduce entropy; Markovity of X, not of Y, gives

    H(Y_j|Y_0,...,Y_(j-1))>=H(Y_j|X_(j-1))=c_N(t).

Consequently, for M transitions and M+1 bits,

    H(Y_0,...,Y_M)>=M c_N(t),
    entropy_rate_per_sample(Y)>=c_N(t).                   (1)

The bound may be strict. Replacing the bit by a fitted two-state transition
matrix or calling c_N its exact entropy rate would be incorrect.

### Pay the discontinuous threshold in NT's actual core

Approximate f=1{z>0} by continuous bounded g_eta in[0,1], equal to f outside
[-eta,eta]. MT static convergence and rho({0})=0 give

    limsup_N ||f_N-J_N g_eta||_2^2<=rho([-eta,eta])->0.

Native/limit contraction and the admitted strong C_b core compare
S_N(t)g_eta with J_N S(t)g_eta. The latter is bounded continuous by the
paid Feller argument. Binary entropy is uniformly continuous on[0,1];
equivalently its binary-coupling bound controls its expectation under L1
approximation. Static convergence therefore gives, after eta->0,

    c_N(t)->c(t), for every fixed t>0.                    (2)

Concavity of b2 and stationarity also give
`c(s+t)=E b2(S(s)S(t)f)>=E S(s)b2(S(t)f)=c(t)`, and similarly for c_N.
The positive small-time coefficient proved below therefore implies
c(t)>0 for EVERY fixed t>0, not only for a shrinking sequence of lags.

This never evaluates an arbitrary L2 representative on native atoms and
never imports convergence of posterior filters or general information.

The same conclusion holds if t_N->t>0. For a nonnegative selfadjoint
generator, spectral calculus and lambda exp(-r lambda)<=1/r give

    ||S_N(t)-S_N(s)||_(2->2)<=|t-s|/min(t,s), s,t>0.       (3)

This is uniform away from zero, independently of the large native fast
eigenvalues. It pays later coarse-lag rounding, NOT arbitrary t_N->0.

## BI-002. A paid short-time entropy boundary layer

Write Phi for the standard normal CDF and

    J_bit=integral_R b2(Phi(u))du,
    C_bit=sqrt(2) rho(0) J_bit,
    rho(0)=2/[12^(1/4) Gamma(1/4)].

Then

    c(t)=C_bit sqrt(t)+o(sqrt(t)), t down to0.             (4)

Here is a direct payment without quoting an unread heat-kernel theorem.
Couple solutions from x<y by the same Brownian path. Their difference is

    (y-x) exp[-(1/3) integral_0^t (Z_y^2+Z_y Z_x+Z_x^2)ds]>0.

Nonexplosion is already paid by NT's quadratic Lyapunov bound. Thus
P_x(Z_t<=0) decreases in x. Stop at |Z|=1, where |b|<=1/3. For
0<=x<=1/2 and t<=3/4, exit requires a Brownian maximum at least
1/(4sqrt(2)); reflection and the union bound give

    P(exit by t)<=4 Phi[-1/(4sqrt(2t))],
    P_x(Z_t<=0)<=Phi[(t/3-x)/sqrt(2t)]+P(exit by t).       (5)

For x>=1/2, monotone coupling reduces the same probability to its value
at1/2, hence bounds it by C exp(-c/t) uniformly even for unbounded x.
Reflection symmetry of the cubic SDE supplies the corresponding negative
half-line bound; no global bounded-drift assertion is used.

At x=sqrt(2t)u with u in a fixed compact set, the stopped drift integral
has size <=t/3. Gaussian sandwich bounds and (5) yield

    P_x(Z_t>0)->Phi(u)

uniformly on that compact set. For tail control use

    b2(p)<=2 sqrt(min(p,1-p)),
    Phi(-u)<=exp(-u^2/2)/2, u>=0.

The first follows from log v<=v-1 applied to sqrt(1/p) and1/(1-p),
then symmetry; the second follows by shifting the Gaussian tail integral.
Thus in the rescaled near region the entropy is bounded by a Gaussian
function of u plus C exp(-c/t). The latter is integrated over an interval
of length O(t^(-1/2)), so it still vanishes. In the far region its rho-
weighted integral is exponentially small. Since rho is bounded by rho(0),
local convergence plus these tails proves (4) with dx=sqrt(2t)du. This
pays the growing-domain remainder rather than asserting a global majorant.

The constant is positive and finite, with the simple analytical bounds

    7/16 < J_bit <=2 sqrt(2 pi).                          (6)

For the lower bound, Phi(1)<7/8 follows from pi>3 and
integral_0^1 exp(-v^2/2)dv<=103/120, since
(103/120)^2/6<(3/8)^2. Hence on[-1,1], p=Phi(u) lies in[1/8,7/8], and
b2(p)>=2p(1-p)>=7/32. The upper bound follows from the two tail bounds
above, integrated over R. No numerical quadrature or useful small-t
threshold is supplied or needed.

For scale context only, the sign-flip probability satisfies

    P(Y_t!=Y_0)~2 rho(0)sqrt(t/pi).

The same local/tail argument pays this, using integral_0^infinity Phi(-u)du
=1/sqrt(2pi). Therefore the observed stationary bit entropy rate per unit
time lies between order t^(-1/2) and order t^(-1/2)log(1/t). We do NOT
claim its exact coefficient or an optimal achievable entropy cost. The
asymptotic coefficient in (4) is for the conditional lower bound c(t).

## BI-003. Source constraint: reuse decimation, not ID's N2 calculation

The [already admitted Fano/power/decimation proof](../fin_intrinsic_entropy_source_review_20261019/SUBSAMPLING_ERROR_FLOOR.md)
applies to any supplied stationary source/readout joining. For a source
with entropy h(T_1)<=H in this CRITICAL comparison clock, sampled at1/m,
integer decimation to lag k/m retains its mismatch e. Conditional chain
rule and binary Fano give

    b2(e)>=c_N(k/m)-(k/m)H.                              (7)

This uses the full stationary history joining, not two marginals or one
finite window. The source readout need not be Markov, and no consistency
between the joinings at different m is required. The rational power rule
is exactly the one already paid in ID; it is not an arbitrary-clock axiom.

Let N->infinity and m_N->infinity. For any fixed t>0, take
k_N=max(1,floor(m_N t)). Then k_N/m_N->t>0, so (2)–(3), NOT a shrinking-
lag native CLT, imply

    liminf_N b2(e_N)>=c(t)-tH.

Define the explicit necessary-bound functions

    D_bit(H)=sup_(t>0)[c(t)-tH]_+,
    e_bit(H)=b2^(-1)(D_bit(H)), inverse on[0,1/2].

Every such source family with one uniform finite H obeys

    liminf_N e_N>=e_bit(H)>0.                            (8)

This also applies directly to the limiting stationary diffusion bit. For
H=0, S(t)f->1/2 in L2(rho): NT's form kernel contains only constants,
so the spectral theorem gives this projection limit. Thus c(t)->log2,
D_bit(0)=log2 and e_bit(0)=1/2. For H>0 positivity follows from (4).

As H grows, the bound function has the asymptotic

    D_bit(H)~C_bit^2/(4H),
    e_bit(H)~C_bit^2/[4H log H].                          (9)

To check the first, use c(t)=(C_bit+o(1))sqrt(t), maximize the two bounding
quadratics in sqrt(t), and exclude t bounded away from0 by c(t)<=log2.
The second uses b2^(-1)(v)~v/log(1/v). This characterizes OUR NECESSARY
LOWER BOUND, not an attained rate-distortion optimum or a constructed
source. No H-independent obstruction is claimed as H becomes unbounded.
The order is: the native limit at fixed H, then the large-H behavior of
this limiting bound. A provider's growing H_N under critical-clock
acceleration is not silently replaced by a constant or covered by a
new uniform-in-H_N native estimate.

## BI-004. Actual iid finite records need their own initial budget

The full stationary-source theorem above concerns the stationary reflected
comparison, not an infinite stationary law of the actual iid parent.
For the actual finite records, use the separate
[finite-record proof](FINITE_RECORD_AND_INITIAL_BUDGET.md). It explicitly
pays native TV-to-block-entropy continuity and includes a source prefix
budget I_N. Finite KS alone cannot bound a supplied finite record.

No source mechanism, apparatus, clock or physical entropy law is adopted.
Positive information entropy coexists with the separately proved vanishing
reversal-KL arrow. Neither fact closes the natural F01 source gate.
