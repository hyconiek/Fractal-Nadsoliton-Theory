# Source probe: complete finite-order correlations do not supply Markov closure

Date: 2026-10-20. Self-contained conditional mathematics. This is a
diagnostic family, **not a newly adopted FIN microscopic law or source**.

## 1. What is being tested

The existing source-bound N2 difference process has positive stationary
probability p and supplied generator

    B=lambda(1 p^T-I), lambda=2,
    P_t=exp(-lambda t)I+(1-exp(-lambda t))1p^T.          (1)

We do not derive (1) from ontology. We test whether even independently
providing all two-time laws of (1), its state-conditioned infinitesimal
increment covariance and stationary symmetries would force the same
**whole history law without assuming Markovity**. The answer is no.

The elementary pairwise-independent block idea already appears as a binary
control in the ID intake. The new local object here is its arbitrary-p,
continuous-clock native-reset specialization, full time-marginal formulas
and precise implication for the current mean-bracket provider gate.
The general phenomenon is known in the literature; no priority is claimed.

## 2. A positive block family with arbitrary prescribed p

Let E be a finite set, p_i>0, sum p_i=1. Choose a bounded real h with

    sum_i p_i h_i=0, |h_i|<=1, v=sum_i p_i h_i^2>0.

For each integer r>=3 and parameter `0<|epsilon|<1`, set

    T_r(i_0,...,i_(r-1))
      =prod_(j=0)^(r-1) p_(i_j)
                    [1+epsilon prod_(j=0)^(r-1) h_(i_j)]. (2)

This is strictly positive, normalized and invariant under any permutation
of its r arguments. Summing out **any** one coordinate kills the correction
because E_p h=0. Therefore every proper marginal of a block is the iid p
law, while its r-fold h moment is `epsilon v^r`.

Take independent, identically distributed blocks `(W_k)_(k in Z)` with law
T_r and concatenate them into Z. Independently choose U uniform on
{0,...,r-1} and put `Y_n=Z_(n+U)`. This defines a bi-infinite sequence law.
A one-step shift cycles U and, at its wrap, shifts the iid block index.
Hence Y is stationary. It is also invariant under reversal of integer time:
reversal reverses iid block order, reverses each exchangeable block, and
permutes the uniformly mixed offsets. These are statements about the
constructed law, not a sourced physical time arrow.

Every collection of at most r-1 **distinct** Y indices has joint law equal
to the corresponding iid p collection: in each contributing block it sees
a proper marginal. A collection of r distinct indices can see a full block
only when the indices are consecutive; in that case alignment occurs for
exactly one offset out of r. Consequently

    E prod_(j=0)^(r-1) h(Y_j)=epsilon v^r/r.           (3)

For r distinct nonconsecutive indices the moment is zero. Repeated indices
reduce the number of distinct variables and agree with the iid mark law.
In particular r=3 gives the explicit prediction failure

    Pr(Y_2=k | Y_0=i,Y_1=j)
       =p_k[1+(epsilon/3)h_i h_j h_k],
    Pr(Y_2=k | Y_1=j)=p_k.                            (4)

Y consists of **refresh marks**, potentially with repeated labels, not the
actual-change jump chain in AH. No access to these marks is added to the
observer. The next construction has a contrast using unmarked state
snapshots only.

## 3. Continuous time and the complete two-time equality

Independently use a two-sided homogeneous Poisson point process with rate
lambda>0. Let N(0)=0; N(t) counts points in (0,t] for t>0 and equals minus
the count in (t,0] for t<0. Define `X_t=Y_(N(t))`. Paths are cadlag with
locally finitely many possible changes; null refreshes are unobserved
construction marks. The Poisson clock is a **supplied comparison resource**,
not independently generated physics.

The entire X law is stationary. Indeed, condition on the Poisson counts;
under a real time shift the common integer offset N(s) can be removed using
integer stationarity of Y. The remaining count differences have the same
joint laws by stationarity of Poisson increments. This proof does not
assume independence of N(s) from overlapping shifted increments.
Reversal of the stationary Poisson intervals, together with reversal and
shift invariance of Y, also preserves the finite-dimensional X history law
(path reversal uses the usual left-limit/right-continuous convention).

For t>0, either no refresh occurs, probability exp(-lambda t), or the
indices 0 and N(t) are distinct. Distinct mark pairs are independent p.
Thus **every epsilon and every r>=3 has the exact pair kernel**

    Pr(X_t=j | X_0=i)
       =exp(-lambda t)delta_ij+(1-exp(-lambda t))p_j.   (5)

This equality is for all i,j and all positive lags, not only a return trace
or selected Fourier features. In particular the pair kernels satisfy the
Chapman--Kolmogorov semigroup identity exactly. That identity for pair
conditionals alone is not the Markov property of the full process.

More generally all finite-time distributions at **at most r-1 observation
times** agree with the native reset process: condition on the count vector.
There are at most r-1 distinct mark indices; the joint mark law is iid p
on those distinct indices, with exactly the same repetitions as for iid
refresh marks. This proves agreement, not merely convergence, for all lags.

The comparison need not be dismissed as a nonergodic artefact. It has a
finite Markov **analysis lift** `H=(current block, phase)`, stationary law
`mu(w,j)=T_r(w)/r`: a refresh advances the phase, and at wrap it draws a
new independent block. This preserves mu by direct inflow balance. After
any n>=r refreshes the current block has law T_r independent of the old
block; its phase is the old phase plus n modulo r. At any fixed t0>0,
let `c=min_(r<=n<=2r-1) Pr(Poisson(lambda t0)=n)>0`. Every target phase has
one such n, so the lift kernel satisfies `P_(t0)(h,.)>=r c mu(.)`.
Since `0<r c<1`, repeated minorization gives geometric total-variation
mixing of the lift. The stationary observed X inherits mixing: condition
a bounded future event on the intervening hidden state and apply that
uniform bound. Thus even stationarity, mixing and time-reversal symmetry
of the observed history do not remove this ambiguity. This lift is not a
physically sourced mediator or a new informational layer.

## 4. A bounded higher-time contrast on the unmarked record

Take r observation times with positive gaps t_1,...,t_(r-1), total T.
For the difference of their h-product moments from the iid-mark reset
process, repeated indices contribute zero difference. Distinct nonconsecutive
indices also contribute zero. The only correction occurs when every
Poisson increment across a gap equals one: then r consecutive mark indices
form a full block with probability 1/r. Independent Poisson increments give

    Delta E[prod_(j=0)^(r-1) h(X_(t_1+...+t_j))]
      =(epsilon/r) v^r lambda^(r-1)
                            prod_(j=1)^(r-1)t_j exp(-lambda T). (6)

It is nonzero. The comparison process is therefore not the Markov reset
law in (1), despite (5). A homogeneous Markov law with positive stationary
p and those pair kernels has its r-time cylinders fixed by their product;
(6) contradicts that factorization. This proves non-Markovity of the
observable X, not a failure of the native FIN model.

For r=3 and t_1=t_2=t,

    Delta E[h(X_0)h(X_t)h(X_(2t))]
       =(epsilon/3)v^3 (lambda t)^2 exp(-2lambda t).    (7)

The observable is bounded by one and uses no attempt marks. Hence its
nonzero difference also bounds the total variation of these snapshot
triples from below by half its absolute value. No apparatus or sample
complexity theorem is inferred.

For the rational control p=(1/2,1/3,1/6), h=(1/2,-1/2,-1/2), epsilon=1/2,
one has v=1/4. At lambda=2 and t=1/2, (7) equals `1/(384 e^2)>1/3456`,
using the elementary e<3. This is a small exact mathematical control,
not a small-N FIN physical witness or an observed number.

## 5. The same state-conditioned covariance, not the same predictable law

For one-hot state observations, (5) directly implies

    lim_(dt down to0) Pr(X_(t+dt)=j | X_t=i)/dt=lambda p_j (j!=i),
    Gamma_i=sum_(j!=i) lambda p_j(e_j-e_i)(e_j-e_i)^T,
    G=sum_i p_i Gamma_i=2lambda[diag(p)-p p^T].          (8)

These are exactly the Markov reset values, including the full matrix,
every conditional entry and the mean state-change intensity
`lambda(1-sum p_i^2)`. But (8) conditions only on the present label. It
does **not** supply a predictable intensity or martingale bracket depending
only on X_t in the full observed-history filtration. Hidden mark correlations
make that stronger conditional-law statement fail. This is the missing
premise, not a contradiction of (8).

The existing formula `Q=-(1/2)diag(p)^(-1)G` still correctly identifies the
unique **reversible Markov candidate**. It cannot certify that an arbitrary
history source with the same G is that Markov process. The original
mean-bracket inverse and AH theorems explicitly assume their Markov class
and remain valid. This research tests an attempted extension outside that
class, not their already closed inverse targets.

## 6. Binding to the actual N2 target and symmetries

For supplied finite nonzero gain in the accepted N2 model,
`p_j=softmax_j[(g/2)A7_(0j)]`, lambda=2. The admitted centered circulant
rank-seven A7 is not constant on its row, so p is not uniform. One can choose
without selecting an orientation or special label

    h_i=p_i-c, c=sum_j p_j^2.

Then E_p h=0, |h_i|<1, and v>0. Every state permutation preserving p also
preserves h, so (2) retains those stationary-state symmetries. Uniform p
is still covered by the general theorem with an explicit nonzero contrast
h, but that choice need not preserve all state symmetries; no stronger g=0
symmetry claim is made here. No actual 144-state labelled-parent lift is
asserted. The exact countercomparison is on the accepted twelve-state
**difference quotient's signatures**, not the full native labelled law.

## 7. Consequence for F01

For any fixed finite correlation order k, choose r=k+1 (r>=3). Then every
full k-time observation law of the native reset target can be preserved
while its higher history law changes. Thus adding a few more covariance
entries, correlation orders or pair-semigroup tests cannot independently
pay source closure outside a justified Markov/state-sufficiency class.

This does not prove that every possible FIN source is impossible or that
all intrinsic evolution is static. It shows that the candidate provider
must export more than a finite hierarchy of signatures **or must derive
the relevant screening/state-completeness theorem**. Imposing Markovity,
maximum path entropy or a Poisson reservoir by fiat would be an additional
premise, not a source found in the currently inspected FIN equations.
No such independent source theorem is proved here: F01 remains OPEN_SOURCE.

## 8. A positive, exact causal-source criterion

The missing condition can be formulated more precisely than "supply more
correlations". Let X be a finite-state cadlag process and F_t its actual
observed-history filtration. For a bounded conservative matrix B, suppose
an independent source derivation pays, for every state j, the martingale

    M_j(t)=1_{X_t=j}-1_{X_0=j}-integral_0^t B_(X_u,j)du. (9)

Then X really has the homogeneous Markov transition law exp(tB), not only
its pair marginals. To prove this, fix s<t and a terminal function f and set
`u(r,i)=(exp((t-r)B)f)_i`. It satisfies `partial_r u+B u=0`.
Multiply (9) by the deterministic C1 coefficients u(r,j), sum over j and
use integration by parts. The drift cancels, so u(r,X_r), s<=r<=t, is a
martingale. All coefficients are bounded on this finite interval. Thus

    E[f(X_t)|F_s]=(exp((t-s)B)f)(X_s).                 (10)

Iteration gives every finite-time cylinder probability, proving the claimed
whole-law identification. This proof imports no unverified uniqueness
theorem. For the counterfamily above, (9) with the native B **cannot** hold
in F_t, because it would contradict the nonzero contrast (6).

Equation (9) is a sufficient **causal source interface**, not an independent
provider found in FIN. Writing B into (9) by fiat is circular. A genuine
provider must derive its full-history drift/compensator from its own admitted
evolution and state/readout map, independently of the desired native B.
Markov closure alone still does not select among the old reversible
activity alternatives. Both causal closure and the remaining kinetic
restriction must be paid at their appropriate levels.
