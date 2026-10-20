# RS: a reversible-source completion of the causal-transfer interface

Date: 2026-10-20. Self-audited conditional mathematics, not a new FIN law,
independent scientific review or a natural source. Projection, reversible
compression and coupling are reused methods. HM already paid the generator
compression identity; CF already paid causal intensity/path transfer. The
new local interface below turns a finite-lag compression defect into a
whole **sampled-record** bound, and states exactly when it implies more.

## 1. Typed inputs and the finite-lag Gram identity (RS-001)

Supply a finite stationary Markov source X, positive stationary weights mu,
a deterministic surjective readout R:S->E, and detailed balance of the
**source**, not only reversal symmetry of its observed history. Let
pi_i=mu(R=i)>0, P=E_mu[.|R], and U_t=exp(tL) in the continuous-time case.
L acts on column functions; transition matrices have row probabilities.
All statements also apply at a fixed step to a supplied reversible
discrete Markov kernel U. No dynamics is constructed from a target here.

The observed stationary pair kernel is

    C_t(i,j)=Pr(R(X_t)=j | R(X_0)=i).

It is reversible for pi but need not be a semigroup. In L2(mu), on the
range of P, self-adjointness gives the exact identity

    P U_(2t) P - (P U_t P)^2
        = P U_t (I-P) U_t P = A_t^* A_t,
    A_t=(I-P)U_t P.                                      (1)

Write f_j=1_{R=j}/sqrt(pi_j), an orthonormal basis of that range. Taking
its trace yields a nonnegative, observable-pair quantity

    delta_t=tr(C_(2t)-C_t^2)
      =sum_x mu_x sum_j
          [U_t 1_{R=j}(x)-C_t(R(x),j)]^2/pi_j >=0.        (2)

The row-coordinate matrix in (1) is similar via diag(sqrt(pi)) to a
symmetric positive-semidefinite matrix; it need not be symmetric in the
ordinary Euclidean metric. The trace in (2) is nevertheless the ordinary
matrix trace. No extra factor pi_i or pi_j is inserted into that trace.
For a discrete kernel replace U_t,U_(2t) by U,U^2.

## 2. Entire sampled record, not just its pairs (RS-002)

Fix t>0 and M>=1. Let P_M be the joint law of
`(R(X_0),R(X_t),...,R(X_(Mt)))`. Let Q_M be the stationary discrete Markov
chain with initial pi and transition C_t. Then

    TV(P_M,Q_M) <= min(1, M sqrt(delta_t)/2).              (3)

Proof. Put a_j(x)=U_t 1_{R=j}(x). At step n, the actual conditional
distribution given the entire observed prefix is
`alpha_j=E[a_j(X_(nt)) | observed prefix]`, by the source Markov property.
Its present row C_t(R(X_(nt)),j) is measurable in that prefix. Conditional
Jensen followed by weighted Cauchy--Schwarz gives

    E sum_j |alpha_j-C_t(R(X_(nt)),j)| <= sqrt(delta_t).

Indeed sum_j pi_j=1, and (2) pays the squared norm after averaging over the
stationary source. Sequential maximal coupling, or telescoping the joint
prefix probabilities, bounds complete-record TV by the sum of the expected
conditional-row TVs. Row TV has the factor 1/2, proving (3). Samples are
not assumed independent and R(X) is not assumed Markov.

If every C_t(i,j)>0, set K_t=max_(i,j) pi_j/C_t(i,j). The chain rule for
finite-alphabet likelihoods and `KL(a||q)<=sum_j(a_j-q_j)^2/q_j` also give

    KL(P_M||Q_M) <= M K_t delta_t,
    TV(P_M,Q_M) <= sqrt(M K_t delta_t/2).                 (4)

This can be better than (3) when K_t is controlled. With an arbitrary
source preparation mu', replacing the stationary record by its actual
record costs at most TV(mu',mu), by source-path contraction. Thus add that
term to (3), when the comparison still starts at pi. No small preparation
error is divided by delta or multiplied by M. No claim about arbitrary
preparations being close is made.

## 3. Exact zero at one lag for a reversible finite CTMC (RS-003)

In continuous time, delta_t=0 at **one** t>0 implies A_t=0. The range of P
is then invariant under U_t. The finitely many distinct eigenvalues of L
are real, and their exponentials at t are distinct and positive. Polynomial
interpolation supplies a real polynomial g with
`g(exp(t lambda))=lambda` at every source eigenvalue. Therefore L=g(U_t),
so the range of P is also invariant under L. Consequently

    sum_(y:R(y)=j) L_xy = qbar_ij
          independently of x with R(x)=i, j!=i.         (5)

These off-diagonal constants are nonnegative and their conservative
diagonal defines a generator. The source short-step equation or CF's
causal criterion now gives a Markov readout with generator qbar for
**every** source preparation. This is strong lumpability, derived here
from stationary pair data plus the source-reversibility hypothesis.

In particular, if C_t=exp(tQ) and C_(2t)=exp(2tQ) for a specified
pi-reversible target Q, the exact readout generator is Q: the same real
spectral-log interpolation makes the reversible logarithm unique. A
nonreversible target logarithm is not covered by this uniqueness claim.
For a reversible discrete source, zero delta only proves invariance under
its one-step kernel, hence Markovity at that supplied step; no continuous
generator is inferred.

This sharpens the scope of ES, rather than invalidating it. Its observed
process is reversal symmetric, but its directed-phase analysis lift is
not a detailed-balanced source. In fact no finite positive-stationary
reversible CTMC realization of that same non-Markov observed law can exist:
its exact pair kernels would make delta_t=0, contradicting (5). This is
an excluded **representation class**, not an independently derived physical
reversibility principle or a proof against non-Markov/field sources.

## 4. Infinitesimal defect and the already paid full-path bound

For a supplied reversible finite L, (2) gives

    V=lim_(t down 0) delta_t/t^2
      =sum_j ||(I-P)L f_j||_(L2(mu))^2.                  (6)

This is the HM generator Gram defect traced over the full readout basis,
not a newly discovered projection identity. Define qbar by conditional
stationary rate averages, and a_ij(x)=sum_(y:R(y)=j)L_xy for j!=R(x).
Weighted Cauchy--Schwarz and (6) imply

    E_mu sum_(j!=R(X)) |a_Rj(X)-qbar_Rj| <= sqrt(V).

CF-002, with the actual stationary preparation, therefore pays

    TV(observed complete paths on [0,T], CTMC(qbar,pi))
        <= min(1,T sqrt(V)).                             (7)

This uses the **infinitesimal** limit V, not a small defect at one coarse
lag. Rates, preparation, horizon and source are inputs. In a triangular
family the relevant integrated bound and preparation must be paid; (7)
does not supply those hypotheses for a native FIN family.

## 5. Essential countercontrols (RS-004)

### Source reversibility cannot be replaced by observed reversal

The ES r=3 block construction has exact pair kernels of a reset Markov
chain at every lag, but a nonzero triple contrast. Its finite directed-phase
lift is not reversible. The independent rational fixture here also works
in discrete time: p=(1/3,2/3), h=(2/3,-1/3), epsilon=1/2. The 24-state
word/phase lift has C_1=C_2=1 p^T and delta=0, yet the observed triple law
differs. This directly rejects (3) with the source hypothesis removed.
Continuous-time ES already supplies the corresponding full-clock example.

### A small coarse-lag defect does not control full continuous paths

Use the reversible three-state source

    L0=[[-1,1,0],[1,-2,1],[0,1,-1]], mu=(1/3,1/3,1/3),
    R(0)=R(1)=A, R(2)=B, pi=(2/3,1/3).

Its stationary mean generator has qbar_AB=1/2, qbar_BA=1. The source
eigenvalues are 0,-1,-3. Its observed nonconstant pair eigenvalue and
doubling defect are exactly

    r(t)=(3/4)exp(-t)+(1/4)exp(-3t),
    delta_t=(3/16)(exp(-t)-exp(-3t))^2.                  (8)

The Markov summary has nonconstant eigenvalue exp(-3t/2). For t>0,
strict convexity of the exponential gives r(t)>exp(-3t/2). The two-time
**joint** laws differ in TV by `(4/9)[r(t)-exp(-3t/2)]`.

Now L_k=k L0 and qbar_k=k qbar. At any fixed sampling lag Delta>0,
delta_Delta tends to zero as k->infinity. But at time t=(log 2)/k the
joint-law TV is the fixed positive number

    (4/9)(13/32-sqrt(2)/4) > 7/300.                    (9)

The rational bound follows from sqrt(2)<283/200. For any fixed T>0 that
time lies in [0,T] eventually, so complete path TV is at least this number.
Here V_k=(3/4)k^2, consistently with (6). This counterfamily has unbounded
rates; it is a general interface countercontrol, **not** an admitted q12
KR rate-capped counterfamily or a refutation of native IBA. It rules out
an unrestricted continuous-path conclusion from the fixed-lag statistic.

## 6. Exact remaining source obligation

No actual independent `(L,R,mu)` has been derived by imposing reversibility
or by calculating C from a selected L. These are conditional compatibility
tests for a provider once available. A target-designed reversible lift can
pass all tests and still be circular as an explanation of the target.
F01/F03/F06 stay OPEN_SOURCE; CSP-001--003 are not activated by this tool.
The physical clock, carrier and readout remain separate source premises.
