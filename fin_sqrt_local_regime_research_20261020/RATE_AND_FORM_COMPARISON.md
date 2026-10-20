# SR-001: exact rates, Gibbs weights and uniform forms

Date: 2026-10-20. Conditional on the objects in [CONTRACT.md](CONTRACT.md).

## 1. Labelled law and actual count quotient

On E^N a single labelled slot a->b at cavity r changes the all-pairs
exponent by h_b-h_a. Indeed, for d=e_b-e_a,

    V(n+d)-V(n)
      =(g/N)d^T A(n-e_a)+(g/(2N))(A_bb-A_aa)
      =h_b-h_a.

The second term vanishes because A has constant diagonal; it is not a
discarded state-dependent self interaction. Labelled rates exp((h_b-h_a)/2)
obey detailed balance with exp(V). Slot-permutation covariance gives the
count rate Q(n,n')=n_a exp((h_b-h_a)/2). Its reverse has departure n_b+1
and the same cavity. Hence

    W(n')/W(n)=[n_a/(n_b+1)]exp(h_b-h_a),
    W(n)Q(n,n')=W(n')Q(n',n).

This pays the invariant unnormalized measure W on the entire finite count
space, not merely its local approximation. Positive rates give exactly
the same legal graph and reflected connected components as HB.

## 2. A global rate comparison, not a shared-Gibbs inference

The cavity has nonnegative entries summing to N-1. For |g_N|<=G,
|h_j|<=H=G M_A. Consequently

    Q(n,n')<=n_a exp(H),
    sum_(n'!=n)Q(n,n')<=RN, R=(q-1)exp(H).

Let Z_h=sum_j exp(h_j). Native HB on the same edge is
B(n,n')=n_a exp(h_b)/Z_h, with its original clock. Therefore exactly

    Q(n,n')/B(n,n')=Z_h exp(-(h_a+h_b)/2),
    a_*<=Q/B<=a^*,
    a_*=q exp(-2H), a^*=q exp(2H).                    (1)

All these constants are independent of N and the local state. Moves with
n_a=0 are absent from both laws, not divided by zero. The comparison is
valid globally and survives deletion of outgoing edges for reflection.
At g=0 the ratio is q, explicitly retaining the clock difference. At
nonzero g it need not be one common clock multiplier.

## 3. Dirichlet forms and gap, with the actual invariant measure

For a reflected reversible count generator L on D and real f,

    E_L(f,f)=(1/2)sum_(x,y in D) pi_D(x)L(x,y)(f(y)-f(x))^2.

Apply (1) edge by edge, never treating Gibbs weights as uniform:

    a_* E_B<=E_Q<=a^* E_B.

Taking the infimum over the same nonconstant f/Var_pi(f) gives

    a_* gap(B_D)<=gap(Q_D)<=a^* gap(B_D).             (2)

MT's paid nonuniform root-path proof supplies
gap(B_D)>=gamma_*(N+1)^(-K), with fixed finite, uncalibrated gamma_*,K.
Thus Q has a polynomial reflected gap in its declared clock. No new
canonical paths or optimized exponent are needed. This is an application
of the inherited theorem, not a newly discovered metastability method.

The independent parity Rayleigh upper obstruction also persists. Each
count jump changes m by at most 2/N. Thus E_Q(m,m)<=2R/N, whereas MT's
static law gives Var_pi_D(m)>=v/sqrt(N) eventually, for some v>0. Hence
gap(Q_D)<=C/sqrt(N), consistent with (2). This does not identify a critical
SDE or turn ten stable directions into a constant full-chain gap.
