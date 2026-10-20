# New consequence of the frozen hypothesis: near-maximal histories and HB

Date: 2026-10-20. Conditional theorem within the unchanged hypothesis. This
is proof support from the same author, not an independent physical source
or a new postulate. All logarithms and rates below use nats.

## 1. Full finite-history law from the actual entropy deficit

Use the finite K, pi and unique stationary maximizing law Q* from the
source proof/audited dual. For any other stationary history law P with the
**same pi and allowed support**, put

    L_K(P)=h(P)+sum_xy c_xy(P)log K_xy,
    delta=L_K*-L_K(P)>=0.

For primary K=0/1 on allowed entries, delta=h(Q*)-h(P). References are
fixed; different coarse or marked conventions cannot be mixed in delta.

For n>=1 let P_[0:n], Q*_[0:n] be the full laws of n+1 states. Then

    D_n=KL(P_[0:n]||Q*_[0:n])
       =H(pi)+n[L_K*-sum c(P)log K]-H_P(X_0,...,X_n),  (1)
    0<=D_n<=n delta,
    TV(P_[0:n],Q*_[0:n])<=min(1,sqrt(n delta/2)).       (2)

Proof: Q* assigns pi_(x0) product P*_(xj,x(j+1)) to each legal word,
with positive probability. Since

    log P*_xy=log K_xy+log z_x+log z_y-log pi_x,

its expected log density under P uses only the fixed one-state margins
and c(P). The cross entropy is the first two terms in (1). Subtracting
P's block entropy proves (1). Stationarity gives

    H_P(X_0,...,X_n)=H(pi)+sum_(j=1)^n H(X_j|X_0,...,X_(j-1)).

Those conditional entropies decrease to h(P), while their first value is
at most L_K*-sum c(P)log K by the already paid flow optimum. Therefore
every summand in D_n lies between zero and delta, proving (2). Pinsker
follows by coarse-graining to a binary event and the elementary binary
relative-entropy curvature bound; it is not an inference from large KL.

Also D_n/n increases to delta. Thus finite-prefix divergences alone give
**lower**, not upper, estimates of delta. An application of (2) needs a
genuine bound on the entropy-rate deficit; a few finite correlations or
samples do not pay it. In a triangular limit n delta->0 is needed on the
actual growing record. No numerical physical error guarantee is supplied.

This supplies a second proof of exact uniqueness: delta=0 forces every
finite-block divergence to vanish, hence the whole history law equals Q*.

## 2. A stopped readout, with the infinite waiting-time tail retained

Consider the two-slot difference readout D and condition on D_0=d, whose
fixed probability is p_d>0. The initial conditional microscopic law is
the same under P and Q*, since both have the same pi. The KL chain rule
for that initial partition gives D_(n|d)<=n delta/p_d, hence conditional
block TV<=sqrt(n delta/(2p_d)).

Let tau=inf{j>=1:D_j!=d}, and define the eventual first-exit destination
as a cemetery value if tau=infinity. Under Q* its difference process is
Markov with diagonal r_d=Pbar*_(dd)<1. Then, even if P has long memory or
a positive never-exit mass,

    TV(Law_P(D_tau|D_0=d),J*_(d,.))
      <=min(1,sqrt(n delta/(2p_d))+r_d^n).              (3)

Proof: maximally couple the conditional microscopic blocks through n, then
extend using their original conditional future laws. Conditional futures
exist on the finite-alphabet path space; this is an analysis coupling, not
a controller. If blocks agree and Q* has exited by n, both first destinations
are already equal. The failure probability is at most the block mismatch
plus Pr_Q*(tau>n)=r_d^n. This includes all unobserved null steps and every
source nonexit; no unbounded stopping time is passed through TV for free.
The cemetery outcome is treated as an output, not discarded.

## 3. The HB incompatibility persists away from the exact maximum

Fix a source state d and destination j for which

    Delta=|J_HB(d,j)-J*(d,j)|>0.

Every stationary allowed-history law P with the same pi and that HB
first-exit probability must satisfy, for each integer n with r_d^n<Delta,

    delta >= (2p_d/n)(Delta-r_d^n)^2 >0.               (4)

This follows by applying (3) to the single destination event and rearranging.
It requires no Markov assumption on P and no choice of its holding clock.
For any fixed nonuniform positive p the source's odds proof gives some
such d,j; r_d<1 gives a finite admissible n. Consequently no sequence of
stationary same-pi laws with delta->0 can have exactly the HB first-exit
law for that fixed p. This is a robustness strengthening, not a new model.

For the existing exact q4 fixture,

    p_2=39/110, J*(2,0)=1/5, J_HB(2,0)=15/71,
    Delta=4/355, r_2=3/13.

Freeze n=4 (no optimization or sweep). Then

    delta >= (39/220)[4/355-81/28561]^2
          =21925107363/1739734008083500 >1/100000.       (5)

The unit is **nats per intrinsic labelled step**. This constant belongs to
the rational q4 fixture, not actual A7 numerical probabilities, a physical
rate or a uniform bound in gain/N. Near the uniform-p boundary the HB
branching discrepancy can vanish; no uniform-away-from-zero gap is claimed.

## 4. What is not rescued or promoted

No hypothetical source was adopted as native FIN, and no self multiplicity
was changed. The new estimate does not justify maximizing entropy as a law
of nature. It shows precisely how that same authorized hypothesis responds
to approximate maximization and why its adverse HB comparison is not just
an exact-equality artefact. Source pi, support, reference and preparation
remain premises; strict F01/F03/F06 and all physical gates remain open.
