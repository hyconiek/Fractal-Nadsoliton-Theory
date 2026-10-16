# Rate-coupling audit and two additional finite-model controls

Date:2026-10-16. Analytical mathematics under supplied Gibbs/detailed-balance
premises, with independent exact rational arithmetic. No physical law is added.

## 1. The derivation is correct, including the general diagonal

For a labelled single-slot change a->b put u=e_b-e_a and r=n-e_a.
Detailed balance gives Phi_i=log(q_xy/q_yx)=V(y)-V(x), hence

    Phi_i=h_i(b)-h_i(a)+(g/N)u^T A r+(g/(2N))(A_bb-A_aa).

This uses symmetry of A and keeps the actual parent denominator N. Changing
a distinct slot c->d changes r by v=e_d-e_c. Its mixed difference is

    C_(ab;cd)=(g/N)u^T A v.

Fields fixed independently of the other slots and the self term cancel.
No constant-diagonal shortcut is needed until specializing A7. This is the
rate counterpart of the already admitted RG-E02 mixed Gibbs odds identity.
The all-N proof is this algebra, not the576 original or1116 new fixtures.

At positive g the A7 form is positive semidefinite, rank7 on the
eleven-dimensional zero-sum tangent space; its null directions are the
four omitted k1/k2 Fourier directions. Feature distances are not spatial
distances. The plaquette **cycle** affinity is zero: four potential differences
around a closed loop telescope. C is a mixed derivative of a potential,
not a circulation, Hamiltonian current or violation of detailed balance.

For exchangeable lumpable labelled kinetics, count rates have the factor
n_a. Their reversed edge has factor n_b+1, so subtract
log[n_a/(n_b+1)] at each edge. The repeated-label fixture's multiplicative
bias is exactly1/3, giving -log3. A count description with nonexchangeable
fields need not close at all. Partial-block MM data do not supply the full
state-conditioned rates required here.

## 2. Exact mathematical number, not a sourced constant

For the same N,g>0, the proposed four-rate products obey

    log R_d=(2g/N)(A00-A0d),
    rho36=log R3/log R6=(1+Lambda6/(Lambda3+Lambda5))/2.

The Fourier cancellations are correct: Lambda4 drops out of both gaps.
The numerator and denominator are positive. At g=0 the logarithms both
vanish and the ratio is undefined. Fixed fields cancel even though raw
individual rates can change. Every symmetric edge activity cancels in
each forward/reverse ratio, not only a common clock factor.

The source certificate is genuinely end-to-end rational: kernel cosines
at rational angles have decreasing alternating terms, root endpoints pay
integer inequalities, harmonic cosines lie in Q(sqrt3), and all divisions
have paid nonzero signs. Public decimals are rounded outward, not inferred
from precision. Its row-gap interval-overlap fixture alone would not prove
the Fourier identity, but that identity is independently exact algebra.

Our independent route uses d^(9/5)=d*(d^4)^(1/5), dyadic bisection,
a degree40 Taylor polynomial with absolute Lagrange remainder, and the
full twelve-entry DFT. Every displayed original spectral/gap/ratio interval
contains the independent raw rational box. The admitted enclosure is

    .774903147882614332681671 <= rho36 <= .774903147882614332681672.

This check does not rely on mpmath.iv or libm for the enclosure. It still
assumes the supplied exact-decimal kernel, A7 projection and correctness
of the analytical interval/Taylor argument and Python integer arithmetic;
it is not formal proof-assistant verification. The separate ordinary
floating rate fixtures are numerical diagnostics only.

The source Barker fixture uses prefactor1/12 rather than the frozen E04
g0-normalized prefactor2/12. This does not affect any affinity here. Our
independent full-potential fixture uses2/12. Neither calculation licenses
an absolute-clock comparison across those source fixtures.

## 3. What kinetic information remains entirely free

Given V, every positive symmetric edge function a_xy=a_yx supplies

    q_xy=a_xy exp[(V(y)-V(x))/2].

It is reversible with the same Gibbs law; conversely every positive
reversible rate has this form with a_xy=sqrt(q_xy*q_yx). The source's
c_xy=pi_x q_xy is a different normalization of symmetric activity:
c_xy=a_xy exp[(V(x)+V(y))/2]/Z. Both parameterizations show the same freedom.
Neither is fixed by the entire family of mixed affinities.

The source's distance-three thinning h3(g)=1/(1+g^2) is an exact surviving
countermodel. It is symmetric on edges, respects slot/D12 symmetry,
retains all tested support, can be implemented by rejection at the same
unit proposal clock, and equals native HB at g0. At g3 it changes the
0->3 versus0->6 outgoing branching ratio by exactly1/10. Thus clock and
symmetry conventions do not remove this freedom. This comparison is not
adopted as a physical FIN law.

### Stationary projected rates do not prove Markov closure

The actual persistent AUX-HB-OU family already supplies another important
boundary. Under its stationary joint law, define the labelled conditional
instantaneous rate

    qbar(x,y)=integral rho_n(x)(dtheta) q_theta(x,y).

Joint edge detailed balance, integrated over theta, gives
pi_x qbar(x,y)=pi_y qbar(y,x). Thus qbar obeys the same RC affinity identities,
where the specified edges exist, although the observed labelled/count
history need not be Markov. Reflection restricts the statement to retained
edges and does not change their detailed balance. It does not grant access
to rate conditioning through a block instrument.

This is a finite stationary identity, not a kappa->infinity theorem and not
an iid-start result. The separate four-state exact fixture has positive
compression defect while its conditional averaged rates satisfy DB. RC
does not supersede MM's multi-time memory discrimination or reconstruct
the observed path from one-time conditional rates.

## 4. A new exact isofingerprint countermodel

Let P_r be the real sine/cosine projector
P_r(i,j)=cos(2pi r(i-j)/12)/6, and put

    A_delta=A7+delta(P3-P5), delta=1/100.

The projectors are orthogonal rank2 with diagonal1/6 and annihilate the
constant and parity vectors. This preserves trace, diagonal, rank7 and
the parity eigenvalue/potential, while changing Lambda3 to Lambda3+delta
and Lambda5 to Lambda5-delta. Positivity is paid by the inherited rational
spectral margins. No family-wide dynamic/MT theorem is asserted here.

At separations3 and6, P3(0,d)=P5(0,d), including d0. Hence BOTH gaps
A00-A03 and A00-A06, and therefore both R3,R6 and rho36, are unchanged.
This is stronger than mere closeness: the chosen ideal rate relation is
exactly blind to a nontrivial supplied operator alternative.

At separation1, P3(0,1)=0 and P5(0,1)=-sqrt3/12, so

    (A_delta,00-A_delta,01)-(A7,00-A7,01)=-delta sqrt3/12 !=0.

Another affinity can detect this change. The full mixed-affinity array's
double-centered reconstruction is already covered by RG-E02; do not
advertise that reconstruction as a new campaign. Our new Q(sqrt3) fixture
checks the projector products exactly. The single rho36 contrast separates
the frozen flat-P7 alternative but is not an A7 uniqueness/source theorem.

## 5. What has and has not entered the bridge

This is **a received and admitted mathematical object**, with a useful
finite-N prediction and a backend-independent arithmetic certificate.
It is not missing as a file, nor is the derivation rejected. However it
removes no Gibbs, reversibility, operator-source or kinetic premise: its
own explicit countermodel preserves the missing time-symmetric freedom.
Its source correctly reports K07_admitted=false and KB-002/003 inactive.

The precise remaining kinetic atom is a target-independent restriction on
symmetric activities or their ratios, not another affinity identity.
Computing activities from the already postulated HB law is a consequence
of that law, not an independent source for it. Maximizing entropy, adding
a reservoir/force/incidence graph or selecting a clock is not a deduction
unless the required new premise is explicitly supplied and approved.

The referenced kinetic-composition submission is hash-bound context. This
intake does not accept that entire separate literature/research package.
All ontology, physical source/units/incidence/conservation/EV, selector,
FIN-M0/B4/task337/task335 and quantum/gravity/ToE gates remain unchanged.
