# A coupling that predicts a kinetic rate relation

Date: 2026-10-16. Status: **SCOPED_ANALYTICAL_DERIVATION_AND_RATIONAL_CERTIFICATE**.
This is separately requested intrinsic mathematics. It is not a new physical
source or an activation of KB-002/003.

## 1. Premises and provenance

Use the existing labelled finite-N FIN parent, N>=2, labels in Z12, with
the supplied symmetric matrix A=A7. Let n(x) be the occupation vector and

    V_N(x) = (g/(2N)) n(x)^T A n(x) + sum_i h_i(x_i),
    pi_N(x) = exp(V_N(x))/Z_N.                                  (1)

The optional one-body fields h_i are fixed independently of the other labels.
The main numerical example uses zero fields and g>0. This potential is
dimensionless; no mechanical-energy interpretation is inferred.

Let q(x,y)>0 on the single-label edges used below, and assume detailed balance

    pi_N(x) q(x,y) = pi_N(y) q(y,x).                             (2)

These are inherited model assumptions, not consequences of the nadsoliton
ontology alone. In particular, stationarity without detailed balance would
not suffice. No heat-bath, Metropolis or Barker acceptance function is
required by the theorem. The process and rates are labelled; hidden-memory
outputs need not have these instantaneous Markov rates.

The static mixed-conditional-log-odds identity is already admitted in
`fin_generative_operational_review_20261014/INTAKE_20261014.md`, under
RG-E02 (the source's fixed M=3 case). The work here extends its statement
to reverse-rate contrasts for any reversible kinetic realization, arbitrary
N and an explicit count correction, and derives the displayed A7 relation.
The underlying finite-difference identity and detailed balance are reused.
No original general statistical-mechanics theorem is claimed.

## 2. One edge affinity and its coupling to another label

For a change a->b at slot i, set u=e_b-e_a. Its rate affinity is

    Phi_i(a->b;x_-i) = log[q(x,x+u)/q(x+u,x)]
                     = V_N(x+u)-V_N(x).                        (3)

Writing r=n-e_a, the exact expansion of (1) gives

    Phi_i = h_i(b)-h_i(a)
          +(g/N) u^T A r +(g/(2N))(A_bb-A_aa).                  (4)

For the supplied circulant A the last term vanishes. Keeping it until this
point avoids silently using the constant-diagonal shortcut for other A.
The original denominator remains N, including in the N-1 cavity.

Now change a distinct slot j from c to d; let v=e_d-e_c. Subtraction of
(4) gives the exact coupling

    C_(ab;cd) = Phi_i(a->b | x_j=d,z) - Phi_i(a->b | x_j=c,z)
              = (g/N) u^T A v
              = (g/N)(A_bd-A_ad-A_bc+A_ac).                     (5)

Here z contains all remaining labels. Their values, fixed one-body fields
and the self-diagonal contribution cancel. The result is exact at finite N,
independently of the initial law. No iid start, large-N limit or stationary
sample approximation is used.

Equivalently, with x_ac=(a,c,z), x_bc=(b,c,z), x_ad=(a,d,z), x_bd=(b,d,z),

    exp(C_(ab;cd))
      = q(x_ad,x_bd) q(x_bc,x_ac)
        / [q(x_bd,x_ad) q(x_ac,x_bc)].                           (6)

Equation (6) predicts a product/ratio of **four actual transition rates**.
It is not a fit of a relaxation curve or an assumed auxiliary field.

## 3. Typed coupling form and current interpretation

Define on the zero-sum occupation tangent space T the bilinear form

    B_N : T x T -> R,       B_N(u,v)=(g/N)u^T A7 v.              (7)

Its elementary entries are exactly (5). Since A7=XX^T,

    B_N(u,v)=(g/N)(X^T u).(X^T v).

Thus B_N is symmetric and positive semidefinite for g>0. Its rank is seven;
within the eleven-dimensional T it has the four omitted k=1,2 Fourier
directions as kernel. In particular

    C_(ab;ab)=(g/N)||X_b-X_a||^2.                              (8)

This is an internal feature-space coupling of configuration changes.
It is not spatial incidence, physical distance, a new layer below the
nadsoliton, or a new fundamental field.

The connection with probability current is explicit. For any instantaneous
distribution rho define c_xy=pi_x q_xy=c_yx. Then

    J_xy(rho)=rho_x q_xy-rho_y q_yx
             =c_xy(rho_x/pi_x-rho_y/pi_y).                      (9)

The coupling fixes changes in the rate-affinity force entering this current;
it does not determine the positive conductance c_xy. At rho=pi, all these
currents vanish. In particular, C is **not** a nonzero cycle affinity:
the affinity around x_ac->x_bc->x_bd->x_ad->x_ac telescopes to zero under
(2). A nonzero mixed response must not be called stationary circulation.

## 4. A concrete parameter-free rate relation for A7

Keep the same N and g>0. For d in {3,6}, select two distinct slots and use
the four configurations (0,0,z), (d,0,z), (0,d,z), (d,d,z). Write

    R_d = q((0,d,z),(d,d,z)) q((d,0,z),(0,0,z))
          / [q((d,d,z),(0,d,z)) q((0,0,z),(d,0,z))].              (10)

For A_r=A_(0,r), circulant symmetry gives

    log R_d = (2g/N)(A_0-A_d).                                 (11)

The existing projected A7 has active eigenvalues Lambda_3, Lambda_4,
Lambda_5, Lambda_6, with the conjugate mode pairs retained. Consequently

    A_d = [2 sum_(k=3,4,5) Lambda_k cos(2 pi k d/12)
                         +Lambda_6 (-1)^d]/12,

    A_0-A_3=(Lambda_3+Lambda_5+Lambda_6)/6,
    A_0-A_6=(Lambda_3+Lambda_5)/3.                              (12)

Both gaps are strictly positive. Therefore

    log(R_3)/log(R_6)
       = (Lambda_3+Lambda_5+Lambda_6)/(2(Lambda_3+Lambda_5))
       =: rho_36.                                             (13)

Equivalently, R_3=R_6^rho_36. The common gain and N cancel. The symmetric
edge activities cancel separately by (2); an overall rate/clock factor
cancels as well. This remains true if the activities depend on the full
configuration, provided (2) holds and all specified edges exist.

For the repository's exact-decimal mathematical input

    w_d=cos(0.18575 d+0.1625)/(1+d^(9/5)),  d=1,...,6,

the fresh rational enclosure is

    rho_36 in [0.77490314788261, 0.77490314788262].              (14)

The tighter output in CERTIFICATE.json is derived from the formula, not
assumed as an input. There is no numerical fitting. A7 is nevertheless
supplied: (14) is a dimensionless model prediction, not a derived constant
of nature. At g=0, R_3=R_6=1 and the ratio of logarithms is undefined;
the positive-g condition must not be dropped. Distances 3 and 6 are label
separations, not spatial lengths.

As one existing countermodel, FLAT_P7_TRACE has equal positive eigenvalues
on the same seven retained modes. Equation (13) gives exactly 3/4 in that
case, regardless of its common spectral scale. The two model predictions
are separated by more than 0.02490314788261. This is a theoretical contrast,
not an empirical discrimination result or a natural selection of A7.

## 5. Count rates require a multiplicity correction

When the labelled kinetics is exchangeable over slots and lumpable to n,
let k_(a->b)(n) be the rate for one specified slot and

    Q(n,n-e_a+e_b)=n_a k_(a->b)(n).

At n'=n-e_a+e_b the reverse multiplicity is n_b+1. Therefore the quantity
corresponding to (3) is

    Phi_(a->b)(n)
       =log[Q(n,n')/Q(n',n)]-log[n_a/(n_b+1)].                  (15)

Using raw count-rate ratios in (5) is generally wrong. The checker preserves
a repeated-label negative control: n_0=2,n_3=1, followed by two 0->3 changes.
The uncorrected mixed contrast has an extra term -log(3). Subtracting the
multiplicity at each edge restores (5) exactly. With slot-dependent fields
or nonexchangeable kinetics a closed count generator need not exist; use
the labelled statement instead.

Likewise, a predetermined block readout is not the full labelled process.
The MM observation instrument does not automatically supply the conditional
rates in (10). No apparatus, conditioning access or finite-record error
budget is inferred here.

## 6. An explicit surviving kinetic freedom, including the g=0 clock

Let q_HB be the existing unit-attempt heat-bath generator. For a single-slot
transition whose internal label distance is d, define only as a countermodel

    h_d(g)=1/(1+g^2) if d=3, and 1 otherwise,
    q_tilde(x,y)=h_d(g) q_HB(x,y).                              (16)

Set the diagonal by the row sum. This is positive on the same support and
reversible with the **same** pi, since h has the same value on an edge and
its reverse. It preserves D12 label symmetry and slot permutation symmetry.
It can be realized at the same unit proposal clock by rejecting a proposed
distance-three change with probability 1-h_3(g). At g=0 it is exactly HB.

It preserves every C and R above. But at g=3, from any state with first
label zero,

    [q_tilde(x,x^(0->3))/q_tilde(x,x^(0->6))]
       = (1/10)[q_HB(x,x^(0->3))/q_HB(x,x^(0->6))].             (17)

Both transitions have positive rates. Thus even the same Gibbs law, complete
family of coupling identities, symmetries, unit attempt convention and g=0
control do not fix same-source branching ratios at positive g.
The diagnostic (16) is not adopted as a new FIN law.

This is the exact boundary of the result: we have derived a coupling and
a rate relation from accepted premises. We have **not** supplied the missing
time-symmetric activity or an independently justified physical kinetic law.
The static mechanism behind (5) was already present in RG-E02. The bridge
therefore remains OPEN_INPUT, without running another no-input KB audit.

## 7. Verification and falsification scope

The proof of (5)–(13) is analytical. The standard-library checker adds:

- 576 exact finite-difference fixtures with nonconstant diagonal and fixed
  slot fields, guarding the signs, self term and denominator;
- an end-to-end rational certificate for the mathematical kernel and (14);
- eight binary64 comparisons of independently evaluated labelled rates
  (HB/MET/BAR and (16), two background contexts), kept separate from proof;
- the multiplicity negative control and the changing branching ratio (17).

For the certificate, d^(9/5) is enclosed as the positive fifth root of d^9
using certified integer inequalities. Cosines of rational kernel angles use
successive alternating-series bounds; harmonic cosines are in
{0,+/-1,+/-1/2,+/-sqrt(3)/2}, with a certified rational sqrt(3) enclosure.
All subsequent interval arithmetic and decimal rounding use Fractions.
Fresh active-spectrum boxes lie inside the actual inherited rational boxes.
No mpmath, libm, point minimizer or unbounded tail enters the certified number.

For exact supplied generator rates, violating (5) refutes the conjunction
of the claimed Gibbs interaction, detailed balance and representation.
It need not identify which premise failed. For measurements, exact equality
is not an acceptance test: rate conditioning, simultaneous uncertainty and
the positive-denominator margin must first be supplied. No empirical gate
or sample-size theorem is activated.
