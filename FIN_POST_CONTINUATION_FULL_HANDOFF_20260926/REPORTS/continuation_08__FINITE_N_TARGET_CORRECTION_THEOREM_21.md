# FINITE-N-TARGET-CORRECTION-THEOREM-21 — exact row-sum classification

Status: **PROVED WITHIN DECLARED HEAT-BATH + ME7 LANE; FULL-TENSOR BASIS REPLAY**.

This result strengthens FINITE-N-CONVENTION-ROBUSTNESS-20. It concerns only
finite-N modifications of the already-declared refresh rule and does not source
the rule, `g`, or a physical clock.

## Setup

Let `Q=12`, `u=1/Q`, retained space `V7`, hidden space `H4`, and let
`P_H` be the Euclidean hidden projector. Consider the target family

    q^B_{j|i}(p) proportional to
      exp[g (A7 p)_j - (g/N) B_{ji}],

where `B` is any fixed real 12x12 matrix. Column constants of `B` are invisible
to the softmax, so first replace `B` by its column-centered representative.
Let

    b = B 1

be its row-sum drift vector. For a retained real probe `phi`, put

    h = P_H(phi^2),          m = <phi^2>_u.

Let `C_B(phi)` denote

    lim N (L_full,B^3 - L_ME7,B^3) (sqrt(N) phi^T p)^4

at uniform equilibrium.

## Theorem

Relative to the empirical-refresh convention `B=0`, the leading quartic defect
is

    C_B(phi) - C_0(phi)
      = -(g/4) m b^T h
      = -3 g <phi^2>_u <P_H b, P_H(phi^2)>_u.       (1)

Consequences:

1. Only the hidden projection `P_H(B 1)` matters.
2. Every doubly-centered correction (`1^T B=0` and `B 1=0`) is exactly
   invisible at this order.
3. More generally, every row-sum drift in the retained space is invisible.
4. The correction map has rank at most four, exactly the dimension of the
   discarded k=1,2 sector.
5. Circulant corrections, the leave-one-out correction `B=A7`, `P_V`, and
   `P_H` are invariant as immediate corollaries, but translation invariance is
   much stronger than necessary.

Combining (1) with GENERAL-QUARTIC-THEOREM-13 gives

    C_B(phi)
      = -12 ||P_H(phi^2)||_u^2
        + 2g <P_H(phi^2), P_H(phi A7 phi)>_u
        - 3g <phi^2>_u <P_H(B1), P_H(phi^2)>_u.     (2)

Column-centering changes `B1` only by a constant-vector component, so the last
term can equivalently use the uncentered `B1` inside `P_H`.

## Proof

### 1. Split the correction into row-sum and doubly-centered pieces

After column centering, decompose

    B = B0 + b 1^T/Q,

where

    1^T B0 = 0,    B0 1 = 0.

The target correction enters the order-`1/N` generator as a common operator
`C_B` in both full and ME7 chains. Its first jump moment is proportional to
`B1=b`. Therefore the `B0` part has zero first jump moment and lowers polynomial
degree by at least two.

The full-minus-ME first-order departure operator is the hidden operator `D`.
For a visible quartic, the only B-dependent third-generator terms contain one
`G0`, one `D`, and one `C_B`, or the B-dependent part of `Delta G2`.

If `B=B0`, degree counting plus `V orthogonal H` kills all such terms:
- `C_B0` lowers degree by at least two;
- the relevant visible `D` terms already lose the degree needed to survive the
  final `G0` evaluation at zero;
- the B-dependent `Delta G2` drift is proportional to `B0 P_H d`; after two
  uniform `G0` contractions it contains the pairing of a hidden vector with
  the retained probe and vanishes.

Hence `B0` contributes exactly zero.

### 2. Rank-one row-sum representative

It remains to take

    B_{ji}=b_j/Q.

Write `a=phi.d` and `H=h.d`. The exact hidden departure operator satisfies

    [D a^4]_{degree 3} = 6 a^2 H.                  (3)

For the row-sum target correction, the degree-lowering-one part on a product of
linear forms is

    [C_b prod_r L_{x_r}]_{n-1}
      = -(g/Q^2) sum_r (b.x_r) prod_{s!=r} L_{x_s}. (4)

The uniform operator gives

    G0(a^2)=2m,
    [G0(6a^2 H)]_{degree 1}=12m H.                 (5)

Among the six permutations of `G0,D,C_b`, exactly three survive. They are the
three orderings in which `C_b` acts after the hidden departure has produced the
hidden linear factor. Each contributes

    -(12 g/Q^2) m (b.h).

The other three vanish because every surviving contraction contains
`phi.h=0`. Summing the three nonzero permutations gives

    -(36 g/Q^2) m (b.h).

For `Q=12` this is exactly

    -(g/4) m b.h,

which proves (1).

The B-dependent second-order target term gives no additional contribution: its
full-minus-ME part is linear in the hidden departure, and all three placements
with two `G0` operators reduce to a retained-hidden inner product.

## Exhaustive full-tensor numerical certificate

`continuation_08_quartic_Bclass.cpp` recomputes all 210 degree-four monomial
coefficients of the third-generator defect.

A basis of the complete 121-dimensional doubly-centered matrix space was used:

    (e_j-e_11)(e_i-e_11)^T,   i,j=0,...,10.

All 121 basis matrices reproduce the baseline full tensor. Worst residual over
all basis elements and all 210 monomials:

    constant coefficient: 0
    coefficient linear in g: 1.13243e-14
    g^2 coefficient: 0.

The 11-dimensional row-sum quotient has numerical rank exactly four. Fourier
row-sum directions k=3,4,5,6 are null to about 3.2e-14, while all four real
k=1,2 directions are nonzero.

For arbitrary random probes and the original noncirculant negative control,
formula (1) agrees with the complete tensor to <1.5e-15.

A sequence-resolved replay with a hidden k=1 row-sum drift finds that only the
three `G1-G1-G0` placements change. Each separately equals

    -(g/12) <phi^2>_u b^T P_H(phi^2)

for Q=12; the other three sequence classes are unchanged to roundoff.

## Methodological boundary

Promoted: exact finite-N-convention classification of the leading quartic
closure tensor inside the declared heat-bath + ME7 model.

Not promoted: a preferred finite-N convention, physical time, measured noise,
source for `g`, laboratory prediction, QW-2191, legacy-to-strict transfer,
`L_total`, SM/GR or ToE closure.
