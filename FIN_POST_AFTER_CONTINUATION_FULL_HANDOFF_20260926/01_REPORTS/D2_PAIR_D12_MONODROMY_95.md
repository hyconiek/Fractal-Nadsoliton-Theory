# D2-PAIR-D12-MONODROMY-95
## The d=2 pair component connects two translated copies of the same large-g support orbit

Date: 2026-09-26

Status:
- global component tracing is numerical pseudo-arclength;
- the lower fold is well resolved numerically;
- the near-k5 internal cluster remains deliberately unresolved at proof level;
- endpoint D12 matching is numerical to ~1e-13.

## 1. Asymptotic endpoint

The large-g pair branch begins at

    S={0,2},
    p=(1/2,1/2),

with full Morse index 1.

Its D12 stabilizer is a single reflection, so the full orbit has size 12.

## 2. First fold

Descending in g reaches a reflection-fixed fold at

    boxed:
    g_lower ≈ 5.20080472903.

The full Hessian there has

    one negative,
    one zero,
    five positive

directions.

Thus the high-g index-1 pair sheet joins an index-2 sheet.

## 3. Near-k5 internal cluster

The index-2 sheet moves into the very flat k=5 regime near

    g5=12/lambda5
      ≈5.220554796949.

Pseudo-arclength sees an additional reflection-fixed turning structure within
about 10^-6 of g5, while a transverse angular eigenvalue is simultaneously
of order 10^-8--10^-9.

Because the k=5 angular stiffness starts only at order r^10, this is consistent
with report 79's degree-12 anisotropy.

No exact ordering of those almost coincident events is claimed here.

## 4. Second lower fold

Leaving the k5 cluster, the same component reaches a second copy of the
well-resolved lower fold at the same gain to numerical precision.

The second fold is related to the first by the half-turn translation T^6.

After it, the branch returns to large g with index 1.

## 5. Translated large-g endpoint

At exact g=100 the return sheet converges to

    S'={6,8},
    p=(1/2,1/2).

But

    S' = T^6 S.

Solving both finite-g stationary points at the same g=100 and applying the
exact D12 action T^6 gives feature-coordinate distance

    approximately 1.1e-13.

So the two large-g ends are not distinct support classes.

They are two representatives of the same D12 orbit.

## 6. Quotient topology

In the unreduced labelled state space the component has two asymptotic ends:

    {0,2}
      -> lower fold
      -> near-k5 internal segment
      -> translated lower fold
      -> {6,8}.

In the quotient by D12, the two ends are identified.

Therefore the component behaves as a stationary-orbit self-connection /
monodromy loop rather than a bridge between different zero-temperature
support classes.

This is qualitatively different from d=1 and d=6:

- d=1 has a transverse daughter that reaches a new size-3 support;
- d=6 has a transverse daughter that reaches another new size-3 support;
- d=2 returns to the same pair-support orbit through the weak k=5 angular
  structure.

## 7. Structural lesson

Support incidence alone does not determine branch topology.

The representation carried by the finite-g soft sector matters:

    d1:
      ordinary Z2 reflection breaking;

    d6:
      D2 -> Z2 transverse breaking plus exact k4×k6 factorization;

    d2:
      interaction with the nearly O(2) k=5 center manifold and its degree-12
      angular anisotropy.

This gives a concrete example of why the same asymptotic Morse index can have
different finite-g genealogies.
