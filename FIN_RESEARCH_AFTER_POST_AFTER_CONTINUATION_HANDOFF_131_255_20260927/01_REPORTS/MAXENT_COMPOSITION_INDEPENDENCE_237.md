# MAXENT-COMPOSITION-INDEPENDENCE-237
## Extending the existing MaxEnt principle to fixed one-unit marginals selects independent composition

Date: 2026-09-26

Status:
exact entropy argument.

For a short interval h, let the one-unit transition kernel be

    P_h
      =
      I+h Q3+O(h^2).

Consider ALL joint two-unit kernels whose two marginals are both P_h.

For fixed marginals,

    H(X',Y')
      <=
    H(X')+H(Y'),

with equality iff X' and Y' are conditionally independent.

Therefore the unique maximum-joint-entropy kernel is

    boxed:
    P_h tensor P_h.

Expand:

    P_h tensor P_h
      =
      I
      +h(Q3 tensor I+I tensor Q3)
      +O(h^2).

Hence the continuous-time MaxEnt composition is

    boxed:
    L_ind
      =
      Q3 tensor I
      +
      I tensor Q3.

In the coupling simplex of report 236 this is exactly

    A=k,
    B=C=0.

## Consequence

The same maximum-entropy logic that naturally selected the single-unit refresh
does NOT generate an interaction when only the marginal laws are supplied.

It selects independence.

So:

    MaxEnt + known marginals
      ->
    no new collective slow scale.

Any nontrivial multi-unit coupling requires at least one additional relational
constraint beyond the one-unit transition law.

This prevents using "maximum entropy" by itself as a hidden source of kappa.
