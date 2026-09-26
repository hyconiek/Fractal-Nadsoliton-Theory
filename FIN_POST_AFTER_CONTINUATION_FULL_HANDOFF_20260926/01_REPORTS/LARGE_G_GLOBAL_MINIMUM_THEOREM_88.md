# LARGE-G-GLOBAL-MINIMUM-THEOREM-88
## The twelve one-label branches are the asymptotic global minima

Date: 2026-09-26

Status:
- analytic large-g variational theorem for the supplied rank-seven mediator;
- physical interpretation remains limited to the FIN probability landscape.

Let x_i denote the i-th row of X7.

Then

    A7_ij = x_i . x_j

and

    p^T A7 p
      = || sum_i p_i x_i ||^2.

All twelve row vectors have the same squared norm

    ||x_i||^2
      = A7_ii
      = 1.271778833994...

and all twelve rows are distinct.

By convexity / the triangle inequality,

    ||sum_i p_i x_i||
      <= sum_i p_i ||x_i||
      = sqrt(A7_ii).

Hence

    boxed:
    p^T A7 p <= A7_ii.

Equality requires all active row vectors to point in exactly the same direction.
Because all have the same norm and are distinct, equality is possible only at
a simplex vertex

    p=delta_i.

Therefore the interaction part of

    V_g(p)
      =
      D(p||uniform)
      -(g/2)p^T A7 p

has exactly twelve asymptotic maximizers of the quadratic term, one for each
label.

By report 85 each vertex generates a unique finite-g interior stationary branch
for sufficiently large g.

Compactness plus the strict interaction gap away from the vertices implies that
for sufficiently large g the global minima of V_g lie in these twelve vertex
neighborhoods.

D12 acts transitively on them, so they have equal energy.

Thus:

    boxed:
    for sufficiently large g,
    the global-minimum set is one 12-member D12 orbit
    converging to the twelve one-label vertices.

The limiting minimum energy is

    V_min(g)
      =
      -(g/2) A7_00
      +log 12
      + exponentially small corrections,

because the isolated one-label branch has no algebraic 1/g correction.

A quantitative identity also holds:

    A7_00 - p^T A7 p
      =
      (1/2) sum_{i,j}
      p_i p_j ||x_i-x_j||^2.

The minimum pairwise row distance is approximately

    1.483486341382.

So concentration away from a vertex carries an explicit O(g) interaction
penalty.

This provides the large-g endpoint of the stable localized branch genealogy:
the finite-g localized minimum ultimately becomes a single-label state in the
zero-temperature limit.
