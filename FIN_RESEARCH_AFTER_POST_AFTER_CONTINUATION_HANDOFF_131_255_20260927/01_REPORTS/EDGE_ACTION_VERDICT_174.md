# EDGE-ACTION-VERDICT-174
## The minimal node-edge-holonomy action fails the joint sparsity/persistence gate

Date: 2026-09-26

The constructive route of report 172 was deliberately minimal:

    node phase theta_x
      +
    binary edge b_xy
      +
    independent link phase a_xy
      +
    edge alignment
      +
    triangle holonomy.

It succeeds at one thing:

    nontrivial gauge-invariant curvature is mathematically available.

But it fails the two physical structure gates.

## Sparsity

Fixed O(1) couplings yield an edge probability that does not scale as 1/n.

Bounded-degree sparsity still needs:
- log(n) tuning; or
- a separate charge/valence constraint.

## Persistence

Fixed O(1) relational barriers imply O(1) relation lifetimes.

The already derived unit dynamics becomes much slower as internal N grows.

Thus the graph anneals before the effective unit state changes.

## Updated conclusion

The missing law cannot be merely

    "add a gauge field on links."

It must generate a relational sector whose:
- support is sparse in n;
- pair identity is persistent on the N-dependent coarse time scale;
- curvature/holonomy can be nontrivial;
- parameters are fixed by one common microscopic rule.

No tested minimal extension currently satisfies all four.

## Next atom

`PAIR-IDENTITY-CONSERVATION-175`

Investigate whether exact reversible bookkeeping can preserve pair identity
while allowing internal link state/phase dynamics.

The key test is whether pair identity can be a genuine conserved relational
charge without freezing an arbitrary graph as external initial data.
