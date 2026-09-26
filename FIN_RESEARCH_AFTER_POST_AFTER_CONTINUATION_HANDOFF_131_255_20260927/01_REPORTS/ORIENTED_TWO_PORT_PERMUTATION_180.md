# ORIENTED-TWO-PORT-PERMUTATION-180
## One-in/one-out conservation gives cycle covers, but MaxEnt does not select one connected cycle

Date: 2026-09-26

Status:
exact permutation-combinatorics theorem.

A stronger local relational constraint is:

    each unit has exactly one outgoing relation
    and exactly one incoming relation.

Then the oriented relation state is a permutation of n units.

Every permutation decomposes into disjoint directed cycles.

So the local charge automatically gives:
- degree two in the underlying undirected sense;
- orientation on each component.

But it does NOT guarantee connectedness.

## Uniform MaxEnt over permutations

There are

    n!

permutations.

The number consisting of one n-cycle is

    (n-1)!.

Therefore

    boxed:
    P(single connected cycle)
      =
      1/n.

As n grows, connected single-cycle incidence becomes LESS likely.

The expected number of cycle components is

    H_n
      =
      1+1/2+...+1/n
      ~ log n+gamma.

So the natural MaxEnt one-in/one-out ensemble fragments into multiple cycles.

## Consequence

A local two-port conservation law is not enough.

The connectedness condition of report 178 needs its own source.
