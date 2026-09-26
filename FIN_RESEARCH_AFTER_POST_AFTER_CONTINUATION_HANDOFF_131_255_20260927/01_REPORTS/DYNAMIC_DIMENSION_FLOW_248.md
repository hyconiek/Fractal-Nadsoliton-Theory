# DYNAMIC-DIMENSION-FLOW-248
## FIN now contains a finite state-space example where the effective number of elementary directions changes with the accessible transition scale

Date: 2026-09-26

Reports 245-247 separate three concepts:

1. algebraic group factorization;
2. chosen elementary transformation set;
3. dynamically accessible transformation set.

They are not the same.

## Current Z12 example

At the lowest nontrivial barrier level:

    accessible generator:
      d3 only;

    geometry:
      three copies of C4;

    local generator count:
      one.

At the connectivity threshold:

    accessible generators:
      d3,d4;

    geometry:
      C3 square C4;

    minimum-bottleneck generator count:
      two.

At the next barrier:

    d5

adds a diagonal/mixed edge but no new independent direction.

Thus the finite state-space picture has the schematic flow

    isolated points
      ->
    1-direction components
      ->
    2-direction connected product
      ->
    same 2-direction geometry with extra diagonals.

## Why this matters

It gives a mathematically explicit version of:

    "dimension can depend on which transformations are dynamically accessible."

The number of elementary directions need not be hard-coded in the vertex set.

It can, in principle, be selected by a cost/barrier hierarchy.

## General multicell implication

Suppose a future multicell FIN law derives:
- a family of commuting reversible transformations;
- their transition costs/barriers.

Then one can define the connectivity threshold

    h_c
      =
      inf{
        h:
        transformations with cost <=h
        generate one connected orbit
      }.

A candidate effective dimension can then be associated with the number of
independent elementary generators needed at h_c.

This is a source strategy, not yet a physical theorem.

## Important caveat

A cyclic group can be connected by one expensive generator or by several
cheaper generators.

Therefore "dimension" is not a pure group invariant.

It is a property of:

    group
      +
    elementary-generator notion
      +
    dynamical cost scale.

This is likely the correct level at which the next FIN geometry campaign must
operate.
