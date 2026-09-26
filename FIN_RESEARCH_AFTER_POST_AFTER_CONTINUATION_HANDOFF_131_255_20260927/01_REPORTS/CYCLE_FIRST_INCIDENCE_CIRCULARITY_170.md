# CYCLE-FIRST-INCIDENCE-CIRCULARITY-170
## A cycle law cannot source adjacency unless an edge-support object exists first

Date: 2026-09-26

Status:
structural no-go.

A holonomy or curvature observable has the form

    H(C)
      =
      product_(e in C) U_e

or an additive analogue around a cycle C.

But a cycle is defined only after one knows which edges belong to the
incidence graph.

Therefore a rule such as

    "physical graphs are those satisfying a cycle-holonomy constraint"

already presupposes a candidate edge support.

## Complete-candidate workaround

One may instead place a link variable on EVERY pair, i.e. start from K_n.

Then all pair relations exist from the beginning.

A cycle action can select phases/fluxes on K_n, but the incidence remains
complete unless a separate activation variable b_xy is introduced.

Once b_xy is introduced, reports 158-167 apply:
- fixed local edge energies give dense activation;
- global edge count gives sparse but exchangeable activation;
- fixed valence gives sparse but nonunique graphs.

So cycle holonomy does not bypass the incidence problem.

It sits one layer above it.

## Logical order

The minimal typed order is:

    pair-support / bond variable
      ->
    link phase / connection
      ->
    cycle holonomy / curvature.

Reversing this order hides the missing support law.
