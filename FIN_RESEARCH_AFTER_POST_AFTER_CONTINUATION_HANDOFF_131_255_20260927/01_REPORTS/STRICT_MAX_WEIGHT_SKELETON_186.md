# STRICT-MAX-WEIGHT-SKELETON-186
## In the declared strict vertex basis, the strongest coupling shell canonically recovers the internal C12 cycle

Date: 2026-09-26

Status:
exact finite-matrix observation in the fixed strict vertex representation.

The six strict shell weights obey

    w1 > w2 > ... > w6

for the current kernel values.

Numerically:


    d=1:
      w_d=0.469985672645020

    d=2:
      w_d=0.192043551690103

    d=3:
      w_d=0.091428614277925

    d=4:
      w_d=0.047029168745650

    d=5:
      w_d=0.024131223363630

    d=6:
      w_d=0.011070817321442


The maximal off-diagonal conductance occurs uniquely at cyclic distance

    d=1.

Keep only pairs carrying this maximal weight.

Every vertex then has exactly

    2

such neighbors, and the resulting support has

    12

edges.

It is exactly the 12-cycle C12.

## Consequence

Given:
- the strict matrix WITH its physical vertex basis;

the nearest cyclic skeleton is recoverable without separately being told which
pairs are nearest.

This is stronger than spectrum-only reconstruction.

It does not contradict the old isospectral obstruction:
an arbitrary isospectral rotation destroys the vertex basis and this
max-weight graph.

## Multicell boundary

The internal strict basis is not yet a set of physical simultaneous units.

So this result cannot simply be exported as spatial incidence.

But it suggests a concrete source template:

    sourced vertex basis
      +
    strongest-coupling relation
      ->
    two-port cycle skeleton.

The missing theorem is now the provenance of the multicell vertex basis /
transformation carrier, not the graph extraction once such a carrier exists.
