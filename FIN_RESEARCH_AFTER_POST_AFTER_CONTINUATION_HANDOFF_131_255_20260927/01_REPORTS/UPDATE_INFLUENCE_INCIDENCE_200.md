# UPDATE-INFLUENCE-INCIDENCE-200
## Incidence defined by actual leave-one-out influence is complete, not local

Date: 2026-09-26

Status:
exact structural result plus numerical local-update replay.

Consider the N labeled microscopic copies inside ONE leave-one-out FIN unit.

For a target copy b, its redraw distribution is

    q_b
      =
      softmax[
        (g/N) A7 m^(-b)
      ].

Change any other copy a from label r to label s.

Then the target logits change by

    boxed:
    delta h
      =
      (g/N) A7(e_s-e_r).

Because A7 has distinct columns, there exist r,s for which delta h is not a
constant vector.

Softmax is invariant only under addition of a constant vector.

Therefore the conditional law of b changes.

Hence EVERY distinct pair

    a != b

has nonzero direct dynamical influence.

## Incidence graph

If adjacency is defined by

    a~b
      iff
    changing a can alter the next-update law of b,

then the graph on microscopic copies is exactly

    boxed:
    K_N.

This graph is derived from the actual microscopic dynamics, not inserted.

But it is complete.

So the most direct operational influence definition recovers mean-field
all-to-all coupling rather than sparse locality.

## Numerical replay

Maximum pairwise total-variation influence over backgrounds:


    N=2:
      I_max=0.797380881024

    N=3:
      I_max=0.707029553403

    N=4:
      I_max=0.584847456572

    N=5:
      I_max=0.500485745235

    N=6:
      I_max=0.439920025168

    N=7:
      I_max=0.371796678586

    N=8:
      I_max=0.347719575628


The pair influence weakens with N but never selects a subset of neighbors.

All source-copy identities remain exchangeable.
