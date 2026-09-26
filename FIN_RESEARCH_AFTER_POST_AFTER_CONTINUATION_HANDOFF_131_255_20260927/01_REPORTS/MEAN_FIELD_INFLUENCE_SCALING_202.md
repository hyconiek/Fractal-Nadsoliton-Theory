# MEAN-FIELD-INFLUENCE-SCALING-202
## Pair influence is O(1/N), while the number of influenced partners is O(N)

Date: 2026-09-26

Status:
exact softmax bound.

For a source-label change r->s, define

    delta
      =
      (g/N) A7(e_s-e_r).

Let

    osc(delta)
      =
      max_i delta_i-min_i delta_i.

For two softmax distributions whose log-density ratio has oscillation Delta,

    TV
      <=
      tanh(Delta/4).

Therefore every one-copy intervention obeys

    boxed:
    I_pair
      <=
      tanh[
        g Delta_A/(4N)
      ],

where

    Delta_A
      =
      max_(r,s)
      osc[A7(e_s-e_r)].

For the current A7,

    Delta_A
      ≈ 3.964066624707.

Hence

    I_pair
      <=
      tanh[
        5.099007360852/N
      ]
      =
      O(1/N).

## Mean-field balance

Each copy has

    N-1

potential partners.

Therefore

    (N-1) I_pair

can remain O(1) as N grows.

This is the characteristic weak-all-to-all scaling of a mean-field model.

It is not finite-degree locality.

## Interpretation

The finite-copy architecture has a controlled thermodynamic limit precisely
because individual influences weaken while their number grows.

That success is mathematically opposite to the desired spatial-locality
mechanism.

So the internal copy index should not be promoted to a physical spatial cell
index.
