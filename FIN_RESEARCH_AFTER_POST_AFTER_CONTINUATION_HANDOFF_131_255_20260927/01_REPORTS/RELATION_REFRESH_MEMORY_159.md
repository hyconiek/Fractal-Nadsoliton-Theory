# RELATION-REFRESH-MEMORY-159
## MaxEnt bond refresh has short O(1) memory but produces an annealed, not persistent, geometry

Date: 2026-09-26

Status:
exact one-edge dynamics for the minimal relation model of report 158.

Let each bond refresh at Poisson rate rho and, at a refresh, redraw from its
Gibbs Bernoulli target with probability p.

For one bond b(t):

    E[b]=p.

The centered autocovariance is exactly

    boxed:
    Cov[b(t),b(0)]
      =
      p(1-p) exp(-rho t).

So the relation-memory time is

    boxed:
    tau_rel=1/rho.

## 1. If the relation uses the same microscopic FIN clock

Set

    rho=1

in dimensionless FIN refresh units.

Then

    normalized bond memory
      =
      exp(-t).

Numerically:

    t=1:
      0.367879

    t=2:
      0.135335

    t=4:
      0.018316

    t=8:
      0.000335.

Thus the adjacency rewires on an O(1) microscopic time.

## 2. Consequence for emergent geometry

A geometry used to define locality for much slower effective unit dynamics
should normally persist over the unit-transition time.

The minimal bond heat bath does the opposite:
- memory is short;
- edges are annealed rapidly;
- no persistent neighborhood identity exists.

At long times the units see the average dense/block coupling, returning to a
mean-field description.

## 3. Slowing the bonds adds a new kinetic parameter

One can choose

    rho << 1

to make relations persistent.

But rho is then another unsourced dimensionless time-scale ratio.

This reproduces, at the relational level, the earlier kinetic nonuniqueness
problem.

## 4. Required improvement

A viable relational theory should generate slow bond persistence
ENDOGENOUSLY, for example from:
- a barrier in a relational potential;
- a conservation law;
- collective metastability of relation states.

The persistence time should follow from the same microscopic law, not from a
new fitted rho.

So the next test is not merely “can bonds exist?”

It is:

    can sparse and persistent bonds arise together without n-dependent
    or separately fitted parameters?
