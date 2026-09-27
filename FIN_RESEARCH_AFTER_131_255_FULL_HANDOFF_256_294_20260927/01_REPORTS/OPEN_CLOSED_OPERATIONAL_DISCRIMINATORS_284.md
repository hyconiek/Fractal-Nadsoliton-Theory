# OPEN-CLOSED-OPERATIONAL-DISCRIMINATORS-284
## Fresh bath, cyclic register and local retained environment are operationally distinguishable at fixed rho

Date: 2026-09-27

Status:
exact finite/infinite model formulas.

This implements the AGENTS requirement to compare open and closed realizations
using return statistics, conservation and size scaling with a common calibrated clock.

Let

    tau=rho t.

Track the coefficient R(t) multiplying the initially tagged state in the one-site
transition kernel.

## 1. Fresh heat bath

For the open full reset process:

    boxed:
    R_fresh(t)
      =
      exp(-tau).

No finite closed-system count conservation is part of the subsystem model.

Long time:

    R_fresh -> 0.

## 2. Global cyclic register

For an L-record deterministic cyclic scheduler Poissonized at rate rho:

    boxed:
    R_reg,L(t)
      =
      P[
        Poisson(tau)=0 mod L
      ]

      =
      (1/L)
      sum_(m=0)^(L-1)
      exp{
        tau[
          exp(2 pi i m/L)-1
        ]
      }.

Global record multiset is exactly conserved.

For L>=3:

    R_reg,L(t)
      =
      exp(-tau)
      +
      O(
        tau^L/L!
      ).

So it agrees with the fresh bath through order L-1 in the short-time expansion.

Long time:

    R_reg,L -> 1/L,

with generally damped oscillatory Fourier contributions.

## 3. Local asynchronous SWAP cycle

For local edge swaps on C_n:

    boxed:
    R_swap,n(t)
      =
      (1/n)
      sum_m
      exp{
        -tau[
          1-cos(2 pi m/n)
        ]
      }.

Global color counts are exactly conserved.

For n>=3:

    boxed:
    R_swap,n(t)
      =
      1
      -
      tau
      +
      (3/4) tau^2
      +
      O(tau^3).

By contrast:

    R_fresh
      =
      1-tau+(1/2)tau^2+...

and for L>=3 the register has the same quadratic coefficient 1/2.

Therefore the local retained environment is distinguishable already at SECOND ORDER
in short-time return probability.

## 4. Example at equal rho, size 12

At tau=1:

    fresh:
      0.3678794412

    cyclic register:
      0.3678794419

    local SWAP cycle:
      0.4657596076.

So the cyclic register is still almost indistinguishable from a fresh bath at this
time, while the retained local environment has a large return excess.

At tau=4:

    fresh:
      0.0183156389

    register:
      0.0189571512

    local SWAP:
      0.2070023459.

## 5. Clean experimental/model discriminators

After calibrating rho independently:

### short-time return curvature

Distinguishes local retained environment from fresh/register.

### exact global count conservation

Distinguishes closed register/SWAP realizations from the open subsystem model.

### finite-size plateau

Closed finite models give:

    R -> 1/L or 1/n.

Fresh bath gives zero.

### approach to plateau

Register:
- complex cyclic modes can yield damped oscillatory structure.

Local SWAP:
- positive sum of real decaying heat-kernel modes.

### size scaling

Local SWAP slowest time:

    ~ n^2.

Global scheduler recurrence/mixing structure scales differently and is controlled by
the Poissonized permutation spectrum.

## Importance

These tests compare distinct CLOSED realizations that share the same one-event reset
intuition.

Therefore agreement with one-unit Q3 is insufficient to validate a particular
multiunit ontology.

The open/closed completion is empirically and mathematically distinguishable.
