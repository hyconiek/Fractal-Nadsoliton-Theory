# LARGE-N-METASTABILITY-260
## Exact D4 quotient extends deep-core capacity to N=12 and the local exponent continues toward B4

Date: 2026-09-27

Status:
exact finite-state potential-theory result for N=12;
asymptotic interpretation remains numerical.

This report continues the parallel one-unit metastability lane independently of the multiunit research.

## 1. Exact N=12 quotient calculation

Full count-state size:

    1,352,078.

Exact D4 quotient orbit count:

    169,532.

The quotient self-adjoint generator residual is below

    3.0e-14.

The Dirichlet CG solve converged with

    info=0.

For the same four-deep-seed Z3 source core and the other two Z3 deep cores as target:

    capacity
      =
      0.000128571234764.

Using the same 3*capacity rate convention as reports 233:

    boxed:
    r_12
      =
      0.000385713704291.

## 2. New local exponent

With

    r_11
      =
      0.000726643165759,

the new finite-size exponent is

    boxed:
    beta_(11->12)
      =
      -log(r_12/r_11)
      =
      0.633340130350.

The mapped inter-Z3 communication barrier is

    B4
      =
      0.662219137127.

The difference is now

    B4-beta
      =
      0.028879006777.

Previous local exponents were:

    6->7:   0.540297
    7->8:   0.563393
    8->9:   0.585834
    9->10:  0.605172
    10->11: 0.620903
    11->12: 0.633340.

The monotone approach toward B4 continues.

## 3. Descriptive trend only

For these six points, the deficit is very well described numerically by

    B4-beta_N
      ~
      0.98 exp(-0.28885 N_upper),

with

    R^2≈0.9919.

This is ONLY a descriptive finite-N fit.

It is not a controlled asymptotic remainder theorem.

## Updated verdict

The evidence that B4 is the large-N switching exponent is materially stronger than after N=11.

But the proof boundary is unchanged:

1. exclude every lower global communication route;
2. derive capacity asymptotics with controlled prefactor/remainder;
3. avoid interpreting monotone numerical convergence as proof.

The one-unit lane is progressing in the expected direction.
