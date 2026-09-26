# SWAP-CYCLE-DIFFUSION-241
## Closing the Z3 refresh environment into a local swap network produces an exact n^-2 collective mode

Date: 2026-09-26

Status:
exact first-moment theorem;
exact finite-state spectral replay in balanced sectors for n=3,6,9.

## 1. Open one-unit law

The established effective generator is

    Q3
      =
      rho(U-I).

At Poisson rate rho the unit receives a fresh uniform trit.

Report 240 shows that the unique reversible symmetric deterministic dilation of
one reset is:

    system trit <-> environment trit.

## 2. Close the environment

Now take n simultaneous Z3 units on a cycle.

On each edge perform a swap at rate

    r_edge=rho/2.

The factor 1/2 is the equal split of one unit's total activity budget rho among
its two equivalent incident edges.

Every elementary state transformation is a bijection.

The full multiset of Z3 labels is conserved.

The product-uniform distribution is stationary, though the process decomposes
into fixed-count sectors.

## 3. Exact density equation

For any single-site observable f and site i,

    d/dt E[f(X_i)]
      =
      r_edge[
        E f(X_(i+1))
        +E f(X_(i-1))
        -2 E f(X_i)
      ].

Thus the one-point density obeys the cycle heat equation EXACTLY.

Fourier mode m has decay rate

    boxed:
    lambda_m
      =
      2 r_edge[
        1-cos(2 pi m/n)
      ]

      =
      rho[
        1-cos(2 pi m/n)
      ].

The slowest nonconstant density mode is

    lambda_1
      ~
      boxed:
      2 pi^2 rho / n^2.

This is the first current multi-unit FIN candidate with a collective time scale
that grows with system size.

## 4. Exact finite-state replay

Set rho=1 and restrict to balanced color-count sectors.

Numerical full-generator gaps are:

    n=3:
      1.500000000000

    n=6:
      0.500000000000

    n=9:
      0.233955556881.

They match

    1-cos(2 pi/n)

to numerical precision.

## 5. Scale separation

Relative to the isolated Z3 relaxation time 1/rho,

    tau_cycle/tau_single
      =
      1/[1-cos(2 pi/n)].

Examples:

    n=12:
      7.464

    n=24:
      29.348

    n=64:
      207.673.

So the hierarchy is asymptotically

    tau_cycle/tau_single
      ~
      n^2/(2 pi^2).

No new continuous kappa was used.

## Boundary

This result is conditional on:
- cycle incidence;
- equal splitting of the local activity budget;
- using the swap dilation as the closed multi-unit realization.

It gives dimensionless diffusion on a graph, not physical SI space.
