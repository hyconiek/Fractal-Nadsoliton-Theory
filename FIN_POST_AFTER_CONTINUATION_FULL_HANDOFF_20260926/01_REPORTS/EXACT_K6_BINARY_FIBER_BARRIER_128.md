# EXACT-K6-BINARY-FIBER-BARRIER-128
## The uniform -> ±k6 pitchfork gives an exact FIN-sourced binary metastable barrier

Date: 2026-09-26

Status:
- exact one-dimensional potential and stationary parameterization;
- exact barrier formula;
- finite-N exponential rate follows from the exact Gibbs realization;
- no Eyring-Kramers prefactor or physical seconds are claimed.

## 1. Exact pure-k6 potential

Let

    J=A6 theta6,
    A6^2=lambda6/12.

Along the pure-k6 line,

    Phi(J,g)
      =
      J^2/[2g A6^2]
      -log cosh(J).

The stationary equation is

    J
      =
      g A6^2 tanh J.

The critical gain is

    boxed:
    g6=1/A6^2
      =5.123427551398616.

For g>g6 there are two stable stationary points ±J and the uniform point J=0
is the intervening index-one saddle.

## 2. Exact parametric branch

At a nonzero stationary point,

    boxed:
    g/g6
      =
      J/tanh J.

So J itself is an exact parameter for the binary daughter branch.

## 3. Exact saddle barrier

Since Phi(0,g)=0, the barrier from either daughter to the uniform saddle is

    Delta Phi_6
      =
      -Phi(J,g).

Using stationarity,

    boxed:
    Delta Phi_6(J)
      =
      log cosh J
      -(J/2)tanh J.

This is positive for J!=0.

No numerical continuation is required.

## 4. Near-critical scaling

Define

    epsilon=(g-g6)/g6.

Then

    J/tanh J
      =
      1+J^2/3+O(J^4),

so

    J^2
      =
      3 epsilon
      +O(epsilon^2).

Also

    Delta Phi_6
      =
      J^4/12+O(J^6).

Therefore

    boxed:
    Delta Phi_6
      =
      (3/4) epsilon^2
      +O(epsilon^3).

Equivalently,

    Delta Phi_6
      =
      (3/4) A6^4 (g-g6)^2
      +....

This is a pitchfork barrier exponent 2, different from the simple-fold
3/2 exponent found earlier.

## 5. Exact finite-N meaning

R7P-111 supplies the N-copy Gibbs measure with rate function V_g, and the
heat-bath chain is reversible with respect to it.

Therefore the ±k6 switching exponent is

    boxed:
    exp[-N Delta Phi_6]

up to subexponential factors.

The finite-size barrier crossover obeys

    N Delta Phi_6=O(1),

hence

    boxed:
    epsilon=O(N^(-1/2)).

So the binary pitchfork fiber has a natural finite-N critical window

    (g-g6)/g6 ~ N^(-1/2).

## 6. Fiber conductance

Let k_f be the effective transition rate from +k6 to -k6 after eliminating the
uniform saddle region.

On the exponential N scale,

    -log k_f / N
      =
      Delta Phi_6 + o(1).

A symmetric two-state chain has Laplacian

    k_f L2.

Therefore, at exponential accuracy, the ST231 fiber coefficient can be
identified with

    boxed:
    mu_f ~ k_f
         ~ exp[-N Delta Phi_6].

This is the first FIN-internal mechanism that supplies BOTH:
- a natural binary child incidence;
- a metastable fiber-rate exponent.

The unresolved piece is the subexponential prefactor / global clock.
