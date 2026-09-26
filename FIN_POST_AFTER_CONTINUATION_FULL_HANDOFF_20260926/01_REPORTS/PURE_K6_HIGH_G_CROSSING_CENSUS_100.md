# PURE-K6-HIGH-G-CROSSING-CENSUS-100
## Exact transverse crossing census and Morse-index staircase on the pure-k6 trunk

Date: 2026-09-26

Status:
- analytic mode-by-mode Hessian formulas on the exact pure-k6 branch;
- one scalar numerical root for the k3c crossing;
- exact k5 and k4 crossing gains;
- exact no-crossing arguments for the k3s and radial-k6 directions.

## 1. Exact pure-k6 branch

Write

    h_j = J (-1)^j,
    q = tanh J.

Stationarity is

    J = (g lambda6/12) q.

The nonzero branch is born at

    boxed:
    g6 = 12/lambda6
       = 5.123427551398616.

The uniform k6 quartic is positive, so this is the supercritical parity branch.

As g->infinity it converges to the six-label support

    {0,2,4,6,8,10}

or its odd translate, with weight 1/6 on each label.

That support has affine rank 5, hence asymptotic Morse index 5.

## 2. Exact transverse Hessian formulas

Under the parity-biased law:

### k3

    h_3c
      = 1/g
        -(lambda3/12)(1+q),

    h_3s
      = 1/g
        -(lambda3/12)(1-q).

### k4

Both real directions are degenerate:

    h_4
      = 1/g-lambda4/12.

### k5

Both real directions are degenerate:

    h_5
      = 1/g-lambda5/12.

### radial k6

    h_6
      = 1/g
        -(lambda6/12)sech^2 J.

So the complete transverse spectrum is available in closed form.

## 3. First secondary crossing: k3c

Solve

    J=(g lambda6/12)tanh J,

    1/g=(lambda3/12)(1+tanh J).

The unique nonzero crossing is

    boxed:
    g_3|6
      ≈ 5.180490619637444.

Below it the pure-k6 branch has index 0.
Above it the branch has index 1.

Report 97 showed that the associated subcritical daughter is exactly the
k3+k6 two-harmonic D3 family.

## 4. Second crossing: exact double k5 event

Because h_5 is independent of the k6 bias,

    boxed:
    g_5|6
      = 12/lambda5
      = 5.220554796949195.

Two eigenvalues cross simultaneously.

Therefore the pure-k6 parent changes

    index 1 -> index 3.

Report 96 classified this as the exact D6 event producing two inequivalent
twelve-state daughter orbits.

## 5. Third crossing: exact double k4 event

Likewise,

    boxed:
    g_4|6
      = 12/lambda4
      = 5.455614632675739.

Again two real directions vanish simultaneously.

The pure-k6 trunk therefore changes

    index 3 -> index 5.

The k4 critical representation of the D6 parent factors through D3.

The cubic coefficient along k4c is numerically

    T3(k4c,k4c,k4c)
      ≈ -0.055490606789568,

strictly nonzero.

Thus this is another cubic D3-type transverse event, not a pitchfork.

Its local daughter amplitude is linear:

    r
      ~ 1.210941466 |g-g_4|.

## 6. No k3s crossing

Using stationarity,

    1/g=(lambda6/12) q/J.

Thus

    12 h_3s
      = lambda6 q/J
        -lambda3(1-q).

Define

    L(J)=q/[J(1-q)].

For J>0,

    d/dJ log L
      = (1+q)/q - 1/J
      >0

because tanh J=q<J.

Also

    lim_(J->0+) L(J)=1.

Since

    lambda3/lambda6 <1,

we have

    L(J)>lambda3/lambda6

for every J>0.

Hence

    boxed:
    h_3s>0

on the entire nonzero pure-k6 branch.

## 7. No radial-k6 recrossing

Similarly,

    h_6
      = (lambda6/12)
        [q/J-(1-q^2)].

For J>0,

    q-J(1-q^2)>0,

because this function vanishes at J=0 and has derivative

    2J q(1-q^2)>0.

Therefore

    boxed:
    h_6>0

everywhere on the nonzero branch.

## 8. Exact Morse staircase

The complete pure-k6 Morse-index sequence is therefore

    g6 < g < g_3|6:
        index 0

    g_3|6 < g < g_5|6:
        index 1

    g_5|6 < g < g_4|6:
        index 3

    g > g_4|6:
        index 5.

So:

    boxed:
    0 -> 1 -> 3 -> 5.

The terminal value 5 exactly matches the support theorem

    index = |S|-1

for the six-label large-g support.

## 9. Structural interpretation

The pure-k6 branch is a central high-g trunk.

It starts as a stable parity-broken state, then sequentially exposes:
- one k3 instability;
- a two-dimensional k5 instability;
- a two-dimensional k4 instability.

The increment in Morse index equals the dimension of each newly unstable
critical representation:

    +1, +2, +2.

This gives a direct bridge between

    representation dimension
and
    growth of saddle complexity.

## 10. Next atom

`PURE-K6-k4-DAUGHTER-101`

Classify the nonlinear daughters of the exact g=12/lambda4 crossing and match
them to the size-4/size-5 high-g support atlas.

Because the cubic D3 invariant is already nonzero, this should be much easier
than the degree-12 k5 problem.
