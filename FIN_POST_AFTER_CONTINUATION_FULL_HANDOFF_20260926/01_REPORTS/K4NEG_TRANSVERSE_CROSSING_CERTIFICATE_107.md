# K4NEG-TRANSVERSE-CROSSING-CERTIFICATE-107
## Interval-assisted certification of the S4/index3 -> S5/index4 symmetry-breaking node

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
**PROVED_INTERVAL_ASSISTED_LOCAL_Z2_CROSSING_WITH_NEGATIVE_EFFECTIVE_QUARTIC**

The parent is the negative-k4 sheet of the exact k4×k6 factorized family.

Using the full accepted strict spectral intervals, solve the three equations:

    parent stationarity in (theta4c,theta6),
    det H_(k3c,k5c)=0.

A radius-4e-8 parametric Krawczyk box gives the strictly interior image

    theta4c:
      [-0.3261397735116484, -0.3261397735104238]

    theta6:
      [1.404596205249227, 1.4045962052519223]

    g:
      [5.76477155708275, 5.764771557084976].

Hence the transverse crossing is unique in that box.

The crossing gain is therefore enclosed near

    boxed:
    g≈5.764771557084.

The parent tangent block (k4c,k6) is positive definite.

At the crossing the full H7 inertia is exactly

    boxed:
    (3 negative, 1 zero, 3 positive).

The soft eigenvalue crosses transversely with

    boxed:
    d lambda_soft/dg
      in [0.037808550201977355, 0.03781120252673251]
      >0.

Reflection/half-turn symmetry makes the critical coordinate odd, so its reduced
potential is even:

    Phi_red
      = Phi0
        +(lambda/2)y^2
        +(D4_eff/24)y^4
        +O(y^6).

Using the exact Schur/Lyapunov-Schmidt formula

    D4_eff
      =
      D4 Phi[v^4]
      -3 b^T H_parent^(-1) b,

interval evaluation gives

    boxed:
    D4_eff
      in [-1.3423440773165296, -1.3422892627692584]
      <0.

So the crossing is rigorously subcritical locally.

Because the slope is positive, the broken-symmetry daughters exist on the
high-g side and gain one negative direction:

    parent index 3
      ->
    daughter index 4.

Numerical continuation from this certified local daughter reaches the exact
large-g support

    S5={2,4,6,8,10}

from report 102.

Thus the local S4/index3 -> S5/index4 edge now has an interval-assisted source
node and an interval-separated quartic sign.
