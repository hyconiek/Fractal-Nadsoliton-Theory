# K5-ANGULAR-EVENT-CERTIFICATION-81
## Interval-assisted certification of the reflection-breaking crossing that seeds the generic g=5 branch

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
**PROVED_INTERVAL_ASSISTED_LOCAL_Z2_TRANSVERSE_CROSSING**
for the parent event;
the nonlinear quartic coefficient is currently high-precision numerical.

---

## 1. Exact reflection decomposition

For the reflection fixing the relevant k=5 parent branch, the seven-dimensional
retained space splits exactly as

    V = V_even direct-sum V_odd,

with

    dim V_even = 3,
    dim V_odd  = 4.

A convenient orthonormal even basis is

    q3+ = (3c-3s)/sqrt(2),
    q4+ = 4c,
    q5+ = (5c+5s)/sqrt(2),

while k6 lies in the odd block.

The parent branch is entirely inside V_even.

---

## 2. Crossing equations

Let s=(s3,s4,s5) denote the even coordinates.

The parent stationary equations are

    F_even(s,g)=0.

Let H_odd(s,g) be the four-dimensional reflection-odd Hessian block.

A reflection-breaking crossing satisfies

    F_even=0,
    det H_odd=0.

The numerical center is

    s3 =  0.157479386899271
    s4 = -0.204940006047350
    s5 = -0.441200652103440

    g  =  5.152672504944274.

---

## 3. Parametric Krawczyk proof

Use the full accepted outward strict spectral intervals and a radius

    2e-8

in every one of `(s3,s4,s5,g)`.

The determinant equation is rescaled by `10^6` for conditioning only.

The Krawczyk image is

    s3 in
      [0.157479386860259,
       0.157479386938073],

    s4 in
      [-0.204940006087089,
       -0.204940006007387],

    s5 in
      [-0.441200652147603,
       -0.441200652059058],

    g in
      [5.152672504916262,
       5.152672504972413].

The image lies strictly inside the source box.

Therefore a unique reflection-symmetric stationary crossing lies in the box,
uniformly over the accepted strict spectral rectangle.

---

## 4. It is transverse, not a fold of the parent branch

The interval determinant of the augmented four-by-four Jacobian is

    boxed:
    [-1.94593e-7,
     -1.93635e-7],

strictly separated from zero.

The even parent Hessian has three separated eigenvalue enclosures:

    [-0.00734264,-0.00734225],

    [ 0.01964185, 0.01964224],

    [ 0.05041040, 0.05041074].

So the parent branch itself is regular in g at the event.

The odd block has:

    one negative eigenvalue
      near -0.0211332,

    one simple zero eigenvalue,

    two positive eigenvalues
      near 0.0053926 and 0.0496177.

Thus the singularity is purely transverse to the reflection-fixed parent.

---

## 5. Crossing direction

Numerically, along the parent branch the critical odd eigenvalue satisfies

    d lambda_odd / dg
      ≈ -0.00472177260.

So as g increases through the crossing,

    lambda_odd:
      positive -> negative.

The parent full Morse index therefore changes

    boxed:
    2 -> 3.

---

## 6. Reduced nonlinear coefficient

Reflection symmetry forces the reduced coordinate y to enter evenly:

    Phi_red
      =
      Phi0
      +(lambda/2)y^2
      +(beta/4)y^4
      +O(y^6).

At the certified crossing center, Lyapunov-Schmidt elimination of the even
stable/unstable complement gives the numerical fourth derivative

    boxed:
    D4 Phi_eff
      ≈ -0.25769283695.

The sign is large and negative.

This numerical coefficient is consistent with, and explains, the observed
SUBCRITICAL reflection-breaking daughter.

A dedicated interval tensor contraction would be needed to promote the quartic
sign itself to interval-certified status; the crossing location/transversality
already is interval-assisted.

---

## 7. Daughter identification

Branch switching along the odd zero mode and pseudo-arclength continuation
gives the generic daughter on the low-g side.

That daughter:

    crosses g=5 as R7P-037 orbit i=14
      with index 3,

    reaches the generic fold
      g≈4.913081079101,

    returns through g=5 as orbit i=5
      with index 2.

The D12 matching errors at g=5 are below `1e-12`.

So the certified crossing is the source node for the previously isolated
generic g=5 component.

---

## 8. Updated local chain

The k=5 genealogy now contains

    reflection parent i=13-side, index 2
          |
          | CERTIFIED transverse Z2 crossing
          | g≈5.152672504944
          v
    generic daughter, index 3
          |
          | atlas i=14 at g=5
          v
    generic fold g≈4.913081079101
          |
    generic branch, index 2
          |
          | atlas i=5 at g=5
          v
    high-g continuation.

This is the first interval-assisted bridge connecting a reflection-symmetric
k=5 branch to a trivial-stabilizer branch in the FIN stationary graph.

---

## 9. Next research atom

`GENERIC-BRANCH-HIGH-G-DESTINATION-82`

Continue the i=5 high-g sheet with robust pseudo-arclength and detect:
- any later symmetry-restoring crossing;
- high-g asymptotic state concentration;
- possible connection to another atlas/symmetry component.

The goal is to determine whether the generic branch is a finite loop or an
asymptotic saddle family.
