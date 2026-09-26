# TWO-HARMONIC-FOLD-CERTIFICATE-73
## Interval-assisted certification of the fold connecting R7P-036 to the MP7-039 parent

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
**PROVED_INTERVAL_ASSISTED_LOCAL_TWO_HARMONIC_FOLD**
within the accepted strict spectral intervals.

This upgrades the numerical edge found in report 72.

---

## 1. Exact reduced system

Use the exact R7P-034 family

    h_j = J cos(pi j/2)+K(-1)^j.

Define

    x = sinh(J)/(cosh(J)+exp(-2K)),
    y = [cosh(J)-exp(-2K)]/
        [cosh(J)+exp(-2K)].

The stationary equations are

    F1 = J-g lambda3 x/6 = 0,
    F2 = K-g lambda6 y/12 = 0.

Let

    A = d(F1,F2)/d(J,K).

A fold satisfies

    det A = 0.

So the augmented fold system is

    G(J,K,g)
      =
      (F1,F2,det A)
      =0.

---

## 2. Parametric Krawczyk certificate

The calculation keeps the full accepted outward intervals for

    lambda3,
    lambda6.

Around the numerical point

    J0 = 0.718505246005629...
    K0 = 0.490908854614046...
    g0 = 4.621196599489125...

use the box with radius `2e-7` in each of `(J,K,g)`.

A midpoint inverse of the augmented 3x3 Jacobian is used only as the Krawczyk
preconditioner; all function and derivative evaluation on the box is interval.

The Krawczyk image is

    J in
      [0.718505245955803,
       0.718505246055457],

    K in
      [0.490908854587306,
       0.490908854640786],

    g in
      [4.621196599487042,
       4.621196599491209].

This image lies strictly inside the radius-`2e-7` source box.

Therefore, uniformly over the accepted `(lambda3,lambda6)` spectral rectangle,
there is a unique fold root in the declared box.

---

## 3. Simple-fold transversality

Let `v` and `w` be unnormalized right/left null directions of the 2x2
stationary Jacobian A.

A convenient algebraic choice for

    A=[[a,b],[c,d]]

at `det A=0` is

    v=(-b,a),
    w=(-c,a).

The parameter and quadratic fold coefficients are

    a_fold = w . F_g,

    b_fold =
      w . D^2_(J,K)F[v,v].

Interval evaluation over the full fold box gives

    boxed:
    a_fold
      in
      [-0.049019977707,
       -0.049019533007]
      <0,

    boxed:
    b_fold
      in
      [0.015748835983,
       0.015749599644]
      >0.

Hence the fold is simple/nondegenerate.

The signs are stated in this unnormalized null-vector convention; rescaling
the null vectors changes magnitudes but not nonvanishing.

---

## 4. Full-X7 transverse inertia

At the fold the two symmetry-related `(k4,k5)` Hessian blocks satisfy

    trace(H45)
      in
      [0.057939668248,
       0.057941155972]
      >0,

    det(H45)
      in
      [-0.007120350091,
       -0.007120181294]
      <0.

Therefore each H45 block contains exactly one negative and one positive
eigenvalue.

The independent `k3-sine` direction satisfies

    H_3sin
      in
      [0.141903259785,
       0.141903398758]
      >0.

The `(k3-cos,k6)` block supplies the single fold zero mode and one positive
hard mode.

Thus the full H7 inertia exactly at the fold is

    boxed:
    (2 negative, 1 zero, 4 positive).

---

## 5. Local branch exchange

Because the fold is simple and the two transverse H45 negative directions stay
strictly separated from zero, the two local stationary sheets have full Morse
indices

    boxed:
    2 and 3.

This rigorously upgrades the numerical identification:

    large two-harmonic sheet
      <-> simple fold
      <-> small two-harmonic sheet.

At g=5 these are precisely the known orbit-size-4 parent states:
- the R7P-036 index-2 counterexample sheet;
- the small index-3 sheet that proceeds toward MP7-039.

---

## 6. Certified structural chain

The following part of the branch graph is now supported by local interval
theorems at both special nodes:

    R7P-036 large parent
        index 2
          |
          | CERTIFIED simple fold
          | g≈4.621196599489
          v
    small two-harmonic parent
        index 3
          |
          | CERTIFIED D3 crossing MP7-039
          | g≈5.171841831943
          v
    D3 branch reorganization.

The long continuation between the two certified neighborhoods is still
numerical; the endpoint events themselves are interval-assisted.

---

## 7. Next certification target

The highest-value new numerical node is now

    g_upper≈5.172231474684,

only `3.8964e-4` above the D3 crossing.

`UPPER-C4-FOLD-CERTIFICATE-74` should certify that fold and thereby validate the
local join

    D3 daughter index 2
      -> upper fold
      -> main localization saddle index 1.

That would create an interval-certified bridge from the D3 event into the
stationary component containing the first localization fold.
