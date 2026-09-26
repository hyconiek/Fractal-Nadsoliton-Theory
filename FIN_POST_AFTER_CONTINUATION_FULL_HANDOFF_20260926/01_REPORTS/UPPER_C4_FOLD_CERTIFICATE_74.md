# UPPER-C4-FOLD-CERTIFICATE-74
## Interval-assisted fold joining the D3 daughter to the main localization saddle

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
**PROVED_INTERVAL_ASSISTED_LOCAL_C4_FOLD**

This upgrades the numerical upper-fold node of reports 69/71.

---

## 1. C4 stationary fold system

Use the reflection-even C4 coordinates

    s=(s3,s4,s5,s6)

and the full supplied four-mode log-partition potential.

The stationary map is

    F(s,g)
      = s/g - E_s[C4].

Its Jacobian is the C4 Hessian

    H = dF/ds
      = I/g - Cov_s(C4).

A fold satisfies

    F=0,
    det H=0.

The augmented system therefore has five variables `(s3,s4,s5,s6,g)`.

---

## 2. Numerical center

The candidate found by pseudo-arclength continuation is

    s3 = 0.064551183234886
    s4 = 0.005557251343608
    s5 = 0.013192678416272
    s6 = 0.418079565489821

    g  = 5.172231474684087.

It lies only

    g-g_D3 ≈ 0.000389642741

above the certified MP7-039 D3 crossing.

---

## 3. Parametric interval Krawczyk

The calculation retains the full accepted outward intervals for

    lambda3,lambda4,lambda5,lambda6.

Use a radius

    2e-8

in every one of the five fold variables.

The determinant equation is rescaled by a fixed positive factor `10^6`
for numerical conditioning; this does not change its zero set.

The Krawczyk map is strictly interior in all five coordinates.

A representative interval image is approximately:

    s3 in
      [0.0645511832272,
       0.0645511832425],

    s4 in
      [0.0055572513144,
       0.0055572513727],

    s5 in
      [0.0131926783476,
       0.0131926784845],

    s6 in
      [0.418079565479,
       0.418079565500],

    g in
      [5.1722314746813,
       5.1722314746869].

Thus, uniformly over the accepted strict spectral rectangle, a unique singular
C4 stationary point lies in the source box.

---

## 4. Augmented Jacobian is nonsingular

Direct interval evaluation of the determinant of the five-by-five augmented
Jacobian gives

    det D(F,10^6 det H)
      in approximately

      [-9.6703e-13,
       -9.5722e-13].

It is strictly negative.

At the same time the C4 Hessian has exactly one fold zero mode, with all other
C4 directions separated from zero.

For a gradient stationary problem with a one-dimensional Hessian kernel,
nonsingularity of the augmented `(F,det H)` Jacobian is the generic simple-fold
condition: both parameter transversality and the quadratic fold coefficient are
nonzero.

So this is a genuine simple C4 fold, not a higher-order cusp.

---

## 5. Certified nonzero C4 complement inertia

Rotate the interval C4 Hessian by the numerical orthonormal eigenbasis at the
fold and apply interval Gershgorin bounds on the three-dimensional complement
of the zero mode.

The three nonzero C4 eigenvalue enclosures are separated as

    negative:
      [-0.00277382,
       -0.00277337],

    positive:
      [0.00745304,
       0.00745347],

    positive:
      [0.01158932,
       0.01158972].

Thus the C4 fold has

    one negative,
    one zero,
    two positive

directions.

---

## 6. Odd/reflection-transverse block

The three sine-sector directions have interval-separated positive eigenvalues:

    [0.00032062,0.00032085],

    [0.01190634,0.01190656],

    [0.05980295,0.05980307].

The smallest is narrow but strictly positive throughout the fold box.

Therefore the full X7 Hessian at the upper fold has exact inertia

    boxed:
    (1 negative, 1 zero, 5 positive).

---

## 7. Local branch-index exchange

Because the fold is simple and every transverse eigenvalue is separated from
zero, the two local sheets have full Morse indices

    boxed:
    2 and 1.

This is exactly the numerical transition observed in pseudo-arclength:

    D3 daughter sheet
        index 2
          |
          | CERTIFIED upper C4 fold
          | g≈5.172231474684
          v
    main saddle sheet
        index 1.

---

## 8. Connection to localization

Numerical continuation of the index-1 sheet from this certified fold reaches
the already certified simple localization fold

    g_f≈3.51564471684,

where it joins the stable localized minimum branch.

Thus the two ends of the relevant saddle segment now have interval-certified
local events:

    upper C4 fold (this report)
          |
       numerical continuation
          |
    lower localization fold (R7P-031/MP7-026).

The long connecting segment is not yet a validated continuation tube, but its
endpoints and local branch indices are certified.

---

## 9. Updated proof-level chain

The following event sequence is now locally interval-assisted:

    two-harmonic fold
      g≈4.621196599489
          |
    small parent index 3
          |
    D3 crossing
      g≈5.171841831943
          |
    daughter index 2
          |
    upper C4 fold
      g≈5.172231474684
          |
    saddle index 1
          |
    localization fold
      g≈3.51564471684
          |
    localized minimum index 0.

The branch segments between isolated event neighborhoods remain numerical
unless already covered by an upstream continuation theorem.

---

## 10. Next high-value certification

The remaining new fold from report 70 is

    g≈4.395526393543

on the second D3 daughter component.

`SECOND-C4-FOLD-CERTIFICATE-75` should apply the same 5D interval machinery
after the appropriate D12 transformation to the signed reflection-fixed C4
representative.

If successful, both D3 daughter triplets will have certified first fold
destinations.
