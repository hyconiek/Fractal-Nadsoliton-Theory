# SECOND-C4-FOLD-CERTIFICATE-75
## Interval-assisted fold on the second D3 daughter component

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
**PROVED_INTERVAL_ASSISTED_LOCAL_C4_FOLD**

This upgrades the numerical fold of report 70.

---

## 1. Signed reflection-fixed C4 representative

After a D12 shift, the second D3 daughter component lies in the real
reflection-fixed C4 chart.

The numerical fold center is

    s3 = -1.330174026969954
    s4 = -0.651860029115752
    s5 =  0.720774033303467
    s6 =  1.069731466438069

    g  =  4.395526393543094.

It joins the g=5 atlas branches later identified as:
- orbit i=9, index 2;
- orbit i=2, index 1.

---

## 2. Five-dimensional fold system

As in report 74, solve

    F(s,g)=0,
    det H_C4(s,g)=0,

with all four strict spectral eigenvalues retained as outward intervals.

Use a radius

    2e-8

in every fold variable.

The determinant equation is rescaled by `10^6` only for numerical
conditioning.

---

## 3. Krawczyk result

The interval Krawczyk image is strictly contained in the source box.

Representative output:

    s3 in
      [-1.330174026970594,
       -1.330174026969313],

    s4 in
      [-0.651860029116126,
       -0.651860029115376],

    s5 in
      [0.720774033303075,
       0.720774033303861],

    s6 in
      [1.069731466437567,
       1.069731466438569],

    g in
      [4.395526393542981,
       4.395526393543206].

Thus the signed-C4 fold is unique in the declared box uniformly over the
accepted strict spectral rectangle.

---

## 4. Generic/simple fold

The determinant of the augmented five-by-five Jacobian is enclosed by

    boxed:
    [0.01253616656,
     0.01253661973]

and is strictly positive.

The C4 Hessian has a one-dimensional kernel, so nonsingularity of this
augmented Jacobian gives the generic simple-fold condition.

This excludes a local cusp/higher-order fold at the certified root.

---

## 5. C4 complement

After removing the fold zero direction, interval Gershgorin bounds give three
strictly positive C4 eigenvalue enclosures:

    [0.06764006,0.06764053],

    [0.14325503,0.14325548],

    [0.15048581,0.15048622].

So the C4 block itself has

    (0 negative, 1 zero, 3 positive)

at the fold.

---

## 6. Odd/reflection-transverse block

The three odd-sector eigenvalues are enclosed by

    negative:
      [-0.12507102,
       -0.12507061],

    positive:
      [0.11452235,
       0.11452274],

    positive:
      [0.16960516,
       0.16960545].

Hence one transverse negative direction persists through the fold.

The full X7 Hessian therefore has exact fold inertia

    boxed:
    (1 negative, 1 zero, 5 positive).

---

## 7. Branch-index exchange

The two local sheets have full Morse indices

    boxed:
    2 and 1.

Thus report 70's numerical branch join is now interval-assisted:

    second D3 daughter
        index 2
          |
          | CERTIFIED fold
          | g≈4.395526393543
          v
    secondary saddle
        index 1.

At g=5 the two sheets match R7P-037 orbits `i=9` and `i=2` respectively.

---

## 8. Consequence for the D3 event

Both nonlinear D3 daughter triplets from report 66 now have certified first
fold destinations:

### `cos(3phi)=+1` triplet
    D3 crossing
      -> index-2 daughter
      -> certified upper fold g≈5.172231474684
      -> index-1 main localization saddle.

### `cos(3phi)=-1` triplet
    D3 crossing
      -> index-2 daughter
      -> certified second fold g≈4.395526393543
      -> index-1 secondary saddle.

The long branch segments remain numerical, but the local endpoints, index
exchanges and fold existence are now interval-assisted.

---

## 9. Next research direction

The next useful task is no longer finding the first folds.

It is:

    STATIONARY-BRANCH-GRAPH-COMPLETION-76

Assign the remaining g=5 atlas orbits to connected components and certify the
lowest-cost additional branch nodes, especially:
- the large negative-energy index-1 orbits;
- orbit-size-3 / stabilizer-8 branches;
- generic orbit-size-24 branches.

The target is a quotient graph of stationary D12 orbits with proof levels.
