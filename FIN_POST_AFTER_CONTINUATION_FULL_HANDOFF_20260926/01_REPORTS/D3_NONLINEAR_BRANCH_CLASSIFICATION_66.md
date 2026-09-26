# D3-NONLINEAR-BRANCH-CLASSIFICATION-66
## The MP7-039 crossing is a genuinely cubic D3 bifurcation, not a hidden pitchfork

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Inputs:
- `fin_rank7_mathphysics_next/proofs/MP7-039_boundary_counterexample_mechanism.md`;
- `fin_rank7_mathphysics_next/results/MP7-039_transverse_crossing.json`;
- accepted strict spectral intervals from
  `R7P-036_stationary_index2_witness.json`.

Status:
- exact D3 invariant-theory reduction;
- interval-separated nonzero cubic coefficient;
- interval-separated crossing slope and complement inertia;
- local nonlinear daughter-branch classification;
- no global phase/minimizer claim.

---

## 1. Certified linear crossing

MP7-039 supplies a unique crossing near

    J* = 0.038434753979
    K* = 0.185272061402
    g* = 5.171841831943

on the exact two-harmonic stationary family.

The critical space is two-dimensional and carries the standard real D3
representation.

Write an orthonormal critical coordinate as

    z = x + i y
      = r exp(i phi).

The base state is D3-invariant.

---

## 2. Critical eigenvector

The cosine `(k4,k5)` Hessian block is

    H45 =
      [[A,C],
       [C,D]]

with zero determinant at the crossing.

Using the accepted crossing/spectral boxes gives a zero-direction of the form

    e1 =
      a e_(4c) + b e_(5c),

and the D3-related sine direction

    e2 =
      a e_(4s) - b e_(5s),

with interval-enclosed coefficients approximately

    a in [0.3899593, 0.3900531],
    b in [0.9207976, 0.9208269].

Thus `(e1,e2)` is the real critical D3 plane.

---

## 3. Cubic coefficient is not zero

For

    Phi(theta,g)
      = ||theta||^2/(2g)
        - log[(1/12) sum_j exp((X theta)_j)],

the third derivative is minus the third cumulant of the corresponding feature
direction.

Interval evaluation over:
- the full MP7-039 crossing box;
- the accepted lambda4,lambda5 intervals;
- the induced critical eigenvector interval

gives

    boxed:
    T3(e1,e1,e1)
      in
      [-0.015475047964,
       -0.015413822050].

Hence the D3-allowed cubic term is rigorously nonzero.

The symmetry-related tensor pattern is

    T111 = t,
    T122 = -t,
    T112 = T222 = 0,

up to the chosen orientation of `(e1,e2)`.

Therefore the reduced potential begins as

    boxed:
    Phi_red
      =
      Phi0
      +(lambda/2) r^2
      +(t/6) r^3 cos(3 phi)
      +O(r^4, |mu| r^3),

where

    mu=g-g*.

This alone excludes a pitchfork normal form.

---

## 4. Crossing slope

MP7-039 certifies

    d(det H45)/dg
      in
      [0.001695735102,
       0.001696038770].

The nonzero eigenvalue of H45 at the crossing is its positive trace, enclosed by

    trace(H45)
      in
      [0.011860844947,
       0.011862158823].

Since

    det = lambda_soft lambda_hard,

and lambda_soft=0 at the crossing,

    lambda_soft'(g*)
      =
      [d(det)/dg]/lambda_hard.

Therefore

    boxed:
    alpha := lambda_soft'(g*)
      in
      [0.1429533298,
       0.1429947679]
      >0.

Thus

    lambda(mu)=alpha mu+O(mu^2).

---

## 5. Daughter branches

Stationarity of the leading D3 normal form gives

    sin(3 phi)=0

and, for r != 0,

    lambda
      +(t/2) r cos(3 phi)
      +O(r^2)=0.

Hence

    boxed:
    r
      =
      -2 lambda/
       [t cos(3 phi)]
      +O(lambda^2).

Because t<0:

### g<g*
`lambda<0`, so the physical positive-r branches use

    cos(3 phi)=-1.

There are three D3-related rays.

### g>g*
`lambda>0`, so the physical positive-r branches use

    cos(3 phi)=+1.

Again there are three D3-related rays, rotated by pi/3 relative to the
below-crossing set.

Thus the crossing exchanges two inequivalent triplets of daughter directions.

---

## 6. Linear amplitude law

The leading branch amplitude satisfies

    r
      ~ C_r |g-g*|,

where

    C_r=2 alpha/|t|.

Using the interval bounds above,

    boxed:
    18.4753
      < C_r <
    18.5541.

So this D3 event has **linear**, not square-root, daughter amplitude scaling.

That is another direct distinction from the simple fold/pitchfork intuition.

---

## 7. Critical-plane Morse index of daughters

At a daughter branch the two critical-plane Hessian eigenvalues are, to leading
order,

    lambda_radial
      = -lambda + O(lambda^2),

    lambda_angular
      = 3 lambda + O(lambda^2).

They always have opposite signs for sufficiently small nonzero mu.

Therefore every sufficiently local daughter branch has exactly **one negative
direction inside the critical D3 plane**.

---

## 8. Noncritical complement at the crossing

Interval evaluation gives:

    H_(3sin) > 0,

while the `(3cos,6)` block has

    det H36
      in
      [-2.16354e-5,
       -2.16269e-5]
      <0.

Hence the noncritical five-dimensional complement contains exactly one negative
direction at the crossing.

The positive H45 hard eigenvalue remains separated from zero.

Therefore, sufficiently close to g*:

### Base two-harmonic branch
- for g<g*: two additional negative critical directions;
  full index = 3;

- for g>g*: both critical directions become positive;
  full index = 1.

### Each daughter branch
- one negative complement direction;
- one negative critical-plane direction;

so the full daughter index is

    boxed:
    index = 2

on both sides of the crossing.

Thus the local Morse pattern is

    index-3 base
        -> D3 crossing
        -> index-1 base,

with index-2 daughter saddles attached on both sides.

---

## 9. Energy splitting

Substitute the leading daughter amplitude into the reduced potential:

    Delta Phi_daughter-base
      =
      (2/3)
      lambda^3/t^2
      +O(lambda^4).

Since

    lambda~alpha mu,

    boxed:
    Delta Phi
      ~ C_E (g-g*)^3,

with nominal

    C_E approximately 8.1684.

Therefore:
- below g*, daughter branches lie below the base branch locally;
- above g*, they lie above it locally.

This is a cubic energy exchange, not the fold's `|g-g_f|^(3/2)` law.

---

## 10. Numerical branch replay

Direct full-seven-coordinate solves at `|g-g*|=10^-5` give critical amplitudes

    g<g*: r/|mu| approximately 18.40
    g>g*: r/|mu| approximately 18.63,

approaching the predicted interval around 18.5 as `mu->0`.

The corresponding energy difference divided by `mu^3` is of the predicted
order-eight magnitude.

This replay is numerical validation of the analytic/interval local
classification, not a replacement for it.

---

## 11. Physical/emergence interpretation boundary

Mathematically, one D3-symmetric relational state loses/gains two transverse
directions and emits three symmetry-related daughter directions.

This is a concrete symmetry-breaking mechanism.

It does NOT identify:
- particles;
- Standard Model families;
- physical spatial dimensions;
- a laboratory bifurcation;
- a physical value of g.

The result is structural only.
