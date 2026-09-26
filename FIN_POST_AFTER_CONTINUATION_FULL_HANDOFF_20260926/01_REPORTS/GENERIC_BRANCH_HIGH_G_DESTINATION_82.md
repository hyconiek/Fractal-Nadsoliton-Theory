# GENERIC-BRANCH-HIGH-G-DESTINATION-82
## The i=5 generic branch is asymptotically a three-label saddle family

Date: 2026-09-26

Status:
- long-range continuation to g=200 is numerical;
- the zero-temperature support equations are solved exactly as a finite linear system;
- asymptotic index statement follows from the support covariance rank, conditional on convergence to that support.

## 1. High-g continuation

Start from the trivial-stabilizer g=5 atlas orbit i=5, index 2.

Fixed-g continuation, cross-checked against the pseudo-arclength branch, was
carried to

    g=200

without an index change.

The branch remains generic: no D12 reflection or rotation stabilizer is
restored.

## 2. Concentrating support

As g increases, the categorical distribution p concentrates on three labels

    S={1,4,8}

using zero-based label numbering.

Representative values are:

    g=20:
      p_S approximately
      (0.34536, 0.24731, 0.35054)

    g=50:
      already within about 5.6e-3 in L1 of the limiting three-support law

    g=100:
      L1 error about 2.70e-3

    g=200:
      L1 error about 1.33e-3.

All other label probabilities are exponentially suppressed.

## 3. Zero-temperature support equations

Stationarity is equivalent to

    p = softmax(g A7 p).

If p_g converges to a finite support S as g->infinity, all labels remaining in
S must have equal leading field:

    (A7 p*)_i = m
    for i in S.

Together with

    sum_{i in S} p_i*=1,

this gives a linear system.

For S={1,4,8} the unique solution is

    boxed:
    p_1* = 0.352590211929591
    p_4* = 0.289510622664754
    p_8* = 0.357899165405655.

The common field value is

    m≈0.481327320872340.

Every outside label lies strictly below this value.

The smallest support-to-outside field gap is

    boxed:
    Delta_gap≈0.477033036060931.

So the three-label support is strongly self-consistent at zero temperature /
large g.

## 4. Limiting retained mean

The limiting retained moment is

    mu*=X^T p*

with

    ||mu*||≈0.693777573054.

Numerically the continued branch satisfies

    ||theta||/g -> ||mu*||

with high accuracy.

## 5. Asymptotic Morse index

The full dual Hessian is

    H_g = I/g - X^T S(p_g) X.

For a distribution supported on three affinely independent feature vectors,
the categorical covariance has rank

    3-1=2.

Thus, as g->infinity:
- two eigenvalues of H_g approach strictly negative limits;
- the remaining five approach zero from the positive side as 1/g.

Therefore the asymptotic Morse index is

    boxed:
    index = 2.

This is exactly the index observed numerically along the full i=5 high-g
continuation.

## 6. Interpretation

The generic branch is not a finite symmetry-breaking loop that must close back
onto a reflection branch.

Its high-g destination is an intrinsically asymmetric three-label mixture.

That explains why the branch can retain trivial stabilizer indefinitely.

## 7. Branch role

The full known genealogy is now:

    certified reflection-breaking crossing
      g≈5.152672504944
        ->
    generic index-3 sheet
        ->
    atlas i=14 at g=5
        ->
    generic fold g≈4.913081079101
        ->
    generic index-2 sheet
        ->
    atlas i=5 at g=5
        ->
    high-g three-label saddle
        with index 2.

No additional symmetry-restoring event is required.

## 8. Next question

The same support analysis should be applied systematically to every large-g
branch.

That motivates `HIGH-G-SUPPORT-MORSE-CLASSIFICATION-83`.
