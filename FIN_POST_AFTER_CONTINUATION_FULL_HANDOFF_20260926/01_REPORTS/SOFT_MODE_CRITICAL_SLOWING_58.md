# SOFT-MODE-CRITICAL-SLOWING-58
## Static softening forces dynamic slowing in the declared maximum-entropy heat-bath lane

Date: 2026-09-26

Status:
- exact local theorem conditional on the heat-bath gradient structure;
- no physical time scale is inferred.

## 1. Linearized operator

At an interior stationary state,

    B=-S H,

where

    S=diag(p)-pp^T

is positive definite on the simplex tangent space and

    H=Hess_tangent V_g.

## 2. Reality of the relaxation spectrum

On the tangent space B is similar to the symmetric matrix

    -S^(1/2) H S^(1/2):

    S^(-1/2) B S^(1/2)
      = -S^(1/2) H S^(1/2).

Therefore all linearized heat-bath relaxation eigenvalues are real.

If H>0, all nontrivial eigenvalues of B are strictly negative.

## 3. Fold/soft-mode theorem

If H has a one-dimensional kernel spanned by v at a simple static fold, then

    S^(1/2) H S^(1/2)

also has a one-dimensional kernel.

Therefore B has exactly one zero relaxation eigenvalue.

As the smallest positive curvature tends to zero, the corresponding relaxation
rate tends to zero.

Thus:

    static soft mode
      -> critical slowing

is exact in this kinetic lane.

## 4. Bounds

Let

    s_min <= eigenvalues(S|_T) <= s_max

and

    h_min <= eigenvalues(H|_T) <= h_max.

Then the positive decay rates of -B obey the Rayleigh bounds

    s_min h_min
      <= rate_min
      <= s_max h_min,

and similarly for the largest rate.

So vanishing h_min necessarily forces rate_min -> 0 as long as S remains
nonsingular on the tangent space.

## 5. Fluctuation counterpart

From report 57,

    Cov ~ H^{-1}

at Gaussian order.

Therefore approaching a soft mode produces simultaneously:

    relaxation time -> infinity
    Gaussian variance along the soft direction -> infinity

within the local approximation.

This is the standard static/dynamic critical pairing, now derived internally
for the declared FIN heat-bath process.

## 6. Boundary

The repository already has certified fold/localized stationary structures, but
this report does not assign laboratory seconds or claim a real physical critical
phenomenon.

It proves the implication only if the maximum-entropy heat-bath dynamics is
adopted as the dynamics of those states.
