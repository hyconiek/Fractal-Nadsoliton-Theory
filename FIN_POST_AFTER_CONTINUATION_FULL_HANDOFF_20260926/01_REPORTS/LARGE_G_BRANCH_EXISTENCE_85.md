# LARGE-G-BRANCH-EXISTENCE-85
## Every isolated strict support generates a unique large-g stationary branch with Morse index |S|-1

Date: 2026-09-26

Status:
- asymptotic implicit-function theorem for the isolated strict support class;
- exact first correction formula;
- applies to all 76 classes of report 84.

Let S be an isolated strict support with

    A_SS p* = mu* 1,
    sum p*=1,
    p*>0,

and strict outside gap

    Delta = min_{j notin S}[mu*-(A p*)_j] >0.

Because the feature rows on S are affinely independent, A_SS is positive
definite on the tangent subspace `sum u_i=0`.

Write epsilon=1/g.

For a finite-g stationary branch,

    log p_i - log p_r
      =
      g[(A p)_i-(A p)_r].

Set

    p_S(g)=p* + u/g + O(g^-2) + O(exp[-g Delta/2]).

The first correction is uniquely determined by

    A_SS u = log p* + c 1,
    sum u=0.

Equivalently,

    [A_SS  -1][u]   [log p*]
    [1^T     0][c] = [   0   ].

Outside the support,

    p_j(g)=O(exp[-g Delta_j]),

with

    Delta_j=mu*-(A p*)_j >0.

Hence each isolated support produces one stationary branch for sufficiently
large g, unique in its support tube up to D12 images.

The dual Hessian is

    H_g = I/g - Cov_{p_g}(X).

At the limiting support the covariance has rank |S|-1 and its kernel is
unchanged by redistributing probability inside the same affine support.
Outside-support corrections are exponentially small.

Therefore, for sufficiently large g,

    boxed:
    Morse index = |S|-1.

So the 76 D12 support classes imply asymptotically:

    1 class  of index 0
    6 classes of index 1
    11 classes of index 2
    22 classes of index 3
    19 classes of index 4
    14 classes of index 5
    3 classes  of index 6.

This is an existence/classification theorem for the large-g stationary
landscape, not a statement that all these branches are already present at g=5.

Energy expansion:

    V_g(p_g)
      =
      -g mu*/2
      +D(p*||uniform)
      +(1/(2g)) u^T A_SS u
      +O(g^-2)
      +O(exp[-g Delta/2]).

Thus large-g branch ordering is first controlled by the common field mu*,
then by support entropy and the 1/g correction.
