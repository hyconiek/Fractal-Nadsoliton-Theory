# MP7-033 — compatibility of 11D probability and 7D mediator fluctuations

Scientific state: **PROVED_ANALYTIC for the identities; Gaussian asymptotics
CONDITIONAL_ON_STABLE_SINGLE_PHASE and on the MP7 finite-copy model**.

Let an interior probability vector `p` be fixed, let

```
Sigma = diag(p)-p p^T,
M7 = X^T Sigma X,
G7 = I_7-g M7.
```

Let `B0` be the 12-by-11 simplex tangent matrix with columns `e_i-e_11`,
`i=0,...,10`.  The probability-space Hessian in this chart is

```
H_chart = B0^T [diag(1/p)-g A] B0.
```

## Determinant identity

Set `D=diag(1/p)` and `H0=B0^T D B0`.  Direct elimination (or the matrix
determinant lemma for a diagonal matrix with one rank-one update) gives

```
det(H0) = 1/prod_j p_j,
B0 H0^{-1} B0^T = Sigma.
```

Since `A=XX^T`, the determinant lemma yields

```
det(H_chart)
 = det(H0) det[I_7-g X^T B0 H0^{-1} B0^T X]
 = det(G7) / prod_j p_j.
```

Thus

```
det(H_chart) prod_j p_j = det(G7).
```

At `g=0` both sides equal one after multiplication by `prod p_j`, fixing the
normalization.

## Leading single-phase covariance

If `H_chart>0` (equivalently `G7>0` on the retained mediator sector) and a
single interior phase is isolated, the local central-limit/Laplace covariance
on the simplex tangent is

```
Cov(p) ~ (1/N)
 [Sigma + g Sigma X G7^{-1} X^T Sigma].
```

This follows from the Woodbury identity applied to the constrained inverse of
`D-gXX^T`.  Multiplying by `X^T` and `X` gives

```
Cov(X^T p) ~ M7 (I_7-g M7)^{-1}/N.
```

`M4` cannot replace `M7` here because the declared finite-copy ensemble uses
the full rank-seven mediator.

## Exact finite-N fluctuation-response identity

Add a microscopic feature source `f` by multiplying the finite-N equilibrium
weight by `exp(N f^T mu)`, where `mu=X^T p`.  Differentiation of the finite
partition sum gives exactly

```
d E_f[mu]/d f^T = N Cov_f(mu).
```

This finite-N identity is exact and distinct from the preceding single-phase
Gaussian approximation.  It is also distinct from response to a source
conjugate directly to the auxiliary mediator `theta`.

## Scope

If `G7` is not positive, the displayed stable-phase Gaussian covariance is not
licensed.  None of these identities sources N, temperature, or dynamics.
