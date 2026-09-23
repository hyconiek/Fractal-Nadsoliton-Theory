# MP7-032 — exact seven-dimensional auxiliary Gaussian representation

Scientific state: **PROVED_ANALYTIC, CONDITIONAL_ON_ADDED_MODEL**.

Assume `g>0` and use the MP7-031 finite-copy model.  Write the j-th row of
`X=X7` as `X_j` and `S=sum_a X_{x_a}`.  The elementary seven-dimensional
Gaussian identity is

```
(N/(2 pi g))^{7/2} int exp[-N||theta||^2/(2g) + theta^T S] dtheta
 = exp[g ||S||^2/(2N)].
```

It follows by completing the square; no saddle approximation is used.
Since `||S||^2=sum_{a,b} A[x_a,x_b]`, summing over labels gives

```
Z_N(g)
 = (N/(2 pi g))^{7/2}
   int_{R^7} exp[-N Phi_g(theta)] dtheta,

Phi_g(theta)=||theta||^2/(2g)
             - log[(1/12) sum_j exp(X_j^T theta)].
```

The integral is finite because the log-sum-exp grows at most linearly in
`||theta||`, while the first term is positive quadratic.

## Exact joint law and conditionals

Before summing the labels, the normalized joint density/mass is proportional to

```
exp[-N||theta||^2/(2g)] prod_{a=1}^N [exp(X_{x_a}^T theta)/12].
```

Hence, conditional on `theta`, the labels are independent with

```
P(x_a=j | theta)=softmax(X theta)_j.
```

Conditional on an empirical occupation vector `p`, completing the square gives

```
theta | p ~ Normal(g X^T p, (g/N) I_7).
```

Therefore the law of total covariance gives the exact equilibrium identity

```
Cov(theta) = (g/N) I_7 + g^2 Cov(X^T p).
```

This is exact for every integer `N>=1` and `g>0` in the declared extension.
At `g=0` the label model is defined directly; the singular Gaussian prefactor
is interpreted only through the `g downarrow 0` limit.

## Small-N normalization checks

For `N=1`, direct label summation gives

```
Z_1(g)=12^{-1} sum_j exp[g A[j,j]/2].
```

The supplied Fourier factor has constant diagonal, so this is
`exp[g A[0,0]/2]`.  Term-by-term Gaussian integration gives the identical
expression.

For any fixed `N`, expanding the Nth power of the finite log-sum factor inside
the integral and integrating each Gaussian exponential term returns exactly
`exp[g||sum_a X_{x_a}||^2/(2N)]`, so the Gaussian and label expressions agree
term by term.  MP7's numerical small-N test checks `N=1,2,3` independently at
floating precision only as an implementation diagnostic.

## Scope

The auxiliary `theta` is an exact mathematical latent variable for this added
equilibrium model.  The identity does not make it a sourced physical field and
does not determine a stochastic trajectory or physical time scale.
