# MP7-031 — finite-copy model and variational limit

Scientific state: **PROVED_ANALYTIC, CONDITIONAL_ON_ADDED_MODEL**.

## Declared extension

Let `X=X7` be the supplied 12-by-7 feature matrix, `A=XX^T`, and let
`u0=(1/12,...,1/12)`.  For an integer `N>=1`, labels
`x_1,...,x_N in {0,...,11}`, and supplied `g>=0`, define

```
Z_N(g) = 12^{-N} sum_{x_1,...,x_N}
         exp[(g/(2N)) sum_{a,b=1}^N A[x_a,x_b]].
```

This is an added dimensionless finite-copy equilibrium model.  `N` is not
identified with a physical bit count, population, volume or time.

## Exact occupation formula

For occupations `n_j>=0`, `sum_j n_j=N`, and `p_j=n_j/N`,

```
sum_{a,b} A[x_a,x_b] = n^T A n = N^2 p^T A p.
```

There are `N!/prod_j n_j!` label sequences with those occupations, hence

```
Z_N(g) = sum_n (N!/prod_j n_j!) 12^{-N}
         exp[(Ng/2) p^T A p].
```

Because every retained nonzero Fourier column has zero uniform mean,
`X^T u0=0`, hence `A u0=0` and

```
p^T A p = (p-u0)^T A (p-u0).
```

With

```
V_g(p)=D(p||u0)-(g/2)(p-u0)^T A(p-u0),
```

the exponential rate of a type is therefore `-V_g(p)`.

## Uniform finite-state variational limit

For the uniform product law on 12 labels, the probability of a type `p` obeys
the standard finite-type inequalities

```
(N+1)^{-12} exp[-N D(p||u0)]
 <= (N!/prod n_j!) 12^{-N}
 <= exp[-N D(p||u0)].
```

These bounds can be obtained without a boundary-nonuniform Stirling formula:
the upper bound follows from the multinomial theorem after choosing the
probability vector `p`; the lower bound follows because there are at most
`(N+1)^12` types and the type `p` is a mode of the multinomial law with
parameter `p` (zero coordinates are omitted, using `0 log 0=0`).

The number of types is at most `(N+1)^12`.  Therefore, if `T_N` is the finite
set of N-types,

```
max_{p in T_N}[-V_g(p)] - 12 log(N+1)/N
 <= (1/N) log Z_N(g)
 <= max_{p in T_N}[-V_g(p)] + 12 log(N+1)/N.
```

The simplex is compact and `V_g` is continuous with the entropy convention
`0 log 0=0`; rational types are dense.  Consequently

```
lim_{N->infinity} (1/N) log Z_N(g)
  = - min_{p in Delta_11} V_g(p).
```

Equivalently, `-N^{-1} log Z_N(g)` converges to the minimum variational free
energy.

## All-pairs convention

The definition includes `a=b`.  For the supplied circulant Fourier factor,
`A[j,j]=||X_j||^2` is independent of `j` because each sine/cosine pair has
constant squared norm and the alternating column has constant square.  If
self-pairs were removed, the exponent would change by the label-independent
constant `-(g/2) A[0,0]`; this changes `Z_N` by a global factor but must not be
silently discarded.

## Scope

This theorem derives the equilibrium variational limit of the **declared
finite-copy extension**.  It does not source `N`, `g`, temperature, a physical
clock, or a microscopic dynamics.
