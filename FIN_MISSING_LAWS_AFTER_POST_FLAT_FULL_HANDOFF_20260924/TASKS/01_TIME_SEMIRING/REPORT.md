# TIME-SEMIRING-01 — sum-product to min-plus

Status: **CONDITIONAL_THEOREM**

## Result
For a finite-N kernel composed by ordinary probabilistic summation, applying
`C_N=-(1/N)log K_N` yields min-plus/dynamic-programming composition in the
large-N Laplace-principle limit. For 12-state N-types the soft-min error is
bounded by `12 log(N+1)/N`.

## Key formulas
\[K_N(x,z)=\sum_yK_1(x,y)K_2(y,z),\qquad C_N=-N^{-1}\log K_N\]
\[C(x,z)=\min_y[C_1(x,y)+C_2(y,z)].\]

## Caveat
Does not source the transition kernel/rates.

## Next question
Source or constrain the microscopic activity.
