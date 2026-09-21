# FR3 — a certified off-face tube around the entire extreme face

Let `M4(J3,J4,J5,J6)` be the physical four-amplitude covariance and let
`rho=1/8192` be the already-certified R7P-068 local-cone radius.

## 1. Uniform face gap away from a cone-contained box

On `J4=J5=0`, put `t=tanh(J3)`, `q=P(even)`, and `y=2q-1`.  Physical
`J6>=0` is exactly `q^2 t^2 <= 2q-1`, i.e. `q^2 t^2<=y`.
The full 4x4 covariance splits into the `(k3,k6)` and `(k4,k5)` 2x2 blocks.
For `tau=sigma-1/150000`, the smaller eigenvalue of each block is uniformly
below `tau`.  Hence `lambda2<=tau` follows from the disjunction that at least
one of the two shifted block determinants is nonnegative.

Exact-rational interval Bernstein subdivision on `(y,t) in [0,1]^2`, with the
nonphysical region excluded by `y-q^2 t^2<0`, proves this disjunction outside

`y in [8191/8192,1]`,
`t in [13973/32768,27953/65536]`.

The cover has zero unresolved leaves at depth 32.  The excluded box is contained
in the face projection of the R7P-068 cone: `1-q<=1/16384`, and the exact
`t_*` interval plus `|d[(1-t)/(1+t)]/dt|<=2` give `|r-r_*|<rho`.

Thus outside the box

`lambda2(M4(J3,0,0,J6)) <= sigma - 1/150000`.

## 2. Global covariance Lipschitz bound

For an exponential-family field parameter `J_k`,

`d Cov(X)/dJ_k = E[(X-mu)(X-mu)^T (c_k-E c_k)]`.

Therefore

`||d Cov/dJ_k|| <= range(c_k) tr Cov`.

Popoviciu on each scaled coordinate gives the global bound

`tr Cov <= lambda3/6 + 3 lambda4/32 + lambda5/6 + lambda6/12 =: T_max`.

Consequently one may take

`L4=(3/2)T_max`, `L5=2T_max`.

If `L4^+ J4 + L5^+ J5 <= 1/150000`, Weyl moves every covariance eigenvalue by
at most the certified face gap.

## 3. The excluded local box after perturbation

Changing `J4,J5` also changes `q`.  Since

`dq/dJk = Cov(1_even,c_k)`, Popoviciu/Cauchy gives
`|dq/dJ4|<=3/8`, `|dq/dJ5|<=1/2`.

Under the same tube condition this drift is much smaller than the half-radius
left between `1-q<=1/16384` and the R7P-068 limit `1-q<=1/8192`.
Also

`1-exp(-3J4/2)<=3J4/2`, `1-exp(-J5/2)<=J5/2`,

and both remain below `rho`.  Thus points above the excluded face box remain
inside the R7P-068 cone.

Combining the two cases yields a certified four-amplitude tube for all
`J3,J6>=0`.  A simple corollary is

`J4,J5 <= 1/600000`.

This result does not transfer to the full seven-coordinate Hessian.
