# FR2 — strengthened determinant reserve on the full extreme face

This continuation strengthens the earlier face reserve `9/50` without changing the completed R7P ledger.

On `J4=J5=0`, write `u=x^2` with `x=q tanh(J3)`, `a=lambda3/6`, and
`eta=1-lambda6 q(1-q)/(3 sigma)`.  The Schur-reduced rank-one update changes only
the k3 scalar.  In the branch where the (4,5) block is the unique supercritical
block and the k3 scalar is subcritical,

`R = det(sigma I-Mtilde)/det(sigma I-Wpar) = (sigma-f)/(sigma-w)`.

The physical `J6>=0` constraint is exactly `q in [(1+u)/2,1]`.  The dangerous
(4,5) block requires `u>=T=t_*^2=1-sigma/a`.

## Monotonicity in q

Parameterize

`u=T+(1-T)z`, `q=(1+u+(1-u)v)/2`, `z,v in [0,1]`.

After removing manifestly positive factors from `dR/dq`, the remaining numerator
is a degree-(4,4) polynomial.  Exact-rational interval Bernstein subdivision
certifies it nonnegative on the square except the single degeneracy corner
`z in [0,1/4096], v in [4095/4096,1]`.  With `W=1-v` the exact corner
polynomial factors as

`N = A z + W^2 Q(z,W)`.

The strict spectral intervals give `A>0`, and a second Bernstein calculation
certifies `Q>0` on `z,W in [0,1/4096]`.  Hence `dR/dq>=0` throughout the
physical dangerous branch wherever the ratio is defined.  The minimum is at
`q=(1+u)/2`, i.e. `J6=0`.

## Uniform reserve

At that endpoint, subtract `19/100` from `R` and substitute
`u=T+(1-T)z`.  The denominator is strictly negative on `z in [0,1]`.
The numerator is certified nonpositive by a six-leaf one-dimensional Bernstein
cover.  Therefore

`R >= 19/100`

on the entire dangerous face.  If instead the k3 scalar itself is already
supercritical, the rank-one update acts in that same scalar channel and cannot
create a second supercritical direction; the threshold-singular scalar case is
handled directly by the block structure.

This is a face theorem only.  It does not prove monotonicity under positive J4
or J5 perturbations.
