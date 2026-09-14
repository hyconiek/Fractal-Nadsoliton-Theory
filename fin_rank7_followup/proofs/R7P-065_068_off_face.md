# R7P-065--068 — off-face four-amplitude continuation

## R7P-065: exact target and stable Schur representation

For the locked four-cosine family with nonnegative fields `J3,J4,J5,J6`, let
`M4=Cov_p(C4)` and

`sigma_*=[2 lambda3(lambda4+lambda5)-lambda4 lambda5]/(24 lambda3)`.

The parity decomposition is exact:

`M4 = W_par + b b^T`,

with the fourth row/column of `W_par` zero and

`b6^2 = lambda6 q(1-q)/3`.

Define

`eta = 1-b6^2/sigma_*`.

Since `q(1-q)<=1/4` and the accepted strict intervals certify
`1-lambda6/(12 sigma_*)>0`, the scalar Schur pivot
`sigma_*-b6^2=sigma_* eta` is strictly positive everywhere. Block congruence
therefore gives

`n_-(sigma_* I4-M4) = n_-(sigma_* I3-Mtilde)`,

where

`Mtilde=W_par[345]+b[345]b[345]^T/eta`.

Thus the target is equivalently `lambda2(Mtilde)<=sigma_*`. This route never
requires the full inverse `(sigma I-W_par)^(-1)` and remains well-defined at
singular intraparity threshold points. The direct `M4` formulation is retained
as an independent numerical cross-check.

## R7P-066: adversarial numerical search

A fixed-seed differential-evolution search was run on the exact compactified
closure, not merely on a finite field box. Three seeds converge to the same
known equality point. The best reported gap `lambda2-sigma_*` is approximately
`2.8e-16`, at roundoff scale, with

`J3 ~= 0.45558363952665`, `J4,J5 ~= 0`, and `J6` large.

Structured finite-`J6` slices and one-sided `J4/J5` perturbations lie below the
ceiling. Finite points are independently replayed with the direct exponential
field implementation. This is **NUMERICAL_SEARCH_ONLY**: saturation of the
search does not prove the global ceiling.

## R7P-067: exact global compactification

Set

`x_k=exp(-J_k)`, `k=3,4,5,6`.

Multiplying every Boltzmann weight by the same factor
`exp(-(J3+J4+J5+J6))` gives

`w_j = product_k x_k^(1-cos(2 pi k j/12))`.

All exponents are nonnegative. For `j=0` every exponent is zero, so `w_0=1`
on the entire closed cube. Consequently the normalized law extends
continuously to `[0,1]^4`; every unbounded field sequence has a convergent
compact-coordinate subsequence, and every closed-cube point is approximated by
finite fields. This is exactly the closure of the nonnegative four-field
family, including simultaneous and multiscale large-field limits.

For compatibility with the already-certified boundary proof use

`r=exp(-2J3), s=exp(-3J4/2), t=exp(-J5/2), y=exp(-2J6)`.

After aggregating reflection-identical cosine states, the even weights are

`1, 2 s t^3, r t^4, 2 r s t`,

and the odd weights are

`2 sqrt(r) s t^(2+sqrt(3)) y`,
`2 sqrt(r) t^2 y`,
`2 sqrt(r) s t^(2-sqrt(3)) y`.

At `y=0` this is exactly the R7P-048--055 boundary cube. The same `t` occurs in
both irrational powers, so no forbidden independent polynomialization is used.
There is no discarded-probability tail error.

## R7P-068: first-order cone result and remaining atom

At the compactified equality point `(r_*,1,1,0)`, R7P-052/R7P-055 control the
entire physical boundary `y=0`. R7P-023 gives the two degenerate threshold
shifts under parity mixing `epsilon=1-q`; both are strictly negative for `M4`.
Hence the projected `K=sigma_* I-M4` mixing derivative is positive definite.

For a physical tangent, write the projected first-order matrix as boundary
part plus the nonnegative parity-mixing part. Boundary safety implies the
second ordered covariance derivative cannot be positive at the equality point.
Adding a strictly negative-definite covariance perturbation from any first-order
parity mixing makes that derivative strictly negative (Weyl inequality on the
two-dimensional threshold eigenspace). Thus no first-order physical tangent can
create two supercritical directions.

The formerly imported reoptimized extreme-face slope is now interval-certified.
Allowing `J3=J3_*+d epsilon` and optimizing the lower of the two threshold
branches gives

`d ~= -0.3362414530085`,

`d lambda2_envelope / d epsilon ~= -0.1312828584000`,

with strict rational interval bounds stored in
`certificates/R7P-023_reoptimized_envelope.json`.

This is not yet the acceptance-grade R7P-068 neighborhood certificate. The
remaining atom is a validated mixed second/higher-order remainder bound that
turns the directional result into an **explicit** 4D radius, especially for
paths whose parity weight enters at higher order than the boundary displacement.
No global R7P-071 theorem is claimed here.
