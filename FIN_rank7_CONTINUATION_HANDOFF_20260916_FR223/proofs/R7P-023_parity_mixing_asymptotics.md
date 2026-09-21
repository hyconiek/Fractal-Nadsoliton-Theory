# R7P-023 — local parity-mixing asymptotics at the double root

On the compactified even-parity boundary (`q=1`) at the threshold balance, the
conditional covariance has a double eigenvalue `sigma`. Put `q=1-epsilon`.
Degenerate first-order perturbation on the two-dimensional threshold eigenspace
is diagonal in the natural channel basis.

For the `k3` branch the shift is

`delta_3 = lambda3 (2 t_star^2 - 1)/6`.

The accepted strict intervals certify `1-2 t_star^2>0`, hence `delta_3<0`.

For the `k4-k5` branch,

`delta_45 = - lambda4 lambda5 t_star^2 /
 [6 sqrt((lambda5-lambda4)^2+4 lambda4 lambda5 t_star^2)]`.

All factors in the quotient are strictly positive, hence `delta_45<0`.
Therefore both threshold branches move downward to first order when finite
odd-parity weight is introduced.

These are **fixed-threshold-point** first-order shifts. The imported number
`0.1312828584` for the slope of the *reoptimized* envelope is not promoted here:
that requires a separate local optimization/implicit-function argument and a
validated remainder.
