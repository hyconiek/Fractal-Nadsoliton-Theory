# Leading Gaussian CLT after hidden elimination

Status: **NEGATIVE_FOR_FRICTION_MEMORY**

## Result
Linearized drift is triangular: hidden fluctuations do not feed back into
visible drift. Therefore leading visible Gaussian CLT remains time-local. Hidden
history survives through the trajectory-dependent covariance and at path-LDP
scale, not as a Gaussian friction kernel.

## Key formulas
\[\frac d{dt}\binom uv=\begin{pmatrix}A(t)&0\\C(t)&-I\end{pmatrix}\binom uv+\text{noise}.\]

## Caveat
Do not claim a generalized-Langevin friction from leading CLT.

## Next question
Inspect the first non-Gaussian correction.
