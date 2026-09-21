# R7P-033 — registered full-7D nonstationary curvature counterexample

Status: **EXACT_PROVED / regression witness**, inherited from the accepted audit and independently wired into this follow-up package.

For
\[
h_j=2\cos(\pi j/2),\qquad p=\operatorname{softmax}(h),
\]
we have `x=tanh(1)>3/4`. The full seven-coordinate covariance contains two orthogonal `(4,5)` trial directions (one cosine, one sine), each with Rayleigh quotient
\[
R={\lambda_4+\lambda_5+2\sqrt{\lambda_4\lambda_5}x\over24}>{313\over960}.
\]
The strict interval enclosures pay `lambda4>2.19`, `lambda5>2.29`, `sqrt(lambda4 lambda5)>2.23`, while six positive Taylor terms give `exp(2)>7` and hence `tanh(1)>3/4`. For every `g>=3.7`,
\[
1/g\le10/37<313/960,
\]
so `H7=I/g-Cov(X7)` is negative definite on this two-dimensional trial subspace.

Therefore `H7` has at least two negative directions at this finite field. The field is **not stationary**, so this refutes only an everywhere-curvature statement, not a stationary-point-only statement.

R7P-036 separately supplies a genuinely stationary index-two witness at exact `g=5`.
