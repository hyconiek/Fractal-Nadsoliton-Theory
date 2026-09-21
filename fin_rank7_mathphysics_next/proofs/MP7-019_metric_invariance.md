# MP7-019 — metric-covariant form of the covariance-curvature theorem

Let the original mediator coordinates have Euclidean quadratic metric and covariance `M`. Under an invertible linear change

`theta = T eta`,

the feature matrix becomes `C'=C T`, the quadratic metric becomes

`G'=T^T T`,

and the covariance becomes

`M'=T^T M T`.

Accordingly the Hessian transforms by congruence:

`H' = G'/g-M' = T^T(I/g-M)T`.

Sylvester's law of inertia therefore makes the Hessian negative/zero/positive counts invariant under every invertible linear reparameterization.

The normalization-invariant spectral statement is expressed through generalized eigenvalues

`M' v = lambda G' v`.

Putting `w=T v` reduces this equation exactly to `M w=lambda w`. Thus the generalized spectrum of `(M',G')` is the original covariance spectrum. Raw eigenvalues of `M'` alone are not invariant when `G'` is discarded.

## Exact rescaling regression

Take `M=diag(1/4,1/8)` and `T=diag(2,1)`. Then

`M'=diag(1,1/8)` and `G'=diag(4,1)`.

The ordinary covariance eigenvalues changed from `{1/4,1/8}` to `{1,1/8}`, while the generalized eigenvalues remain exactly `{1/4,1/8}`. The Hessian inertia computed from `G'/g-M'` is therefore the same as from `I/g-M`.

For a nonlinear chart, the Hessian contains additional second-derivative-of-chart terms proportional to the gradient. Only at a stationary point do those terms vanish and the Hessian reduce to a Jacobian congruence. This is the same issue paid explicitly in MP7-037's polar-coordinate identity.

Therefore `67/250` is meaningful only together with the supplied feature normalization/metric; it is not a normalization-free physical constant.
