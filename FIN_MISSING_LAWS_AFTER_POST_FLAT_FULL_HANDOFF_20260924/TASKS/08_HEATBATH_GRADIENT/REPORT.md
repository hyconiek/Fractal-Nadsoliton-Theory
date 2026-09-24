# Heat-bath exact gradient flow of U_g

Status: **EXACT_WITHIN_DECLARED_HEATBATH**

## Result
Contracting to retained mean gives U_g(mu). Declared heat-bath mean dynamics
is exactly a gradient flow with a positive mobility equal to an average Fisher
covariance along the natural-coordinate segment.

## Key formulas
\[U_g(\mu)=I(\mu)-\frac g2\|\mu\|^2,\quad \nabla U=\theta(\mu)-g\mu\]
\[\dot\mu=-M(\mu)\nabla U,\quad M=\int_0^1\nabla^2\psi(\theta+s(g\mu-\theta))ds.\]

## Caveat
Dissipative gradient flow is not inertial mechanics.

## Next question
Exploit this lane for controlled instanton/reaction-coordinate tests.
