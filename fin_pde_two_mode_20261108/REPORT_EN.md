# Shape and mixed radiation channels: auxiliary PDE experiments

Date: 2026-11-08. Parent 5c923c2f50437599971a90c25c23fce20364f0d1.

Six nonlinear Fourier/RK4 runs test two explicitly different physical
preparations at w=1.5 and g4=g6=1. Shape-only starts eta=epsilon g,
psi=0. Mixed starts eta=epsilon g, psi=epsilon f. All physical velocities
and independent continuum data start at zero. These preparations are
outside the historical AW prepared-data theorem; no extension of that
theorem is claimed. The full PDE and canonical reconstruction are retained.

Both preparations are run at epsilon=0.15 and 0.1 on dx=0.1,dt=0.02;
epsilon=0.1 is repeated on dx=0.05,dt=0.01. All use T=100,L=130 and
physical flux controls at R=10,20. Consequently the two amplitudes do
not have identical slow-time endpoints. This is a finite-window experiment,
not an asymptotic convergence theorem.

## Estimators and results

The unchanged CD coefficients at w=1.5 are
k_a=0.2836716043, k_b=0.0347670893, k_x=0.2588921676.
Shape-only preserves psi=0 in the solver and has I=0. Its effective k_b
is the least-squares slope of 1/J(t).

For mixed data, define X(t)=integral I J and subtract the fixed predicted
self terms from the action losses:

    Y_I=I(t0)-I(t)-k_a integral I^2,
    Y_J=J(t0)-J(t)-k_b integral J^2.

Fit each Y against X with a free intercept on t in [25,100]. This avoids
differentiating oscillatory action data. It is conditional on the supplied
self coefficients, not an unbiased simultaneous identification of all
three parameters. The fits are diagnostics, never inputs to the PDE.

| Preparation | epsilon | dx | Effective coefficient(s), window [25,100] |
| --- | --- | --- | --- |
| Shape | 0.15 | 0.10 | k_b=0.0351085195 |
| Shape | 0.10 | 0.10 | k_b=0.0349138105 |
| Shape | 0.10 | 0.05 | k_b=0.0349136872 |
| Mixed | 0.15 | 0.10 | k_x from I=0.266470304; from J=0.266563637 |
| Mixed | 0.10 | 0.10 | k_x from I=0.262822226; from J=0.262746814 |
| Mixed | 0.10 | 0.05 | k_x from I=0.262822278; from J=0.262746678 |

At epsilon=0.1 on the finer mesh, k_b is about 0.422% above theory;
k_x is about 1.518% and 1.489% above theory from the two action equations.
Both deviations decrease between the tested amplitudes, but remain clearly
larger than the mesh-refinement changes. The positive mixed contribution
and agreement of the two estimates support the channel structure while
leaving its quantitative finite-amplitude correction unresolved.

Changing the analysis window to [50,100] gives k_b=0.0348908796 and
k_x=0.2628979685 / 0.2628458820 on the fine epsilon=0.1 runs. No preferred
window was selected to erase the discrepancy. Quadrature/sample sensitivity
also matters: using every second stored sample changes the mixed estimates
to 0.2629710271 / 0.2627843521 on [25,100]. This is not a rigorous error
bar, and the fitting samples are temporally correlated. No statistical
confidence interval or fitted asymptotic exponent is claimed.

## Numerical checks and limits

On the finer mesh, relative total-energy drift is 7.50e-9 (shape) and
1.31e-8 (mixed). The maximum sampled local energy-flux residual divided
by initial energy is 4.07e-8 (shape) and 2.47e-8 (mixed), rounded upward.
Physical background and interaction terms are retained. The fine residual
need not decrease by the time-integrator order because local quadrature,
spatial approximation and shared diagnostics also contribute.

These computations share RK4 and the canonical map with earlier runs.
They are not independent-team replication, empirical FIN validation or
proof of AW for new data. Fourier collocation has no explicit dealiasing;
resolution refinement is a check, not a certified aliasing bound.

## Next step

Resolve the mixed-channel discrepancy with a smaller-amplitude experiment,
matched slow-time windows and varied relative modal phase, keeping the CD
coefficients fixed. Accumulate action-product integrals at solver stages
or finer diagnostic spacing to separate quadrature effects. Compare the
full two-action ODE trajectories and report reconstruction, transient and
finite-amplitude contributions separately. Do not refit the theoretical
coefficients merely to force agreement. The FIN-to-action nonlinear
matching problem remains separately open; routine self-audit chaining
remains stopped.
