# Fourier replication and physical local energy flux

Six full nonlinear Fourier-collocation evolutions test w=1.25,4/3,1.75
at epsilon=0.1, g4=g6=1, T=100, L=130. Resolutions (dx,dt)=(0.1,0.02)
and (0.05,0.01) are compared with the previous fourth-order finite-difference
campaign. The spatial method is different; RK4 and the canonical observable
reconstruction are shared. This is not independent-team replication.

## Damping comparison

| w | Fourier dx | Effective k_a | Difference from FD dx=0.1 | Difference from theory |
| --- | --- | --- | --- | --- |
| 1.25000000 | 0.10 | 0.408588859 | -0.00386% | -0.2184% |
| 1.25000000 | 0.05 | 0.408588854 | -0.00386% | -0.2184% |
| 1.33333333 | 0.10 | 0.364456503 | -0.00434% | -0.2458% |
| 1.33333333 | 0.05 | 0.364456505 | -0.00434% | -0.2458% |
| 1.75000000 | 0.10 | 0.192256518 | -0.00573% | -0.4302% |
| 1.75000000 | 0.05 | 0.192256516 | -0.00573% | -0.4302% |

The same slope of 1/I over t in [25,100] is used, with the same sampled
times. The damping law is not inserted into either field solver.

## Physical balance

At R=10 and R=20 compute the actual flux
j=-phi_t phi_x-psi_t psi_x, with phi=K+eta. Integrate
j(R)-j(-R) at every RK4 stage with its RK4 weights. Compute the energy
inside [-R,R] independently by composite Simpson quadrature, subtracting
the static kink density. The full potential and the linear terms
K' eta_x+2K(K^2-1)eta are retained. The tested residual is

    E_R(t)-E_R(0)+integral_0^t [j(R)-j(-R)] ds.

It is normalized by the initial total perturbation energy. The stored
maximum is over diagnostic times spaced by 0.5, not a continuous-time
supremum. All six runs show positive final net outward flux and all
reported balance residuals are below 1e-7. The fine-grid numbers and
full sampled budgets are in results.json. Flux was not defined as energy
loss; it is separately evaluated from physical derivatives.

## Limits and next work

The physical balance and cross-discretization agreement support genuine
energy transport rather than numerical dissipation in these cases. They
do not prove AW, identify all three damping coefficients, or establish
empirical FIN. The solver has no explicit dealiasing: the doubled spatial
resolution is a numerical check, not a rigorous aliasing bound. No absorber
is used, and the box is larger than T plus a margin. New domain checks and
a different temporal integrator remain possible additional controls.

Next use separately labelled two-mode initial preparations to identify
k_b and k_x and compare their fixed theoretical values, retaining energy
and flux checks. Such runs test the supplied model beyond the historical
AW preparation and are not automatically covered by that theorem. Continue
the FIN nonlinear matching problem separately; numerical agreement does
not derive this PDE from the fundamental core. Routine self-audit chaining
remains stopped and independent review of AW remains open.
