# MP7-029 — controlled critical slowing under supplied dynamics

Scientific state: **PROVED_ANALYTIC / INTERVAL-ASSISTED, CONDITIONAL_ON_ADDED_DYNAMICS**.

This task adds dynamics explicitly.  None of the rates below is a physical clock derived by FIN.

## 1. Deterministic gradient flow

Supply

`dot s = -L grad Phi(s)`

with constant symmetric positive definite mobility `L`.
At a stationary point the linearization is

`delta dot s = -L H delta s`.

As proved in MP7-028, this is similar to the symmetric matrix

`-L^(1/2) H L^(1/2)`.

### Unit mobility

For `L=I`, the stable branch's slow relaxation rate is exactly the smallest positive Hessian eigenvalue.  MP7-026 therefore gives, for every

`0 < epsilon=g-g_fold <= 1e-7`,

`0.214905205154 sqrt(epsilon) <= gamma_I <= 0.234915431887 sqrt(epsilon)`.

Hence the linear relaxation time obeys

`4.25685 / sqrt(epsilon) <= tau_I=1/gamma_I <= 4.65322 / sqrt(epsilon)`

(up to the displayed outward rounding).

This is now a controlled square-root slowing law, not a log-log regression.

### Scalar mobility demonstrates rate nonuniqueness exactly

For `L=alpha I`, `alpha>0`, the linearized generator is simply `-alpha H`. Therefore

`gamma_alpha = alpha gamma_I`,  `tau_alpha=tau_I/alpha`.

In particular `L=I` and `L=2I` have the same equilibria and the same Hessian stability classification, while every deterministic linear relaxation rate differs by a factor of two.  The square-root exponent comes from the fold potential; the prefactor additionally requires the supplied mobility.

### General constant SPD mobility

Let `v` be the normalized fold null vector.  Standard simple-eigenvalue perturbation of

`L^(1/2) H L^(1/2)`

at the fold gives the leading slow rate

`gamma_L ~ [sqrt(2 |a| b) / (v^T L^(-1) v)] sqrt(epsilon)`.

Thus even at fixed static FIN coefficients the leading prefactor is mobility-dependent.
For the illustrative diagonal choice

`L_D=diag(1,2,3,4)`,

using the R7P-031 fold vector gives

`v^T L_D^(-1) v ~= 0.538497738`,

so the nominal leading coefficient is about

`0.417662507`.

This last decimal is an asymptotic diagnostic, not a separately interval-certified finite-epsilon tube.  The fully controlled finite-epsilon statement above is the `L=alpha I` family.

## 2. Constant-mobility dual Langevin extension

On the full mediator space add the stochastic model

`d theta_t = -L grad Phi(theta_t) dt + sqrt(2 L/N) dW_t`,

where `L` is constant SPD and `N` is the supplied finite-copy parameter from MP7-031/032.

The Fokker--Planck probability current for a density `rho` is

`J = -L grad Phi rho - (L/N) grad rho`.

For

`rho_N(theta) = Z_aux^(-1) exp[-N Phi(theta)]`,

`grad rho_N = -N rho_N grad Phi`, hence `J=0` identically.  Thus `rho_N` is stationary and the process is reversible with respect to this density under the standard full-space diffusion assumptions.

Normalization is finite because the feature set is bounded and

`Phi(theta) = ||theta||^2/(2g) - log(mean_j exp(X_j theta))`

has a positive quadratic term and at most linear growth in the log-sum-exp term.  Therefore `Phi(theta)->+infinity` quadratically as `||theta||->infinity`.

By MP7-032 this auxiliary equilibrium density is exactly the Gaussian-field representation of the explicitly added finite-N label model, after the known normalization constant is included.  Equality of equilibrium measures does **not** identify microscopic trajectories.

## 3. Same equilibrium law, different kinetics

For any two constant SPD mobilities `L1` and `L2`, the Langevin models

`dtheta=-Li grad Phi dt + sqrt(2 Li/N)dW`

have the same stationary density proportional to `exp(-N Phi)` but generally different relaxation spectra and paths.
The pair `L1=I`, `L2=2I` is an exact example: the equilibrium density is unchanged while deterministic relaxation rates are doubled and the diffusion covariance per unit supplied time is also doubled.

Therefore the static potential plus its finite-N Gibbs law does not determine a unique kinetic prefactor or physical time unit.

## 4. Scope

The controlled `epsilon^(-1/2)` divergence is a conditional statement about the supplied gradient dynamics near the certified local fold.  It is not a Kramers switching-time theorem, does not identify a dominant global escape saddle, and does not source a physical mobility, noise amplitude, or clock.
