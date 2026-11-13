# Transverse stability of the conditional uniform A11 embedding

## Result

All eleven omitted ring directions have positive lowest squared frequencies
at the static kink throughout w in [1.25,1.75], for the frozen A11 kernel and
the explicitly imposed homogeneous coefficients of the preceding bridge.
The minimum is approximately 0.098022788 at w=1.25 in the alternating mode
j=6. This is conditional linear transverse stability of the static kink,
not nonlinear stability or stability of the radiating trajectories.

## Exact operator reduction

Use the preceding [state-map package](../fin_core_state_map_bridge_20261113/REPORT_EN.md),
with Phi=tanh(x), Psi=0, m=w-2-kappa0, y=3, c4=12 and c6=144.
Let M=(K+K^T)/2 and

    kappa_j = sum_(d=1)^11 K_sg(d) cos(2*pi*j*d/12),
    Delta_j = kappa0-kappa_j.

The Phi/Psi mixed Hessian vanishes at Psi=0. Quartic and sextic carrier
terms contribute zero to this Hessian. Each real ring mode therefore obeys

    z_tt + L_j z = 0,
    L_j = -d_x^2 + 4 - 6 sech^2(x) + w - Delta_j.

The uniform mode j=0 is retained by RP; j=1,...,11 are transverse.
Set B=d_x+2 tanh(x). On the real line,

    L_j = B* B + w - Delta_j,
    B sech^2(x) = 0.

The factorization is a quadratic-form identity on H1, with its standard
self-adjoint realization. Since sech^2 is square integrable, the bottom
of the spectrum is exactly w-Delta_j, not merely bounded below by it.
Consequently w>Delta_j is sufficient and necessary for a strictly positive
channel; w=Delta_j gives a localized zero mode, and w<Delta_j gives an
exponentially growing linear solution with rate sqrt(Delta_j-w).
These claims are exact conditional on the specified kernel eigenvalues.
The factorization and localized eigenfunction identities are symbolically
checked. A matrix eigensolver independently checks the finite Fourier formula.

## Frozen-kernel evaluation

The kernel is cos(0.18575*d+0.16250)/(1+d^1.8), d=1,...,11.
The decimal tuple is treated as the supplied mathematical input, without
physical uncertainty or a claim that it was selected by strict dynamics.

| Mode pair | Real multiplicity | Delta_j | Lowest value at w=1.25 |
|---|---|---|---|
| 1,11 | 2 | 0.381991005 | 0.868008995 |
| 2,10 | 2 | 0.763546856 | 0.486453144 |
| 3,9 | 2 | 0.968667352 | 0.281332648 |
| 4,8 | 2 | 1.079935894 | 0.170064106 |
| 5,7 | 2 | 1.135237903 | 0.114762097 |
| 6 | 1 | 1.151977212 | 0.098022788 |

Every eigenvalue increases with w with slope one, so checking the lower
endpoint classifies the entire declared interval. No transverse negative
or zero mode is found there. The overall positivity threshold is
w_critical=max_j Delta_j, approximately 1.151977212. At this threshold
the alternating carrier direction is marginal. Below it the explicit
localized instability provides a negative control for the criterion.

Numerical values are double precision, not outward-rounded interval bounds.
The analytical criterion is exact; the decimal sign classification is a
numerical evaluation with a sizeable margin, not a formal interval certificate.
results.json retains all numbers and multiplicities. The Phi sector has the
usual translational zero mode; this is not a transverse carrier instability.
No strict spectral gap for the entire unreduced system is claimed.

## What remains unresolved

The positive static transverse quadratic energy bounds the associated linear
perturbations in their energy norm. It does not prove asymptotic decay:
localized positive-frequency bound modes can persist. It does not establish
nonlinear orbital stability, a basin of attraction, or stability of the
nonstationary mixed-mode preparations already simulated. All kinetic,
coefficient-equality, state-map and physical sourcing premises remain added
assumptions. AW is not extended to the full 13-field system.

For a general uniform embedded trajectory (phi(x,t),psi(x,t)), the omitted
linear operators instead contain

    L_j(t) = -d_x^2 + w-2-Delta_j + 6 phi(x,t)^2
             + 3 psi(x,t)^2 + 5 psi(x,t)^4.

The coefficients 3 and 5 follow from Psi_i=psi/sqrt(12), c4=12 and c6=144.
The Phi cross-coupling still lies in the uniform direction, so these channels
remain decoupled at linear order in a fixed ring basis. Time dependence can
exchange energy with the background even when instantaneous operators are
positive. The static result does not exclude parametric amplification.

The next bounded step is to evolve a transverse perturbation on an already
specified radiating RP background, starting with the weakest j=6 channel
and a static-background control. Measure perturbation energy and amplification
with a step-refinement check; call it a finite-time test, not a global stability
theorem. No new same-author AW audit or parameter refit is required.

No kernel-source derivation, selector closure, legacy-role transfer, physical
units, empirical confirmation, photons, SM/GR or ToE closure follows.
