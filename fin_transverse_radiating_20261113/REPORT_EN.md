# Transverse perturbations on a radiating RP background

This bounded experiment evolves the weakest omitted ring channel j=6 on
the conditional uniform A11 embedding. It tests finite-time amplification;
it does not establish nonlinear or all-time stability.

## Setup

Use the frozen A11 kernel shift Delta_6 from the preceding static-kink
package. Fix w=1.5, epsilon=0.1, T=100, periodic domain [-130,130).
The driven background is the full nonlinear RP solution with initial
eta=epsilon sqrt(3/2) sech(x)tanh(x), psi=epsilon sqrt(3/4)sech^2(x),
and zero velocities. The static control uses eta=psi=0.
The background equations are unchanged from the previous PDE campaigns.

Each infinitesimal transverse perturbation solves

    z_tt = z_xx - [w-2-Delta_6 + 6 phi^2 + 3 psi^2 + 5 psi^4] z,
    phi = tanh(x)+eta.

Two columns are evolved simultaneously: (z,z_t)=(f/sqrt(gap),0) and
(0,f), f=sqrt(3/4)sech^2(x), gap=w-Delta_6. Their normalization is a basis
choice in a linear variational problem, not a finite physical perturbation.
They have no back reaction on the background. All RK4 stages use the
corresponding evolving background, avoiding frozen-step interpolation.
Fourier collocation is used without explicit dealiasing.

The fixed positive reference energy is

    E0[z,v] = integral(v^2 + z_x^2 + [4-6 sech^2(x)+gap]z^2)/2.

It is conserved on the static background but need not be conserved on the
driven one. Using this fixed norm avoids treating an instantaneous possibly
varying quadratic form as an invariant. A two-by-two Gram matrix measures
both columns and their cross term. Its generalized eigenvalues relative to
the initial Gram matrix give minimum and maximum energy gains over ALL
initial phase combinations within this two-dimensional localized seed space.
This is not a supremum over arbitrary spatial perturbations. Energy gain
and amplitude gain differ: amplitude in this norm is the square root.

## Controls and interpretation

Four runs comprise a static control at dx=0.1, dt=0.02 and driven cases
(dx,dt)=(0.1,0.02), (0.1,0.01), (0.05,0.01). Thus time refinement is
separated from spatial refinement. Diagnostics are stored every 0.5 time
unit; maxima refer to those samples, not a certified continuous-time maximum.
No extrapolated continuum error bound is claimed. Raw CSV columns preserve
both optimal gains and gains of the two individual seeds.

The numerical summary is in RESULTS.md and results.json. A gain greater
than one indicates finite-time transfer into the transverse perturbation
norm. It is not by itself evidence of exponential instability; bounded
modulation or transient growth can also produce it. Conversely a modest
gain over T=100 cannot establish long-time stability, other preparations,
other ring channels or nonlinear stability of the 13-field model.

This study inherits the conditional homogeneous coefficients and wave
kinetics of the preceding bridge. It does not derive them from strict FIN,
extend AW, resolve the selector, or supply physical units or experimental
evidence. No historical scientific inputs are edited.

## Next bounded step

Compare amplification across initial background phase and w, especially
near w=1.25 where the static transverse gap is smallest. Preserve the
same seed-space norm and static controls. Only a sustained growth pattern
with resolution and horizon controls would motivate a parametric-instability
claim; do not infer it from one finite-time maximum.
