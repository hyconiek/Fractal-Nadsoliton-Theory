# Background phase and transverse-gap comparison

This campaign extends the finite-time j=6 variational test, not the static
stability theorem. Nine runs comprise one static control at w=1.25 and four
driven preparations, each on two joint space/time resolutions.

## Fixed experiment

The background is the full nonlinear RP model with epsilon=0.1, T=100,
L=130, eta(0)=epsilon sqrt(3/2)sech(x)tanh(x), eta_t(0)=0,
psi(0)=epsilon cos(theta)f and psi_t(0)=-sqrt(w)epsilon sin(theta)f,
where f=sqrt(3/4)sech^2(x). Values w=1.25,1.5 and theta=0,pi/2
are compared. Phase refers to the background psi field, not the transverse
seed phase. Physical quadratic background actions match between phases;
nonlinear canonical initial actions need not match exactly.

The transverse equation, static reference energy, two-column Gram matrix
and seed-space phase optimization are unchanged from
[the preceding experiment](../fin_transverse_radiating_20261113/REPORT_EN.md).
The gap is w-Delta_6, Delta_6=1.1519772119538363. Each initial seed pair
is normalized in its own positive static energy norm. Therefore comparisons
across w concern relative gain, not equal raw displacement amplitudes.
The variational columns do not back-react on the nonlinear background.

Coarse dx=0.1,dt=0.02 and fine dx=0.05,dt=0.01 test their combined effect.
The previous campaign separately checked time refinement at w=1.5, phase=0;
that separation is not independently repeated for every new preparation.
The Fourier/RK4 solver has no explicit dealiasing. Sampling interval is 0.5;
reported maxima are over stored samples. The static control verifies the
energy norm at the smaller gap. Earlier static control at w=1.5 is retained
in its original package. All raw traces and run metadata are included.

## Interpretation boundary

RESULTS.md and comparison.json contain refined gains and full-trace
refinement differences. Each maximum is optimized only over the two
localized initial displacement/velocity seeds, not arbitrary spatial data.
The optimal seed can differ between times. No logarithmic slope, Floquet
exponent or all-time instability claim is inferred from these finite runs.
Static positivity does not make a driven energy norm invariant, and a
smaller static gap need not produce a larger gain over a fixed time interval.

This remains a conditional A11/RP construction. Coefficients, wave kinetics,
preparation and clock are supplied; no strict FIN derivation, AW extension,
selector resolution, units, experimental confirmation or ToE closure follows.

## Replay

From the repository root with Python and NumPy:

    OPENBLAS_NUM_THREADS=1 python fin_transverse_phase_gap_20261113/run.py
    python fin_transverse_phase_gap_20261113/analyze.py

The checker verifies trace shape, finite positive Gram eigenvalues, initial
normalization, stored maxima, static conservation and joint refinement.
It does not prove a continuum error bound or nonlinear stability.
