# MP7-030 — local mechanism, global crossing status, and kinetic nonuniqueness

Scientific state: **DONE_SYNTHESIS**.

## Static landscape

The certified simple fold lies near `g=3.51564471684` and has one soft mode.
MP7-025 proves `a<0`, `b>0`; MP7-026 pays the Lyapunov--Schmidt remainder and
shows for

`0 < epsilon=g-g_fold <= 10^-7`

exactly two local stationary branches with

- `|xi|/sqrt(epsilon) in [1.8619243,1.9186339]`,
- local minimum--saddle energy splitting divided by `epsilon^(3/2)` in
  `[0.5126868,0.5586507]`,
- soft Hessian magnitude divided by `sqrt(epsilon)` in
  `[0.2149052,0.2349155]`.

These are controlled finite-error laws, not fitted exponents.

The previously local equal-energy event is now globally classified by
MP7-016/017:

`g_eq in [3.7183448971203875,3.7183448991203876]`.

For every smaller positive g the uniform state is the unique global minimizer.
At `g_eq`, the uniform state and twelve D12-related localized minima coexist
globally; the index-one saddle remains more than `0.04655` above them.  The
localized branch crosses with strictly negative energy slope, so above the
event it immediately beats the uniform energy.

Thus the fold/spinodal and the global coexistence event are distinct:

- fold: birth of a local minimum+saddle pair;
- coexistence: later equality of localized and uniform global energies.

## Static response

MP7-022 derives the response formulas under two explicitly different source
conventions.  MP7-026 implies the stable-branch susceptibility along the soft
direction diverges like `epsilon^-1/2` under a source conjugate to the dual
coordinate, with the precise coefficient depending on the selected
normalization/source observable.

## Added deterministic dynamics

Under the supplied gradient law

`dot s = -L grad Phi`,

MP7-028 shows that constant SPD mobility preserves the stationary set and
Hessian-inertia stability count while changing rates and paths.

For unit mobility, MP7-029 gives the controlled slow rate

`0.214905205154 sqrt(epsilon) <= gamma <= 0.234915431887 sqrt(epsilon)`,

hence

`4.25685/sqrt(epsilon) <= tau <= 4.65322/sqrt(epsilon)`.

For `L=alpha I`, rates are multiplied exactly by alpha.  Therefore the
square-root exponent is a property of the static fold geometry, whereas the
clock prefactor is not fixed by the potential.

## Added stochastic dynamics / finite N

The explicitly supplied Langevin law with noise covariance `2L/N` is
reversible with stationary density proportional to `exp(-N Phi)`.  MP7-031--035
connect this to the finite-copy equilibrium extension, but do not identify the
stochastic trajectory with a microscopic FIN dynamics.

After MP7-017 globally closes the coexistence set, MP7-035 gives the global
large-N rounding window

`g = g_eq + c/N`,

with localized-family/uniform weight ratio tending to

`exp(A-c DeltaV')`.

The equal-mass point is shifted by

`g_N = g_eq + c_*/N + o(1/N)`,

`c_* in [1.34676582095488,1.34676616197633]`.

## Nonconclusions

The static potential does not derive a unique dynamics, physical time unit,
noise amplitude, temperature, copy count N, or escape-rate prefactor.  The
local index-one saddle is not automatically the dominant global transition
state for any supplied kinetics.  No D12 orbit member is physically selected.
