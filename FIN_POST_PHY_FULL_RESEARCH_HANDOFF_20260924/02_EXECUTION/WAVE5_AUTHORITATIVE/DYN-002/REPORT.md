# DYN-002 — Derive coarse mobility and noise from a declared microscopic update

Status: **CONDITIONAL_KINETIC_COEFFICIENTS_DERIVED__FULL_MARKOV_CLOSURE_REMAINS_UNSOURCED**.

## Declared microscopic law

This task must not pretend that FIN statics selects kinetics. We therefore freeze one explicit reversible local update law only as a conditional model. On a collision-resolved phase graph with effective energy `E(x)`, let coordinate `x_i` propose `x_i -> x_i +/- eps` with rates

`q_i^+ = M0/(beta eps^2) exp[-beta (E(x+eps e_i)-E(x))/2]`,
`q_i^- = M0/(beta eps^2) exp[-beta (E(x-eps e_i)-E(x))/2]`.

The rate ratio satisfies detailed balance exactly:

`q_i^+(x) / q_i^-(x+eps e_i) = exp[-beta (E(x+eps e_i)-E(x))]`.

This law is **not** derived from FIN; it is the declared kinetic primitive demanded by the task.

## Small-step coarse generator

Taylor expansion gives, uniformly on a smooth compact region,

`L_eps f = sum_i M0[-(partial_i E)(partial_i f) + beta^-1 partial_i^2 f] + O(eps)`.

Thus the limiting conditional SDE is

`dx_i = -M0 partial_i E dt + sqrt(2 M0/beta) dW_i`.

The mobility and noise are therefore not free once this microscopic rule is frozen:

`mobility = M0`, `diffusion = M0/beta`.

Changing the microscopic attempt scale `M0 -> c M0` multiplies **both** drift and diffusion by `c`, leaves the Gibbs stationary law `exp(-beta E)` unchanged, and rescales time. This reproduces the repository's clock-gauge no-go in an explicit coarse derivation rather than by assertion.

## Projection to the locked common phase

Inside the PHA-001 amplitude box, the transverse phase Hessian has certified gap

`g_perp >= 1.028201984333153e-5`.

Under the declared overdamped mobility, the linearized transverse relaxation rate is at least `M0 g_perp`, so a necessary averaging condition is

`t >> tau_perp = 1/(M0 g_perp)`

while the common phase and graph geometry change little over that interval. In the units of this chosen update law, `tau_perp <= 97257.1/M0` from the certified lower bound. This is a dimensionless kinetic consequence, not seconds.

The common phase `alpha` then inherits a gradient/noise generator obtained by projecting the full reversible generator onto the locked tangent. Its numerical mobility depends on the microscopic mobility metric and therefore on the chosen update law. It is not a new universal FIN constant.

## Why the full phase/amplitude Markov closure is not yet certified

Two independent obstructions remain:

1. the programme does not provide a quantitative radial slow/fast separation covering the entire intended amplitude dynamics; and
2. REF-002 proves that eliminating inertial interior degrees of freedom can create frequency-dependent response memory. If such hidden modes are present, their memory kernel must remain explicit rather than being folded into a fitted Markov mobility.

Therefore this task derives a valid **conditional** coarse generator for the declared reversible jump law and a phase-locking timescale condition, but it does not certify a universal memoryless phase+amplitude dynamics.

## Held-out attempt-rate law

Replacing `M0` by `3 M0` leaves the invariant distribution and all static FIN quantities unchanged while multiplying relaxation, mobility and diffusion by three. This is the requested clean separation between static invariants and kinetic factors.

## Allowed conclusion

Given an explicit reversible microscopic rule, coarse mobility/noise and a conditional phase projection can be derived. FIN still does not select that rule, an absolute clock, or a universal Markov closure.
