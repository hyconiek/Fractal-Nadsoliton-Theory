# FIN ST2802–ST2891 — Five-Axis Operational Completion

## Main result

The smallest declared operational model containing all currently separated
obstructions has five parameters:

\[
z=(\theta,q,a,c,j),
\]

representing equilibrium ternary coupling, prism refinement, kinetic activity,
clock scale and signed circulation.

A corresponding observable vector was constructed from:

1. a stationary triangle moment,
2. the prism-renormalized coupling,
3. an untimed jump probability,
4. a calibrated holding time,
5. a signed cycle current.

At the declared interior point \((0.3,0.2,0.1,2,0.05)\), the Jacobian has

\[
\boxed{\operatorname{rank}J=5},
\qquad
\det J=-3.3901090721157604\times10^{-4}.
\]

Thus the five source axes are locally identifiable and independent in this
model. No smooth one-parameter law can cover an open neighbourhood of the full
five-dimensional operational family.

## Record hierarchy

| Available records | Rank |
|---|---:|
| stationary state and prism equilibrium | 2 |
| plus untimed jump sequence | 3 |
| plus calibrated waiting times, no arrow | 4 |
| plus signed current, no clock | 4 |
| complete record | 5 |

Different observations remove different gauge fibres. A timestamp cannot
replace an arrow instrument, and a signed current cannot calibrate a clock.

## Interpretation

The result refutes the idea that the present mathematics is missing only one
number. Within the declared smooth operational class, at least five typed
objects remain:

- state coupling \(\theta\),
- vertical/refinement coupling \(q\),
- kinetic activity shape \(a\),
- clock scale \(c\),
- arrow polarity/current \(j\).

This is a scoped local-rank theorem, not an absolute lower bound over all
future mathematical formulations. A deeper theorem could reduce the rank by
coupling several axes, but that coupling must be derived rather than asserted.

## Operational package

The minimal record must distinguish preparation, registration, calibrated
clock, signed cycle instrument and frozen analysis. Suggested event fields are:

`run`, `layer`, `configuration_before`, `configuration_after`, `waiting_time`,
`timestamp`, `cycle_orientation`, and `clock_calibration`.

Hashing a file does not create independent custody, and a synthetic record is
not physical evidence.

## Gate

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{strict physical rows}.}
\]

The highest-value next theorem is one that genuinely couples two or more of
the five axes and thereby lowers the rank. Without such a theorem, the missing
bridge must be treated as a typed operational package rather than one hidden
scalar or one spectral identity.

No selector closure, physical clock, units, apparatus evidence, Standard
Model, gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
