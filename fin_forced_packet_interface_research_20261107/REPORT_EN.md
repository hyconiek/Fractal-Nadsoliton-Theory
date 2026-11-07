# Forced histories and the actual cubic spectral split

The causal packet estimate now includes packets born at different source
times. Its sufficient condition is an integrated source-time atomic mass
M_T=O(epsilon). Under that condition the full triangular cubic comparison
has integrated squared trace O(epsilon^2), including mixed terms and the
initial transients of each quadratic response.

Applying this condition to the whole physical source requires a spectral
residual. The localized f^3 source profile has positive spectral mass in
1/2<=|k|<=1, whereas the broad near-threshold atoms have exponentially
small mass there. An explicit lower bound quantifies the required growth
of their coefficient mass. The f^3 response already has the separate BT
localized estimate, so this calculation guides the decomposition of the
source before applying the packet interface.

A further estimate applies directly to CS's actual cubic forcing. For any
fixed 0<kappa<=1, the output with |k|>=kappa has integrated squared trace
at most C_S epsilon^2/kappa. Combining this with CS's local-continuum
inequality locates the remaining critical trace estimate in |k|<kappa.
Constants deteriorate when kappa shrinks with epsilon. The intermediate
bands and the threshold core are the next research obligation.

## Numerical checks

The unitary Darboux transform of the projected f^3 profile agrees with
direct quadrature of its differential expression. Its norm on the stated
annulus is approximately 0.0873428493040 at both tested quadrature orders.
For widths in [R,1.5R], mass squared 16/3, and an allowed annular error
of half that norm, the necessary coefficient masses obey these bounds.

| R | Necessary coefficient mass, lower bound |
|---|---|
| 4 | 0.06969 |
| 8 | 19.88 |
| 16 | 3.724e11 |
| 32 | 1.296e53 |

These bounds describe the normalized spatial profile. They do not estimate
the magnitude of the final actual nonlinear field.

A momentum-forced spectral oscillator with a closed-form resonant solution
checks the causal source-time integral. Its relative errors at 256, 512
and 1024 time subintervals are 3.241e-6, 8.100e-7 and 2.025e-7.
The reduction by approximately four matches second-order quadrature.
All 31 recorded checks passed, including predecessor provenance.

The analytical claims use the supplied RP model and the predecessor
lemmas. Independent review and formal proof verification are outstanding.
The actual near-threshold cubic estimate and modal-action errors are open.
