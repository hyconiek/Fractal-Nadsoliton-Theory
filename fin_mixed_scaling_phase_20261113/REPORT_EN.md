# Mixed-channel amplitude and phase campaign

## Outcome

At matched slow time tau=epsilon^2 t in [0,1], reducing epsilon from
0.1 to 0.075 lowers the mixed-coefficient discrepancy for both tested phases
and both action estimates in the primary window. The discrepancy remains
positive. This supports a finite-amplitude explanation but neither establishes
its asymptotic order nor proves that its limit is zero.

The frozen theoretical coefficient is kx=0.258892167626. No coefficient was
retuned. Primary fits use tau in [0.25,1] and common regression nodes spaced
approximately one physical time unit apart:

| epsilon | psi phase (rad) | T | kx from I | kx from J | I deviation | J deviation |
|---|---|---|---|---|---|---|
| 0.075 | 0.000000 | 177.777778 | 0.2610087569 | 0.2612718734 | 0.8176% | 0.9192% |
| 0.075 | 1.570796 | 177.777778 | 0.2610574778 | 0.2614233594 | 0.8364% | 0.9777% |
| 0.1 | 0.000000 | 100.000000 | 0.2629731243 | 0.2627854226 | 1.5763% | 1.5038% |
| 0.1 | 1.570796 | 100.000000 | 0.2627943184 | 0.2633094301 | 1.5072% | 1.7062% |

The second window tau in [0.5,1], all integration-grid controls, and all
regression-sampling controls are retained in comparison.json. They are
sensitivity diagnostics, not confidence intervals or independent replicates.

## What was held fixed and what changed

This is the same full nonlinear two-field RP model as the preceding package,
with w=1.5 and g4=g6=1. Let f=sqrt(3/4) sech^2(x) and
g=sqrt(3/2) sech(x)tanh(x). Initial data are eta=epsilon g,
eta_t=0, psi=epsilon cos(theta) f, psi_t=-sqrt(w) epsilon sin(theta) f,
with theta=0 or pi/2. The quadratic physical actions are phase matched;
canonical initial actions need not coincide because the nonlinear coordinate
map depends on the initial phase. Actual canonical actions are used, never
replaced with the quadratic initial values.

Fourier collocation uses dx=0.1, periodic [-L,L), L=ceil(T)+30,
T=epsilon^(-2), and RK4 with steps=round(T/0.02), dt=T/steps.
Thus both the observation horizon and domain grow at smaller amplitude.
No explicit dealiasing is applied. The old inverse canonical map is reused.
The physical outward flux is integrated at RK4 stages at radii 10 and 20;
it is not defined by subtracting local energies. Approximately 0.1-spaced
canonical diagnostics are retained, including the exact final time.

## Estimator separation

The candidate law is I'=-ka I^2-kx IJ, J'=-kb J^2-kx IJ. For each
window we integrate IJ and form I(t0)-I-ka integral(I^2), and its J analogue.
Each is regressed on integral(IJ) with a free intercept, keeping ka and kb
at their theoretical values. These are conditional effective estimates;
self-channel corrections can bias the inferred mixed coefficient.

To isolate quadrature sensitivity, integration uses every 1st, 2nd, 5th or
10th dense sample while regression nodes stay fixed. The largest absolute
kx change between approximately 0.1 and 1.0 integration spacing is
6.7821608e-06, or 0.002620% of theory. It is too small in these controls
to explain the observed primary discrepancy. Separately, with dense
integration fixed, changing regression-node spacing from approximately
0.1 to 1.0 changes kx by up to 0.00055082667
(0.212763% of theory). Thus decimating an entire trace mixes
two effects; it is not a pure quadrature test. The finest diagnostic interval
itself has not been halved in this campaign.

## Numerical integrity and limits

Maximum relative global-energy drift is 7.359336e-07; the maximum recorded
local-energy-plus-outward-flux residual relative to initial energy is
2.8766586e-07. These are consistency diagnostics, not rigorous error bounds.
The preceding campaign supplied spatial/time refinement at epsilon=0.1;
this new epsilon=0.075 horizon has not received its own mesh refinement.
Two amplitudes and two phases cannot determine an asymptotic exponent,
exclude a nonzero limit, or cover all phases. No stochastic error bars are
appropriate for these deterministic runs. Shared RK4 and canonical code
mean this is not an independent-team replication.

## Scope and next decision

Status: FINITE_AMPLITUDE_TREND_OBSERVED; MIXED_CHANNEL_DISCREPANCY_OPEN.
AW remains a conditional theorem for its historical preparation, pending
independent review. These auxiliary shape/mixed initial data do not extend
its theorem. This campaign is numerical work on a supplied classical PDE,
not a derivation from the strict FIN core or an empirical test of FIN.
No photons, SI units, natural action, SM/GR or ToE closure is inferred.

The next bounded numerical step is a third amplitude epsilon=0.05 at the
same slow-time window, with a smaller-amplitude resolution/time-step control
before any extrapolation. Use the separated quadrature and sampling tests;
do not launch another chain of same-author proof audits. The FIN bridge
remains a separate priority: a strict-core nonlinear response tensor and
state/clock/locality map must constrain the PDE couplings. Linear spectral
matching alone cannot select g4 and g6, as the existing bridge counterfamily
already demonstrates; the observed damping trend does not discharge that gap.
