# Critical cubic trace: improved bound and a sharp limitation

Date: 2026-11-04. Conditional analytical component result;
actual critical-time damping remains open.

The actual cubic-forcing contribution to VD's compact threshold trace
satisfies `integral_0^(S epsilon^-2) |q_ge3|^2 <= C_S epsilon`.
This improves the O(1) consequence of the previous pointwise estimate.
A finite-time free Klein-Gordon observation estimate is applied to the
actual QP forcing, without assuming a slowly varying amplitude.

The stronger O(epsilon^2) target cannot follow from these norm bounds
alone. An explicit even, zero-initial-data comparison family has phase
energy O(epsilon), forcing O(epsilon^3), and integrated compact trace
between epsilon/80 and 4epsilon/5 through time epsilon^-2. The new order
is therefore sharp for that forcing class. The comparison forcing is NOT
the actual nonlinear FIN forcing and is never inserted into the model.
This identifies a missing source-specific cancellation, not a failure
of the prepared solution or a universal no-go for damping.

## Continuity correction

The initial CR lookup selected ST but missed CD and VD, already present
in the 2026-10-21 source commit. CD derives all three frozen radiative
channels; VD supplies an actual continuum virial with one compact-mean
defect. Those documents are now read and bound. No channel calculation
is claimed as new. VD's NEXT_TRACE_AND_MODAL_GATE.md is the predecessor.
CR's single-channel candidate remains a provisional diagnostic.

## Remaining work

Derive near-threshold cancellation for the actual transformed cubic
forcing. Separately control the localized-history trace q_loc and the
full three-channel modal action discrepancy. Even an O(epsilon) whole
continuum budget would not make the leading absolute action-error bound
o(epsilon^2). Positive reference coefficients are not actual damping.

## Validation and scope

The proof is in PROOF_OR_OPEN_LEMMA.md. The standard-library checker
tests the differentiated comparison equation, initial conditions,
scaling, oscillatory bounds, mutation controls and bounded floating
previews. Source hashes are checked separately. These are not PDE
simulations, outward interval certificates or independent review.

All new content is in English and publication is explicitly authorized.
Historical AGENTS and scientific source files are unchanged. No action,
parameter, initial wave, bath, experiment, physical clock or source
selection is added. F01 and all physical interpretation gates remain open.
