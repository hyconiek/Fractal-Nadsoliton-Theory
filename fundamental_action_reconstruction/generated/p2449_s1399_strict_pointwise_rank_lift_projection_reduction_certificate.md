# P2449/S1399 strict pointwise rank-lift projection-reduction certificate

Status: `PASS_STRICT_POINTWISE_PROJECTION_REDUCTION_MATCHES_CENSUS_NO_SELECTOR_SOURCE_THEOREM`

## Linear-algebra reduction

Projection vector norm: `0.050792746818405954`.
Maximum determinant-vs-projection identity error on audited samples: `3.4694469519536142e-17`.

## Stationary root reconstruction

Zero-projection roots: `[0.21223939395855995, 4.056893978208217]`.
Stationary-factor roots: `[0.7852889044322362]`.
All reconstructed roots match P2448: `True`.
Best reconstructed root: `{'root_d': 0.7852889044322362, 'residual': 2.553296116203363e-17, 'iterations': 34, 'final_width': 5.81756864903582e-14, 'root_family': 'stationary_factor', 'normalized_rank_lift_volume': 0.022155470302414403}`.

## Guardrail

This exports the determinant-to-projection reduction and an analytic stationary-factor replay only.  It exports no exact interval root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.
