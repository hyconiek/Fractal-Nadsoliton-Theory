# P2764/S1714 field-curvature normalization compatibility obstruction

Status: `P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE`

## Finite normalization orbit
- distinct_normalized_coefficient_triples=5
- all_coefficients_normalization_sensitive=True

## Rows
- lambda_sm_eff: term=lambda_sm_eff * phi^4; sensitive=True; accepted=False
- kappa_gr_eff: term=kappa_gr_eff * R; sensitive=True; accepted=False
- epsilon_mix_eff: term=epsilon_mix_eff * phi^2 * R; sensitive=True; accepted=False

## Decision
The finite normalization orbit shows that positive scalar/curvature normalizations change all three candidate coefficients. Current artifacts do not export a common field-curvature normalization theorem, and P2762/P2763 keep reference-cell and sign provenance open.

## Recommendation
Do not promote field or curvature rescaling conventions to physical coupling provenance.  The next honest proof-grade move is the remaining variational-insertion atom: a nonproxy EOM residual insertion test explicitly conditional on future reference-cell/action-density, sign, and field/curvature normalization theorems.  If that cannot be supplied, preserve the P2697-P2764 no-closure certificate.
