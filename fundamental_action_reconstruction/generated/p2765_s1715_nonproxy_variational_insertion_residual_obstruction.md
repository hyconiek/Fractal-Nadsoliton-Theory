# P2765/S1715 nonproxy variational-insertion residual obstruction

Status: `P2765_NONPROXY_VARIATIONAL_INSERTION_RESIDUAL_OBSTRUCTION_NO_CLOSURE`

## Residual support
- p1563_formal_eom_keys=['gr_sector', 'psi_sector']
- p1866_eom_export_keys=['eom_A_proxy_1d', 'eom_phi_proxy_1d']
- missing_nonproxy_residual_row_count=4

## Rows
- lambda_sm_eff: accepted=False; required=['scalar_4d_euler_lagrange_residual']; proxy_hits=['eom_phi_proxy_1d']
- kappa_gr_eff: accepted=False; required=['metric_4d_einstein_residual']; proxy_hits=[]
- epsilon_mix_eff: accepted=False; required=['mixed_scalar_metric_4d_residual', 'metric_backreaction_residual']; proxy_hits=[]

## Decision
P1563/P1866 provide formal/proxy EOM support, but the required 4D nonproxy scalar, metric, and mixed residual rows are absent; P2762-P2764 provenance prerequisites also remain open.

## Recommendation
Do not promote P1563/P1866 formal or proxy EOM exports to a nonproxy variational-insertion theorem.  The honest next move is a post-P2761-to-P2765 provenance-state reconciliation: either supply a genuinely new theorem fixing at least one of reference-cell/action-density, sign, field/curvature normalization, or 4D nonproxy residual rows, or preserve the P2697-P2765 no-closure certificate.
