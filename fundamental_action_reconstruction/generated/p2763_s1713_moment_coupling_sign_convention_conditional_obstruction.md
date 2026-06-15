# P2763/S1713 moment-to-coupling sign-convention conditional obstruction

Status: `P2763_MOMENT_COUPLING_SIGN_CONVENTION_CONDITIONAL_OBSTRUCTION_NO_CLOSURE`

## Sign rows
- lambda_sm_eff: source_sign=1; coupling_sign=1; rectification=False; theorem_exported=False
- kappa_gr_eff: source_sign=-1; coupling_sign=1; rectification=True; theorem_exported=False
- epsilon_mix_eff: source_sign=-1; coupling_sign=1; rectification=True; theorem_exported=False

## Finite branch family
- branch_family_count=4
- all_branches_preserve_magnitudes=True
- unique_branch_selected_by_current_artifacts=False

## Decision
A finite four-branch sign family preserves the P1562 coupling magnitudes, while current artifacts do not export a rule selecting the positive kappa/epsilon branch. Because P2762 leaves reference-cell/action-density normalization open, even a future sign rule would have to be explicitly conditional rather than physical-coupling closure.

## Recommendation
Do not promote the positive P1562 kappa/epsilon signs to physical coupling signs.  The next honest proof-grade move should target the remaining provenance atom that is not merely sign choice: field/curvature normalization compatibility with the still-open canonical reference-cell/action-density theorem.  If that object cannot be supplied, preserve the P2697-P2763 no-closure certificate.
