# P2443 S1393: strict moment supplemental-constraint rank-lift certificate

Status: `PASS_STRICT_MOMENT_SUPPLEMENTAL_CONSTRAINT_RANK_LIFT_FRONTIER_NO_SOURCE_THEOREM`

## Finite facts

- Inherited base rank: `3`.
- Inherited nullspace dimension: `1`.
- Candidate count: `12`.
- Rank-lifting candidate count: `12`.
- Minimal rank-lift antichain size: `1`.
- Rank-lifting candidates by kind: `{'raw_moment': ['M0', 'M1', 'M2', 'M3'], 'kernel_sample': ['K_at_d_0', 'K_at_d_0.25', 'K_at_d_0.5', 'K_at_d_1', 'K_at_d_2', 'K_at_d_5', 'K_at_d_10', 'K_at_d_20']}`.

## Hard limits

- Rank-lifting a Jacobian row is not proof that the row is an admissible strict observable/source constraint.
- Raw moments and pointwise kernel samples remain candidate constraints until a source theorem licenses them.
- No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.

## Next honest step

Choose one rank-lifting singleton and prove it is an admissible strict observable/source constraint or a lawful gauge-fixing condition for the P2442 null direction.

## Gatekeepers

`{'rg_audit_ran': True, 'inherited_rank_three': True, 'inherited_one_null_direction': True, 'twelve_candidates': True, 'all_candidates_rank_lift': True, 'minimal_singleton_rank_lift': True, 'raw_moments_rank_lift': True, 'kernel_samples_rank_lift': True, 'no_observable_source_export': True, 'no_coefficient_theorem_export': True, 'no_value_generator_export': True, 'no_qw2191_discharge': True, 'no_ltotal_export': True, 'no_toe_export': True, 'fingerprint_stable': True}`
