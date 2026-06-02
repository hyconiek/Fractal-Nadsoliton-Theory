# P2444 S1394: strict moment rank-lift conditioning certificate

Status: `PASS_STRICT_MOMENT_RANK_LIFT_CONDITIONING_FRONTIER_NO_SOURCE_THEOREM`

## Finite facts

- Inherited candidates: `12`.
- Robust rank-lift candidates: `6`.
- Best-conditioned candidate: `K_at_d_1`.
- Best normalized volume: `0.020925092943529396`.
- Weakest-conditioned candidate: `K_at_d_20`.
- Robust candidates by kind: `{'kernel_sample': ['K_at_d_1', 'K_at_d_0', 'K_at_d_0.5', 'K_at_d_2', 'K_at_d_0.25'], 'raw_moment': ['M0']}`.

## Hard limits

- Numerical rank-lift conditioning is not an admissibility theorem for a supplemental constraint.
- The best-conditioned singleton is still only a candidate until a strict observable/source theorem licenses it.
- No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.

## Next honest step

Try to prove strict source admissibility for the best-conditioned rank-lift singleton, or explain why it is a gauge choice rather than a physical observable.

## Gatekeepers

`{'rg_audit_ran': True, 'twelve_candidates_inherited': True, 'twelve_conditioned_rows': True, 'best_candidate_k_at_d_1': True, 'robust_candidate_count_six': True, 'weakest_candidate_k_at_d_20': True, 'best_volume_above_threshold': True, 'weakest_volume_below_threshold': True, 'no_observable_source_export': True, 'no_value_generator_export': True, 'no_qw2191_discharge': True, 'no_ltotal_export': True, 'no_toe_export': True, 'fingerprint_stable': True}`
