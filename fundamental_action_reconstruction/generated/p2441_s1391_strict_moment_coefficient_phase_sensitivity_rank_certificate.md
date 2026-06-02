# P2441 S1391: strict moment coefficient phase-sensitivity rank certificate

Status: `PASS_STRICT_MOMENT_COEFFICIENT_PHASE_SENSITIVE_NO_GENERATOR`

## Finite facts

- Coefficients: `['lambda_sm_eff', 'kappa_gr_eff', 'epsilon_mix_eff']`.
- Parameters: `['omega', 'phi', 'beta', 'eta']`.
- Jacobian rank: `3`.
- Phi column norm: `14.221800887991094`.
- Phi sweep rows: `4`.
- Every phi sweep row changes a coefficient: `True`.
- P2440 phase-null inherited: `True`.

## Hard limits

- P1562-style moment coefficients are locally phase-sensitive at the current strict tuple.
- A phase-null coefficient ansatz cannot replace the moment route without a separate phase-invariance theorem.
- This sensitivity certificate is not a Standard-Model/GR physical-value generator.
- No QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.

## Next honest step

Construct a phase/topology-sensitive strict coefficient map, or prove an explicit phase-invariance theorem identifying which coefficients may lawfully ignore phi.

## Gatekeepers

`{'rg_audit_ran': True, 'three_coefficients': True, 'four_parameters': True, 'jacobian_rank_three': True, 'phi_column_nonzero': True, 'phi_sweep_four_rows': True, 'phi_sweep_changes_coefficients': True, 'p1562_replay_close': True, 'p2440_phase_null_inherited': True, 'p2440_not_replacement': True, 'no_phase_invariance_theorem_export': True, 'no_value_generator_export': True, 'no_qw2191_discharge': True, 'no_ltotal_export': True, 'no_toe_export': True, 'fingerprint_stable': True}`
