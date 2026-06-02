# P2442 S1392: strict moment coefficient local-identifiability nullspace certificate

Status: `PASS_STRICT_MOMENT_COEFFICIENT_LOCAL_NULLSPACE_NO_GENERATOR`

## Finite facts

- Input Jacobian rank: `3`.
- Nullspace dimension: `1`.
- Max linear null residual: `6.38378239159465e-16`.
- Max coefficient delta under null perturbations: `4.161198759788931e-06`.
- Max kernel sample delta under null perturbations: `0.00013878834148423058`.
- Extra constraint/gauge theorem required: `True`.

## Hard limits

- Three moment coefficients cannot locally identify four strict kernel parameters without an extra premise.
- The null direction is not proven to be a gauge redundancy by this certificate.
- Null-direction first-order coefficient stability is not SM/GR physical-value generation.
- No QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.

## Next honest step

Add one independent strict observable/source constraint for the null direction, or prove a gauge-redundancy theorem that quotients it before promoting moment coefficients.

## Gatekeepers

`{'rg_audit_ran': True, 'rank_three_inherited': True, 'one_dimensional_nullspace': True, 'null_residual_small': True, 'null_perturbation_rows_two': True, 'null_perturbation_coefficients_stable': True, 'null_perturbation_changes_kernel_samples': True, 'not_locally_injective': True, 'extra_constraint_required': True, 'no_coefficient_theorem_export': True, 'no_value_generator_export': True, 'no_qw2191_discharge': True, 'no_ltotal_export': True, 'no_toe_export': True, 'fingerprint_stable': True}`
