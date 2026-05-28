# P2302/S1252 — provider-lift policy-lock candidate

- Status: `OPEN_BRIDGE_OBLIGATION_WITH_CONDITIONAL_G1_G3_POLICY_LOCK_CANDIDATE_TRACE`
- Conditional candidate found: `True`
- Minimal provider_lift_per_step: `0.0068`
- Conditional gap rows: `[('G1_reduction_certainty', 'CONDITIONALLY_CLOSED_UNDER_PROVIDER_MARGIN_LIFT'), ('G2_nonlinear_trajectory_realism', 'CLOSED_FROM_P2301'), ('G3_operational_policy_rule', 'CONDITIONALLY_CLOSED_UNDER_PROVIDER_MARGIN_LIFT')]`
- Strict closure status: `HELD_OPEN_UNTIL_PROVIDER_TO_MARGIN_BRIDGE_PROVEN`
- Theorem fingerprint: `7c8da2c6d1100354192abf32e99e4a0f94941e5c28c7c68602c980544e762c34`

## Guardrail statement
P2302 finds a numerical conditional G1/G3 policy-lock candidate, but it does not prove the provider-to-margin bridge from P2300 coefficients. Therefore strict Task-3 closure remains open; no QW-2191, selector, legacy-kernel, or ToE closure is claimed.

## Next honest step
Prove or refute the provider-to-margin bridge: derive the P2302 provider_lift_per_step bound directly from the P2300 canonical spatial-EOM coefficients without adding a selector premise; only then update strict G1/G3.
