# P2318 S1268: selector-lane to Task-3 margin bridge obligation audit

Status: `OPEN_SELECTOR_LANE_TO_TASK3_MARGIN_BRIDGE_FIELDS_MISSING`

## FAR bridge grep result
- Bridge/interface evidence hits: `22`.

## Required lift
- P2302 required `provider_lift_per_step`: `0.0068`.

## Lane scale witnesses
- Diagonal/local max defect: `25.5735`; target scale: `0.0002659`.
- Shannon max defect: `35`; target scale: `0.000194286`.
- The same lane data can pass, fail, or vanish under different unexported scale/sign/response maps.

## Missing bridge fields
- Missing fields: `signed_scalar_lift_per_step, time_step_normalization, margin_response_functional, orientation_to_policy_sign_rule, p2281_replay_semantics, admissibility_theorem_no_selector_premise`.
- Task-3 G1/G3 update: `NONE`.

## Theorem fingerprint
`239d954eeb53c29beaa205a1d8e8d708a90b5cd5eb293473d1be6cd885ffabf5`
