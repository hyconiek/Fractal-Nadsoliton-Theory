# P691 Current Strict `T174`: Sign-Fixed Directed State Full-Edge Coherence Under Oriented Transition Sign-Lift (Audit Probe)

Status: `P691_CURRENT_STRICT_T174_SIGN_FIXED_DIRECTED_STATE_EDGE_COHERENCE_UNDER_ORIENTED_TRANSITION_SIGN_LIFT_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit that the oriented edge sign-lift exported by `F691` removes edgewise sign flips for the exported **sign-fixed** directed representative (`F690`)
on every overlap edge of the strict global transition object on `C_v1`.

Concretely, this probe checks:

```text
for every edge pairi_to_pairj:
  (s_ij * O_ij) u_i_fixed ≈ u_j_fixed
```

## Inputs

1. `F691` oriented sign-lift object:
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_sign_fixed_directed_state_strict_convention_v1.json`
2. `F690` sign-fixed directed state representative:
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json`

## Output

Writes:

- `fundamental_action_reconstruction/generated/p691_current_strict_t174_sign_fixed_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe.json`
- `fundamental_action_reconstruction/generated/p691_current_strict_t174_sign_fixed_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- strict-core physical sign/orientation datum,
- kernel-alone/global `QW-2191` discharge,
- operator-level groupoid identities (`N512` boundary),
- ToE closure.

