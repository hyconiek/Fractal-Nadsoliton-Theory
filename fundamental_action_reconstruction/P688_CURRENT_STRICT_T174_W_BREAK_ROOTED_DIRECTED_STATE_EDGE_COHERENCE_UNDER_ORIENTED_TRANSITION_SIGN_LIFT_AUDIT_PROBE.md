# P688 Current Strict `T174`: W‑Break Rooted Directed State — Full‑Edge Coherence Under Oriented Transition Sign‑Lift (Audit Probe)

Status: `P688_CURRENT_STRICT_T174_W_BREAK_ROOTED_DIRECTED_STATE_EDGE_COHERENCE_UNDER_ORIENTED_TRANSITION_SIGN_LIFT_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit that the exported oriented edge sign‑lift (`F688`) removes the edgewise sign flips found under the base axis-only global transition object (`P686/N686`)
for the exported `w_break`‑rooted directed representative section (`F684`).

Concretely, this probe checks:

```text
for every edge pairi_to_pairj:
  (s_ij * O_ij) u_i ≈ u_j
```

where `O_ij` is the exported axis-only transition operator from `F469` and `s_ij` is the exported sign-lift from `F688`.

## Inputs

1. `F688` oriented sign-lift object:
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_w_break_rooted_directed_state_strict_convention_v1.json`
2. `F684` directed state representative:
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json`

## Output

Writes:

- `fundamental_action_reconstruction/generated/p688_current_strict_t174_w_break_rooted_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe.json`
- `fundamental_action_reconstruction/generated/p688_current_strict_t174_w_break_rooted_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- strict-core physical sign/orientation datum,
- kernel-alone/global `QW-2191` discharge,
- operator-level groupoid identities (`N512` boundary),
- ToE closure.

