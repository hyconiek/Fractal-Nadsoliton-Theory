# P689 Current Strict `T174`: Oriented Edge Sign‑Lift Gauge‑Equivalence Across Directed Representatives (Audit Probe)

Status: `P689_CURRENT_STRICT_T174_ORIENTED_EDGE_SIGN_LIFT_GAUGE_EQUIVALENCE_ACROSS_DIRECTED_REPRESENTATIVES_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After `F688/P688/N688` exported a **convention-layer** oriented edge sign‑lift chosen to make the exported `w_break`‑rooted directed representative
transport without sign flips on all overlap edges, audit whether the induced edge sign‑lift pattern depends on the chosen directed representative
only by a chart‑level `Z2` gauge (0‑cochain) relift:

```text
u_i -> t_i * u_i,    t_i ∈ {±1}
```

Concretely, this probe:

1. extracts the `F688` edge sign‑lift pattern `s_ij^(A)`,
2. computes the edge sign‑lift pattern `s_ij^(B)` induced by the exported premise-based directed representative
   `SelectorState_global_C_v1_directed_strict_v1`,
3. checks whether there exists `{t_i}` such that for every overlap edge:

```text
s_ij^(B) = t_i * s_ij^(A) * t_j
```

## Inputs

1. `F688` oriented sign‑lift object:
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_w_break_rooted_directed_state_strict_convention_v1.json`
2. `F469` base strict global transition object:
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_strict_v1.json`
3. Premise-based directed representative:
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_strict_v1.json`

## Output

Writes:

- `fundamental_action_reconstruction/generated/p689_current_strict_t174_oriented_edge_sign_lift_gauge_equivalence_across_directed_representatives_audit_probe.json`
- `fundamental_action_reconstruction/generated/p689_current_strict_t174_oriented_edge_sign_lift_gauge_equivalence_across_directed_representatives_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- strict-core physical sign/orientation datum,
- kernel-alone/global `QW-2191` discharge,
- `Aut(Z_12)`-invariant sign canonicity (`N462`),
- operator-level groupoid identities (`N512` boundary),
- ToE closure.

