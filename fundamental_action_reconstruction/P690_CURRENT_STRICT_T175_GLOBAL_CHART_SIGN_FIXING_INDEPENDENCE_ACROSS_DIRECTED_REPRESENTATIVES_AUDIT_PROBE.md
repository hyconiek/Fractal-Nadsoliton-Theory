# P690 Current Strict `T175`: Global Chart Sign Fixing Independence Across Directed Representatives (Audit Probe)

Status: `P690_CURRENT_STRICT_T175_GLOBAL_CHART_SIGN_FIXING_INDEPENDENCE_ACROSS_DIRECTED_REPRESENTATIVES_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit that the `T175` chart-level sign-fixing rule (a `Z2` 0-cochain) exported by `F690` yields the **same** sign-fixed directed representative on `C_v1`
when applied to multiple already exported directed representatives.

Concretely, this probe:

1. loads the sign-fixed directed representative exported by `F690`,
2. applies the same sign-fixing rule to the exported `w_break`-rooted directed representative (`F684`),
3. checks that the resulting sign-fixed vectors match `F690` on all `{pair1..pair5}` (numeric tolerance),
4. checks that the enforced positivity condition holds: `<w_i, u_i_fixed> > 0` for every chart used by the rule.

## Inputs

1. `F690` sign-fixed directed state object:
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json`
2. `F684` directed representative (w_break-rooted, convention scope):
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json`
3. `F647` strict witness provider payload weights:
   - `fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`

## Output

Writes:

- `fundamental_action_reconstruction/generated/p690_current_strict_t175_global_chart_sign_fixing_independence_across_directed_representatives_audit_probe.json`
- `fundamental_action_reconstruction/generated/p690_current_strict_t175_global_chart_sign_fixing_independence_across_directed_representatives_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- strict-core physical sign/orientation datum,
- `Aut(Z_12)`-invariant sign canonicity (`N462`),
- kernel-alone/global `QW-2191` discharge,
- operator-level groupoid identities (`N512` boundary),
- ToE closure.

