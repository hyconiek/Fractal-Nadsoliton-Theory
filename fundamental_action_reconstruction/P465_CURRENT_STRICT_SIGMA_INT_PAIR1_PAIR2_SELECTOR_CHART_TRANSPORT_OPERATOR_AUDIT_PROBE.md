# P465 Current Strict Sigma-int Pair1/Pair2 Selector Chart-Transport Operator Audit Probe

Status: `P465_DRAFT_EXPECTED_EXECUTED_BY_p465_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Audit that the lane-scoped chart-transport operator exported by `F461`:

1. is orthogonal and involutive on the declared `n=12` carrier,
2. transports the canonical `pair1` plane to the canonical `pair2` plane and vice versa,
3. transports a finite family of reduced/tangent directions `u_{m,θ}` / `τ_{m,θ}` as predicted by the `SO(2)` transition angle `alpha_12`,
4. yields sign-gauge-invariant projector transport on that sampled family.

This probe is strictly lane-scoped and is **not** a global selector atlas/gluing object (`H41`) and does **not** discharge `QW-2191` (`H40`).

## Inputs

1. `O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1` (`F461`).

## Outputs

- `fundamental_action_reconstruction/generated/p465_current_strict_sigma_int_pair1_pair2_selector_chart_transport_operator_audit_probe.json`
- `fundamental_action_reconstruction/generated/p465_current_strict_sigma_int_pair1_pair2_selector_chart_transport_operator_audit_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. global selector atlas/gluing data,
2. global discharge of `QW-2191`,
3. sign-sensitive physical orientation derived,
4. strict-core selector closure / admissible `S_sel_int`,
5. ToE closure.

