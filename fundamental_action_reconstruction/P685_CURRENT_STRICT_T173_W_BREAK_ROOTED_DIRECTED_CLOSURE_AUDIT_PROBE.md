# P685 Current Strict `T173` `w_break`-Rooted Directed Closure Audit Probe (No False‑PASS)

Status: `P685_CURRENT_STRICT_T173_W_BREAK_ROOTED_DIRECTED_CLOSURE_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit that the exported convention-scoped global directed representative (`F684`) together with the promoted global output channels (`F660`)
induces a **chartwise directed closure output** that glues in the fixed output basis `(o_+,o_-)` on `{pair1..pair5}`.

This probe remains strictly below any claim of:

- a strict physical orientation datum,
- `Aut(Z_12)`-invariant sign canonicity,
- kernel-alone/global `QW-2191` discharge,
- ToE closure.

## Inputs

- `generated/selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json` (`F684`)
- `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json` (`F660`)
- `generated/selector_output_sign_lift_global_c_v1_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json` (`F683`) (consistency cross-check)

## Output

- `generated/p685_current_strict_t173_w_break_rooted_directed_closure_audit_probe.json`
- `generated/p685_current_strict_t173_w_break_rooted_directed_closure_audit_probe_summary.json`

