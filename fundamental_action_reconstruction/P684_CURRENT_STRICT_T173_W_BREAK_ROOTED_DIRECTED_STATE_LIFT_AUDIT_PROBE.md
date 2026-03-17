# P684 Current Strict `T173` `w_break`-Rooted Directed State Lift Audit Probe (No False‑PASS)

Status: `P684_CURRENT_STRICT_T173_W_BREAK_ROOTED_DIRECTED_STATE_LIFT_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit whether the repo can construct a **global directed (sign-sensitive) vector representative section** on `C_v1`
by fixing the residual `Z2` sign on `pair1` using the exported strict-core reflection-breaking weight payload
`w_break_by_x` (from the constructed internal source-object witness provider `F647`) and propagating that choice
to other charts via the exported rooted chart-transport operators `O_1m`.

This probe is **convention-scoped**:

- it does not claim `Aut(Z_12)`-invariant sign canonicity,
- it does not claim any strict physical orientation datum,
- it does not claim `QW-2191` discharge,
- it does not claim ToE closure.

## Inputs

- `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json` (`w_break_by_x`)
- `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` (`u_m`, `coords`, `A_m`)
- `generated/o1m_pair1_pairm_selector_chart_transport_operator_*.json` (rooted transports `O_12,O_13,O_14,O_15`)

## Output

- `generated/p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe.json`
- `generated/p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe_summary.json`
