# F688 Current Strict `T174`: Global Oriented Transition (Edge Sign‑Lift) Export Packet (No False‑PASS)

Status: `F688_CURRENT_STRICT_T174_GLOBAL_ORIENTED_TRANSITION_EDGE_SIGN_LIFT_EXPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Export an explicit **edgewise oriented lift** (`α mod 2π`) of the already exported strict global transition object on `C_v1` (`F469`) by
adding the missing sign choice `s_ij ∈ {±1}` on each overlap edge.

This export is **strict_convention** only: it is a tracked gauge/section choice, not a strict physical sign datum and not a global
`QW-2191` discharge.

## Inputs reused

1. `F469/N515`
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_strict_v1.json`
2. `F684`
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json`

## Exported artifacts

Executed by:

```bash
python3 fundamental_action_reconstruction/f688_current_strict_t174_global_oriented_transition_edge_sign_lift_export_packet.py
```

Exports:

1. global oriented transition edge sign‑lift object (strict_convention):
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_w_break_rooted_directed_state_strict_convention_v1.json`
2. summary:
   - `fundamental_action_reconstruction/generated/f688_current_strict_t174_global_oriented_transition_edge_sign_lift_export_packet_summary.json`

## Meaning (no false‑PASS)

This packet means only:

1. the repo now exports an explicit **edgewise** sign‑lift `s_ij` that upgrades the axis-only transition representatives to an oriented lift,
2. the lift is constructed to transport the already exported directed representative section without sign flips,
3. this is a convention layer and does **not** upgrade any directed sign into strict-core physics.

## Hard limits

This packet does **not** claim:

1. a strict-core physical sign/orientation datum,
2. `Aut(Z_12)`-invariant canonicity,
3. operator-level groupoid/cocycle identities (`N512` boundary),
4. kernel-alone/global `QW-2191` discharge,
5. ToE closure.

