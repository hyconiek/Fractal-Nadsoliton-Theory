# F691 Current Strict `T174`: Global Oriented Transition (Edge Sign-Lift) From Sign-Fixed Directed State Export Packet (No False-PASS)

Status: `F691_CURRENT_STRICT_T174_GLOBAL_ORIENTED_TRANSITION_EDGE_SIGN_LIFT_FROM_SIGN_FIXED_DIRECTED_STATE_EXPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Re-export a `T174`-class global oriented transition edge sign-lift on `C_v1`, but anchored to the **sign-fixed** directed representative exported by `T175` (`F690`).

This removes dependence on an arbitrary starting directed representative: the input directed vectors are first stabilized by the exported chart sign-fixing (0-cochain),
then the edgewise signs `s_ij in {+-1}` are exported so that the stabilized directed representative transports without sign flips on every overlap edge.

This export is `strict_convention` only: it is explicit gauge/section data, not a strict physical sign datum and not a kernel-alone/global `QW-2191` discharge.

## Inputs reused

1. `F469/N515`
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_strict_v1.json`
2. `F690`
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json`

## Exported artifacts

Executed by:

```bash
python3 fundamental_action_reconstruction/f691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_export_packet.py
```

Exports:

1. global oriented transition edge sign-lift object (strict_convention):
   - `fundamental_action_reconstruction/generated/selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_sign_fixed_directed_state_strict_convention_v1.json`
2. summary:
   - `fundamental_action_reconstruction/generated/f691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_export_packet_summary.json`

## Hard limits

This packet does **not** claim:

1. a strict-core physical sign/orientation datum,
2. `Aut(Z_12)`-invariant canonicity,
3. operator-level groupoid/cocycle identities (`N512` boundary),
4. kernel-alone/global `QW-2191` discharge,
5. ToE closure.

