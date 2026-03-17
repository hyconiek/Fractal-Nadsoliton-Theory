# F690 Current Strict `T175`: Global Chart Sign Fixing From Strict-Core Payload Weights Export Packet (No False-PASS)

Status: `F690_CURRENT_STRICT_T175_GLOBAL_CHART_SIGN_FIXING_FROM_STRICT_CORE_PAYLOAD_WEIGHTS_EXPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Export an explicit **chart-level** sign-fixing (a `Z2` 0-cochain) rule and the resulting **sign-fixed directed representative** on `C_v1`,
constructed deterministically from already exported strict-core per-site payload weights.

This is a `strict_convention` gauge-fixing layer. It does not claim any strict physical sign datum and does not imply kernel-alone/global `QW-2191` discharge.

## Inputs reused

1. `F647` strict witness provider payload weights:
   - `fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
2. Premise-based directed representative on `C_v1`:
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_strict_v1.json`

## Exported artifacts

Executed by:

```bash
python3 fundamental_action_reconstruction/f690_current_strict_t175_global_chart_sign_fixing_from_strict_core_payload_weights_export_packet.py
```

Exports:

1. sign-fixed directed selector state object (strict_convention):
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json`
2. summary:
   - `fundamental_action_reconstruction/generated/f690_current_strict_t175_global_chart_sign_fixing_from_strict_core_payload_weights_export_packet_summary.json`

## Meaning (no false-PASS)

This packet means only:

1. a deterministic chart sign-fixing rule from exported strict-core payload weights is now exported as explicit data, and
2. the resulting sign-fixed directed representative on `{pair1..pair5}` is exported as a convention-layer state object.

## Hard limits

This packet does **not** claim:

1. strict-core physical sign/orientation datum,
2. `Aut(Z_12)`-invariant sign canonicity (`N462`),
3. kernel-alone/global `QW-2191` discharge,
4. operator-level groupoid/cocycle identities (`N512` boundary),
5. ToE closure.

