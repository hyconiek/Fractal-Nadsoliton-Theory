# P466 Current Strict Pair12 Chart‑Glued Projector Operator Section Audit Probe

Status: `P466_DRAFT_EXPECTED_EXECUTED_BY_p466_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Audit that the exported two‑chart projector operator section (`F462`) is numerically consistent:

1. `O_12` is orthogonal and involutive (`F461`),
2. the full-carrier projector gluing holds:
   `A_2 ≈ O_12 A_1 O_12^T`,
3. the pair-plane 2×2 gluing holds under `G(alpha_12)`:
   `A_2(c2,s2) ≈ G12 A_1(c1,s1) G12^T`,
4. residual sign flips `u -> -u` do not change `A_1` nor `A_2` (projector-level gauge).

This probe is **lane-scoped** and does not promote any global selector atlas/closure claim.

## Inputs

1. `A_1_pair1_orientation_projector_operator_strict_core_v1` (`F456`)
2. `O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1` (`F461`)
3. `A_2_pair2_orientation_projector_operator_strict_core_v1` (`F462`)

## Outputs

- `fundamental_action_reconstruction/generated/p466_current_strict_pair12_chart_glued_projector_operator_section_audit_probe.json`
- `fundamental_action_reconstruction/generated/p466_current_strict_pair12_chart_glued_projector_operator_section_audit_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. a global selector atlas / overlap-domain declaration / cocycle data,
2. global discharge of `QW-2191`,
3. a sign-sensitive physical orientation datum,
4. strict-core selector closure / admissible `S_sel_int`,
5. ToE closure.

