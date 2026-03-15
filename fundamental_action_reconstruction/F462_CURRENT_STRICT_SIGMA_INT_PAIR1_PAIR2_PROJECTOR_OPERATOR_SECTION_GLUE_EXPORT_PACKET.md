# F462 Current Strict Sigma-int Pair1/Pair2 Projector Operator Section Glue Export Packet

Status: `F462_DRAFT_EXPECTED_EXECUTED_BY_f462_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

After:

- `F456` exporting the minimal strict-core projector operator on `pair1`:
  `A_1(pair1) := |u_1><u_1|`,
- `F461` exporting the lane-scoped chart-transport operator `O_12` between `pair1↔pair2`,

the next narrow, honest B3→H39 continuation move is to **glue the projector-level operator data across the two charts**
in a sign-gauge-safe way, without claiming a global selector atlas or a sign-sensitive physical orientation datum.

Concretely, this packet exports:

1. the corresponding strict-core projector operator on `pair2`:
   `A_2(pair2) := |u_2><u_2|`,
2. one explicit two-chart **operator section** with the gluing law:

```text
A_2(pair2)  =  O_12 A_1(pair1) O_12^T
```

in the declared sigma-int corridor scope on the `n=12` real Fourier carrier.

## Strict inputs

1. `ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1` (`F451/N489`)
   - provides the lane-scoped unit representatives `u_1 ∈ span{c1,s1}` and `u_2 ∈ span{c2,s2}`.
2. `A_1_pair1_orientation_projector_operator_strict_core_v1` (`F456`)
   - provides `A_1(pair1) := |u_1><u_1|` as a strict-core exported operator object.
3. `O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1` (`F461`)
   - provides the explicit lane-scoped orthogonal chart-transport operator `O_12`.

## Construction

Define the `pair2` operator:

```text
A_2(pair2) := |u_2><u_2|.
```

Define the two-chart glued section on `{pair1, pair2}` by the compatibility law:

```text
A_2(pair2)  =  O_12 A_1(pair1) O_12^T.
```

All objects are **projector-level** and therefore residual-sign-gauge invariant (`u -> -u`).

## Exported artifacts

- `fundamental_action_reconstruction/generated/a_2_pair2_orientation_projector_operator_strict_core_v1.json`
- `fundamental_action_reconstruction/generated/a_12_pair12_chart_glued_orientation_projector_operator_section_strict_core_v1.json`
- `fundamental_action_reconstruction/generated/f462_current_strict_sigma_int_pair1_pair2_projector_operator_section_glue_export_packet_summary.json`

## Hard limits (no false pass)

This packet does **not** claim:

1. a global selector atlas / overlap-domain declaration / cocycle data (`H41` remains open globally),
2. a global selector transition/gluing object discharging `QW-2191` (`H40` remains open globally),
3. a sign-sensitive physical orientation datum (this is projector-level / sign-free),
4. strict-core selector closure / admissible `S_sel_int`,
5. ToE closure.

