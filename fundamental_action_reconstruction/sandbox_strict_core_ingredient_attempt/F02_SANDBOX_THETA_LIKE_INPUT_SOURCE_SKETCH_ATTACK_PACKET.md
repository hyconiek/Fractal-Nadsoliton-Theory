# F02 Sandbox Theta-Like Input Source Sketch Attack Packet

Status: `F02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Refine the rho-oriented scaffold into something that names the minimal
theta-like input structure demanded by the real residual-orientation target
slot, while still remaining below any actual phase export.

## Theta-like input source sketch

Define:

```text
theta_like_input_source_sketch_v0 :=
(
  required_inputs_theta_1_theta_2,
  orientation_slice_class_level_candidate_present,
  basis_pair_export_skeleton_packet_ready,
  strict_core_actual_theta_export_absent,
  axiom_augmented_theta_source_branch_present_non_strict
)
```

where:

1. `required_inputs_theta_1_theta_2`
   - reuses the explicit input requirement from `R1`,
2. `orientation_slice_class_level_candidate_present`
   - reuses the class-level candidate from `C47`,
3. `basis_pair_export_skeleton_packet_ready`
   - reuses the minimal basis-pair export skeleton from `C48`,
4. `strict_core_actual_theta_export_absent`
   - reuses the negative strict-core verdict from `C35`,
5. `axiom_augmented_theta_source_branch_present_non_strict`
   - reuses the axiom-augmented source-branch verdict from `C35`.

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v2 :=
(
  rho_int_orientation_request_slot_v1,
  theta_like_input_source_sketch_v0
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v2 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v2,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not derive actual theta-like inputs.

It supplies something narrower:

1. the rho slot now knows exactly that the target-slot language needs two
   theta-like inputs,
2. the rho slot now knows that class-level orientation-slice and basis-pair
   scaffolds already exist,
3. the rho slot now also carries the exact split:
   strict-core actual theta export absent,
   axiom-augmented theta source branch present but non-strict.

So the slot is no longer only target-slot-aligned. It is now
target-slot-aligned plus theta-source-sketch-aware.

## What this still does not claim

`F02` does not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.
