# F03 Sandbox Non-Placeholder Strict-Core Theta-Source Skeleton Attempt Packet

Status: `F03_SANDBOX_NONPLACEHOLDER_STRICT_CORE_THETA_SOURCE_SKELETON_ATTEMPT_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Replace the purely sketch-level theta input note by one explicit
strict-core-only skeleton attempt that is richer than a placeholder but still
honestly below any real source discharge.

## Strict-core-only theta-source skeleton attempt

Define:

```text
theta_source_skeleton_attempt_v0 :=
(
  theta_formula_class_present,
  local_representative_class_present,
  orientation_slice_candidate_class_present,
  basis_pair_export_skeleton_present,
  actual_theta_source_rule_missing,
  populated_basis_pair_instance_missing
)
```

with the following meanings:

1. `theta_formula_class_present`
   - reuse `C33`:
     `theta_i = atan2(<s_i,u_i>, <c_i,u_i>)`,
2. `local_representative_class_present`
   - reuse `C34`:
     `u_i(theta_i)=cos(theta_i)c_i+sin(theta_i)s_i`,
3. `orientation_slice_candidate_class_present`
   - reuse `C47`:
     `S_orient_cand(theta_1,theta_2)=span{u_1(theta_1),u_2(theta_2)}`,
4. `basis_pair_export_skeleton_present`
   - reuse `C48`,
5. `actual_theta_source_rule_missing`
   - reuse the strict-core absence from `C50`,
6. `populated_basis_pair_instance_missing`
   - keep explicit that no actual populated `u_1,u_2` instance is exported.

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v3 :=
(
  rho_int_orientation_request_slot_v2,
  theta_source_skeleton_attempt_v0
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v3 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v3,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not produce a packet-ready strict-core theta source.

It produces something narrower:

1. one explicit strict-core-only skeleton attempt that names all formula/class
   ingredients already present on the strict side,
2. one equally explicit record of the two missing ingredients that still block
   actual sourcehood:
   - no actual theta-source rule,
   - no populated basis-pair instance.

So the sandbox no longer stops at a source sketch.
It now contains one non-placeholder strict-core theta-source skeleton attempt.

## What this still does not claim

`F03` does not claim:

1. actual `theta_1`, `theta_2`,
2. a packet-ready strict-core minimal source skeleton in the sense denied by
   `C50`,
3. actual populated `u_1`, `u_2`,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
