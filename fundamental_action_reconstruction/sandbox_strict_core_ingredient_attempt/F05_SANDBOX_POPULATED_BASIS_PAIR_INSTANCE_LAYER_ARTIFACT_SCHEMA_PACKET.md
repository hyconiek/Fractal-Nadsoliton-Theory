# F05 Sandbox Populated Basis-Pair Instance Layer Artifact Schema Packet

Status: `F05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ARTIFACT_SCHEMA_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Attack the missing populated basis-pair instance layer directly by replacing
its generic absence with one explicit packet-ready artifact schema.

## Packet-ready artifact schema for the missing instance layer

Define:

```text
PopulatedBasisPairInstanceLayerArtifactSchema_v0 := {
  object,
  required_inputs,
  population_rule,
  orientation_slice_role,
  conditional_theta_serialization_rule,
  current_absence,
  forbidden_claims
}
```

with the following values:

1. `object`

   ```text
   object := populated_basis_pair_instance(theta_1,theta_2)
   ```

2. `required_inputs`

   ```text
   required_inputs := [theta_1, theta_2]
   ```

3. `population_rule`

   ```text
   u_1 := cos(theta_1)c_1 + sin(theta_1)s_1
   u_2 := cos(theta_2)c_2 + sin(theta_2)s_2
   ```

4. `orientation_slice_role`

   ```text
   S_orient_cand := span{u_1,u_2}
   ```

5. `conditional_theta_serialization_rule`

   ```text
   if actual populated u_1,u_2 are present,
   then
     theta_1 = atan2(<s_1,u_1>, <c_1,u_1>)
     theta_2 = atan2(<s_2,u_2>, <c_2,u_2>)
   ```

6. `current_absence`

   ```text
   current_absence := [
     no_strict_core_supplied_actual_theta_1_theta_2,
     no_actual_populated_basis_pair_instance,
     no_packet_ready_strict_core_source_skeleton
   ]
   ```

7. `forbidden_claims`

   ```text
   forbidden_claims := [
     no_actual_basis_pair_export,
     no_actual_theta_export,
     no_actual_internal_orientation_datum,
     no_actual_E_orient,
     no_admissible_S_sel_int,
     no_strict_core_selector_closure,
     no_ToE_closure
   ]
   ```

## Direct layer attack package

Define:

```text
populated_basis_pair_instance_layer_attack_v0 :=
(
  PopulatedBasisPairInstanceLayerArtifactSchema_v0,
  artifact_schema_packet_ready,
  persisted_instance_absent,
  upstream_input_provider_absent
)
```

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v5 :=
(
  rho_int_orientation_request_slot_v4,
  populated_basis_pair_instance_layer_attack_v0
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v5 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v5,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not create an actual populated basis-pair instance.

It creates something narrower:

1. one explicit packet-ready artifact schema for the missing instance layer,
2. one explicit separation between:
   - artifact schema present,
   - actual populated instance absent,
   - upstream theta supplier absent,
3. one sharper direct statement of what the missing layer still needs before
   any strict-core instance export could occur.

So the sandbox now attacks the populated basis-pair instance layer directly,
not only the theta-source side feeding it.

## What this still does not claim

`F05` does not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual basis-pair export,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
