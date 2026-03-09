# F17 Sandbox Actual Populated Basis-Pair Instance Semantic Blocker Attack Packet

Status: `F17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Attack the populated-instance blocker directly by replacing the flat absence
label with one explicit populated-instance gate audit and one dedicated
candidate file below actual population.

## Executed additive move

The following removable sandbox file now exists:

```text
fundamental_action_reconstruction/sandbox_strict_core_ingredient_attempt/generated/strict_core_populated_basis_pair_instance_candidate_pair1_pair2.json
```

It records one candidate-only populated-instance template below actual
population.

## Semantic gate refinement

Define:

```text
actual_populated_basis_pair_instance_semantic_gate_v0 := {
  minimal_basis_pair_export_skeleton_present,
  conditional_populated_instance_schema_present,
  dedicated_population_candidate_file_present,
  theta_input_slots_written,
  actual_theta_input_values_present,
  actual_populated_basis_pair_instance_present,
  instance_gate_verdict
}
```

with:

1. `minimal_basis_pair_export_skeleton_present := yes`

   Reason:
   `C48` already exports the minimal basis-pair export skeleton.

2. `conditional_populated_instance_schema_present := yes`

   Reason:
   `C49` already exports the conditional populated-instance schema.

3. `dedicated_population_candidate_file_present := yes`

   Reason:
   the new sandbox candidate JSON file now exists.

4. `theta_input_slots_written := yes`

   Reason:
   the candidate file explicitly contains `theta_1` and `theta_2` slots.

5. `actual_theta_input_values_present := no`

   Reason:
   both `theta_1` and `theta_2` remain `null`.

6. `actual_populated_basis_pair_instance_present := no`

   Reason:
   without actual theta inputs, the candidate file still does not populate
   actual `u_1`, `u_2`.

7. `instance_gate_verdict := candidate_template_present_below_actual_theta_inputs`

## Populated-instance interpretation

This means the populated-instance blocker is no longer best described as:

```text
total absence of populated basis-pair instance layer
```

It is now more sharply described as:

```text
the populated-instance lane is structurally present below actual theta inputs,
and the live negative point is the absence of actual theta values
```

## Candidate file refinements

The populated-instance candidate file now stores:

1. required inputs,
2. formulas for `u_1`, `u_2`,
3. explicit `theta_1`, `theta_2` slots,
4. explicit `u_1`, `u_2` slots,
5. `population_status := candidate_template_only_not_populated`.

The provider candidate file now records the same boundary under:

```text
actual_populated_basis_pair_instance_status
```

with:

1. `minimal_basis_pair_export_skeleton_present := true`
2. `conditional_populated_instance_schema_present := true`
3. `dedicated_population_candidate_file_present := true`
4. `theta_input_slots_written := true`
5. `actual_theta_input_values_present := false`
6. `actual_populated_basis_pair_instance_present := false`
7. `instance_gate_verdict := candidate_template_present_below_actual_theta_inputs`

## Refined provider attack package

Define:

```text
strict_core_theta_input_provider_attack_v11 :=
(
  StrictCoreThetaInputProviderArtifactSchema_v5,
  provider_emission_gate_audit_v0,
  provider_emission_four_clause_attack_bundle_v0,
  theta_output_clause_positive_support_v0,
  strict_to_axiom_theta_source_bridge_artifact_schema_v1,
  theta_source_supply_boundary_direct_attack_v0,
  bridge_discharge_status_refined_v0,
  actual_strict_core_theta_source_supply_semantic_gate_v0,
  actual_populated_basis_pair_instance_semantic_gate_v0,
  provider_emission_admissibility_failed,
  provider_not_emitted
)
```

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v17 :=
(
  rho_int_orientation_request_slot_v16,
  strict_core_theta_input_provider_attack_v11
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v17 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v17,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

`F17` does not populate an actual basis-pair instance.

It still makes one real narrowing move:

1. the populated-instance blocker is no longer a flat `false`,
2. it now has one explicit gate audit,
3. a dedicated populated-instance candidate file now exists,
4. the live negative point is sharply isolated at:
   actual theta input values are absent,
5. therefore actual populated instance still remains absent.

## Loop visibility

This attack also makes the dependency loop explicit:

1. theta-source supply was narrowed down to missing populated instance,
2. populated-instance realization is now narrowed down to missing actual theta
   inputs.

The sandbox therefore exposes one real circular frontier rather than hiding
it.

## What this still does not claim

`F17` does not claim:

1. actual populated basis-pair instance,
2. actual `u_1`, `u_2`,
3. actual `theta_1`, `theta_2`,
4. actual strict-core theta source supply,
5. actual bridge discharge,
6. actual provider emission,
7. actual `E_orient`,
8. admissible `S_sel_int`,
9. strict-core selector closure,
10. ToE closure.
