# F09 Sandbox Carrier-Directory Precondition Attack Packet

Status: `F09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Clear the single missing carrier-directory precondition identified by `F08`
without creating any provider file or provider instance.

## Executed additive move

The following directory now exists:

```text
fundamental_action_reconstruction/sandbox_strict_core_ingredient_attempt/generated/
```

This move is:

1. additive,
2. sandbox-local,
3. removable by deleting the sandbox subtree,
4. below any file-level provider creation,
5. below any provider emission.

## Carrier-directory attack package

Define:

```text
carrier_directory_precondition_attack_v0 :=
(
  target_carrier_directory_specified,
  target_carrier_directory_created,
  no_provider_file_created,
  no_provider_instance_emitted
)
```

with:

```text
target_carrier_directory_specified := yes
target_carrier_directory_created := yes
no_provider_file_created := yes
no_provider_instance_emitted := yes
```

## Creation-status field refinement after directory creation

Define:

```text
creation_status_refined_v1 := {
  filename_path_convention_specified,
  minimal_template_content_specified,
  carrier_directory_present,
  non_destructive_creation_admissibility,
  created_file_status
}
```

with the following values:

1. `filename_path_convention_specified := yes`
2. `minimal_template_content_specified := yes`
3. `carrier_directory_present := yes`
4. `non_destructive_creation_admissibility := yes_partial`
5. `created_file_status := file_not_created`

Meaning:
the directory-level precondition is now cleared, but file-level provider
creation still has not happened.

## Refined provider-object candidate

Define:

```text
PairIndexedStrictCoreThetaInputProviderCarrierCandidate_v2 := {
  semantic_name,
  pair_scope,
  target_outputs,
  filename_path_convention,
  minimal_template_content,
  creation_status,
  emission_status
}
```

with:

```text
creation_status := creation_status_refined_v1
```

and all other fields inherited from
`PairIndexedStrictCoreThetaInputProviderCarrierCandidate_v1`.

## Refined upstream provider schema

Define:

```text
StrictCoreThetaInputProviderArtifactSchema_v3 := {
  provider_object,
  target_outputs,
  downstream_consumer_contract,
  support_lane,
  strict_absence,
  fallback_contrast,
  forbidden_claims
}
```

with the refinement:

```text
provider_object := PairIndexedStrictCoreThetaInputProviderCarrierCandidate_v2
```

and all other fields inherited from `StrictCoreThetaInputProviderArtifactSchema_v2`.

## Refined upstream provider attack package

Define:

```text
strict_core_theta_input_provider_attack_v3 :=
(
  StrictCoreThetaInputProviderArtifactSchema_v3,
  carrier_directory_precondition_cleared,
  no_provider_file_created,
  provider_instance_absent
)
```

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v9 :=
(
  rho_int_orientation_request_slot_v8,
  strict_core_theta_input_provider_attack_v3
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v9 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v9,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not create a provider file.

It creates something narrower:

1. the missing carrier-directory precondition is now actually cleared,
2. the provider-carrier lane is now blocked only at file-creation and
   provider-emission level,
3. the old `directory absent` blocker no longer applies.

So the provider-carrier frontier is tighter without becoming a real created
provider carrier.

## What this still does not claim

`F09` does not claim:

1. actual provider file creation,
2. actual provider emission,
3. actual `theta_1`, `theta_2`,
4. actual populated `u_1`, `u_2`,
5. actual internal orientation datum,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. ToE closure.
