# F07 Sandbox Provider Object Field Refinement Packet

Status: `F07_SANDBOX_PROVIDER_OBJECT_FIELD_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Refine the single `provider_object` field from `F06` into one explicit
pair-indexed provider-carrier candidate.

## Provider-object field refinement

Define:

```text
PairIndexedStrictCoreThetaInputProviderCarrierCandidate_v0 := {
  semantic_name,
  pair_scope,
  target_outputs,
  filename_path_convention,
  minimal_template_content,
  creation_status,
  emission_status
}
```

with the following values:

1. `semantic_name`

   ```text
   semantic_name := strict_core_theta_input_provider_for_pair1_pair2
   ```

2. `pair_scope`

   ```text
   pair_scope := [pair1, pair2]
   ```

3. `target_outputs`

   ```text
   target_outputs := [theta_1, theta_2]
   ```

4. `filename_path_convention`

   ```text
   filename_path_convention :=
     sandbox_strict_core_ingredient_attempt/generated/
     strict_core_theta_input_provider_pair1_pair2_candidate.json
   ```

5. `minimal_template_content`

   ```json
   {
     "provider_object": "strict_core_theta_input_provider_for_pair1_pair2",
     "pair_scope": ["pair1", "pair2"],
     "target_outputs": ["theta_1", "theta_2"],
     "downstream_consumer_contract": "PopulatedBasisPairInstanceLayerArtifactSchema_v0",
     "status": "candidate_only"
   }
   ```

6. `creation_status`

   ```text
   creation_status := carrier_not_created
   ```

7. `emission_status`

   ```text
   emission_status := provider_not_emitted
   ```

## Refined upstream provider schema

Define:

```text
StrictCoreThetaInputProviderArtifactSchema_v1 := {
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
provider_object := PairIndexedStrictCoreThetaInputProviderCarrierCandidate_v0
```

and all other fields inherited from `StrictCoreThetaInputProviderArtifactSchema_v0`.

## Refined upstream provider attack package

Define:

```text
strict_core_theta_input_provider_attack_v1 :=
(
  StrictCoreThetaInputProviderArtifactSchema_v1,
  provider_object_field_refined,
  carrier_candidate_not_created,
  provider_instance_absent
)
```

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v7 :=
(
  rho_int_orientation_request_slot_v6,
  strict_core_theta_input_provider_attack_v1
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v7 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v7,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This refinement does not create an actual provider.

It creates something narrower:

1. the `provider_object` field is no longer only a generic label,
2. it is now one explicit pair-indexed provider-carrier candidate,
3. the candidate has:
   - semantic name,
   - pair scope,
   - output scope,
   - filename/path convention,
   - minimal template content,
   - explicit not-created / not-emitted status.

So the missing provider object is now semantically tighter without becoming a
real emitted provider.

## What this still does not claim

`F07` does not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual provider carrier creation,
4. actual provider emission,
5. actual internal orientation datum,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. ToE closure.
