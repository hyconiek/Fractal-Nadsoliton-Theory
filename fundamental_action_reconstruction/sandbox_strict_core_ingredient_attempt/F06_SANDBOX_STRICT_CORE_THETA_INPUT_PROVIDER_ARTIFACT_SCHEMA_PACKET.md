# F06 Sandbox Strict-Core Theta/Input Provider Artifact Schema Packet

Status: `F06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ARTIFACT_SCHEMA_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Attack the missing upstream strict-core theta/input provider directly by
replacing its generic absence with one explicit packet-ready artifact schema
bound to the downstream `F05` contract.

## Packet-ready artifact schema for the missing upstream provider

Define:

```text
StrictCoreThetaInputProviderArtifactSchema_v0 := {
  provider_object,
  target_outputs,
  downstream_consumer_contract,
  support_lane,
  strict_absence,
  fallback_contrast,
  forbidden_claims
}
```

with the following values:

1. `provider_object`

   ```text
   provider_object := strict_core_theta_input_provider_for_pair1_pair2
   ```

2. `target_outputs`

   ```text
   target_outputs := [theta_1, theta_2]
   ```

3. `downstream_consumer_contract`

   ```text
   downstream_consumer_contract := PopulatedBasisPairInstanceLayerArtifactSchema_v0
   ```

   Meaning:
   the provider would have to supply exactly the missing inputs required by
   `F05`.

4. `support_lane`

   ```text
   support_lane := strict_core_only_candidate_provider_lane
   ```

5. `strict_absence`

   ```text
   strict_absence := [
     no_packet_ready_strict_core_minimal_source_skeleton,
     no_actual_theta_1_theta_2_export,
     no_provider_instance_emitted
   ]
   ```

6. `fallback_contrast`

   ```text
   fallback_contrast := [
     axiom_augmented_theta_source_branch_exists,
     but_is_not_identified_with_strict_core_provider
   ]
   ```

7. `forbidden_claims`

   ```text
   forbidden_claims := [
     no_actual_strict_core_theta_provider,
     no_actual_theta_export,
     no_actual_basis_pair_export,
     no_actual_internal_orientation_datum,
     no_actual_E_orient,
     no_admissible_S_sel_int,
     no_strict_core_selector_closure,
     no_ToE_closure
   ]
   ```

## Direct upstream provider attack package

Define:

```text
strict_core_theta_input_provider_attack_v0 :=
(
  StrictCoreThetaInputProviderArtifactSchema_v0,
  artifact_schema_packet_ready,
  provider_instance_absent,
  downstream_contract_bound_to_F05
)
```

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v6 :=
(
  rho_int_orientation_request_slot_v5,
  strict_core_theta_input_provider_attack_v0
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v6 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v6,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not create an actual upstream strict-core theta/input
provider.

It creates something narrower:

1. one explicit packet-ready artifact schema for the missing provider layer,
2. one explicit binding from that provider layer to the downstream `F05`
   populated-instance contract,
3. one explicit separation between:
   - provider schema present,
   - provider instance absent,
   - axiom-augmented fallback present but not identified with strict core.

So the sandbox now attacks the upstream strict-core provider directly rather
than only its downstream consumer.

## What this still does not claim

`F06` does not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual basis-pair export,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
