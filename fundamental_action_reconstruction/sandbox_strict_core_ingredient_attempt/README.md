# Sandbox Strict-Core Ingredient Attempt

Status: `SANDBOX_ONLY_REMOVABLE_STRICT_CORE_INGREDIENT_ATTEMPT`
As of: `2026-03-09`

## Purpose

This folder is a removable sandbox for one explicit attempt to build a
genuinely new strict-core ingredient instead of continuing the
`strict_extension_only` clause-lift ladder.

It is intentionally isolated from the official theorem lane:

1. no files in this folder are referenced from the main `T/F/P/N` chain,
2. no files in this folder are referenced from the top-level release notes,
3. no claim in this folder counts as current strict-core discharge.

If this route turns out to be unproductive, the whole attempt may be removed
simply by deleting this folder.

## Files

1. `T00_SANDBOX_STRICT_CORE_INGREDIENT_ATTEMPT_SCOPE.md`
   - scope and guardrails of the sandbox attempt,
2. `F00_SANDBOX_STRICT_CORE_INGREDIENT_CANDIDATE_PACKET.md`
   - one explicit sandbox-only candidate object,
3. `P00_SANDBOX_STRICT_CORE_INGREDIENT_CANDIDATE_PROBE.md`
   - clause-by-clause diagnostic against the `F29` admission contract,
4. `N00_SANDBOX_STRICT_CORE_INGREDIENT_STATUS_NOTE.md`
   - the strongest honest current reading of the sandbox attempt.
5. `T01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_SCOPE.md`
   - repo-consistent scope for attacking the rho orientation slot,
6. `F01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_PACKET.md`
   - refinement of the rho slot into a target-slot-aligned request scaffold,
7. `P01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_PROBE.md`
   - check that the rho refinement stays below discharge,
8. `N01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the rho-slot attack.
9. `T02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_SCOPE.md`
   - scope for attacking theta-like inputs for the refined rho slot,
10. `F02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_PACKET.md`
   - refinement of the rho slot into a theta-source-sketch-aware scaffold,
11. `P02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_PROBE.md`
   - check that the theta-like input refinement stays below discharge,
12. `N02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the theta-like input attack.
13. `T03_SANDBOX_NONPLACEHOLDER_STRICT_CORE_THETA_SOURCE_SKELETON_ATTEMPT_SCOPE.md`
   - scope for writing a strict-core-only theta-source skeleton attempt,
14. `F03_SANDBOX_NONPLACEHOLDER_STRICT_CORE_THETA_SOURCE_SKELETON_ATTEMPT_PACKET.md`
   - explicit non-placeholder strict-core theta-source skeleton attempt,
15. `P03_SANDBOX_NONPLACEHOLDER_STRICT_CORE_THETA_SOURCE_SKELETON_ATTEMPT_PROBE.md`
   - check that the theta-source attempt stays below actual sourcehood,
16. `N03_SANDBOX_NONPLACEHOLDER_STRICT_CORE_THETA_SOURCE_SKELETON_ATTEMPT_STATUS_NOTE.md`
   - strongest honest status after the non-placeholder theta-source attempt.
17. `T04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_SCOPE.md`
   - scope for writing one narrower conditional strict-core theta-source rule
     candidate above the skeleton attempt,
18. `F04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_PACKET.md`
   - explicit conditional strict-core theta-source rule candidate,
19. `P04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_PROBE.md`
   - check that the rule candidate sharpens the sandbox without contradicting
     `C50`,
20. `N04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_STATUS_NOTE.md`
   - strongest honest status after the conditional rule-candidate attack.
21. `T05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ATTACK_SCOPE.md`
   - scope for a direct attack on the missing populated basis-pair instance
     layer,
22. `F05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ARTIFACT_SCHEMA_PACKET.md`
   - packet-ready artifact schema for the missing populated basis-pair
     instance layer,
23. `P05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ARTIFACT_SCHEMA_PROBE.md`
   - check that the layer attack is sharper than `F04` while remaining below
     actual instance export,
24. `N05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ARTIFACT_SCHEMA_STATUS_NOTE.md`
   - strongest honest status after the direct populated-instance-layer attack.
25. `T06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ATTACK_SCOPE.md`
   - scope for attacking the upstream strict-core theta/input provider with
     `F05` as downstream contract,
26. `F06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ARTIFACT_SCHEMA_PACKET.md`
   - packet-ready artifact schema for the missing upstream strict-core
     theta/input provider,
27. `P06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ARTIFACT_SCHEMA_PROBE.md`
   - check that the provider attack is sharper than `F05` while remaining
     below actual theta export,
28. `N06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ARTIFACT_SCHEMA_STATUS_NOTE.md`
   - strongest honest status after the direct upstream provider attack.
29. `T07_SANDBOX_PROVIDER_OBJECT_FIELD_REFINEMENT_SCOPE.md`
   - scope for refining the single `provider_object` field from `F06`,
30. `F07_SANDBOX_PROVIDER_OBJECT_FIELD_REFINEMENT_PACKET.md`
   - pair-indexed provider-carrier candidate refining the `provider_object`
     field,
31. `P07_SANDBOX_PROVIDER_OBJECT_FIELD_REFINEMENT_PROBE.md`
   - check that the field refinement is stronger than the generic `provider_object`
     label while remaining below provider emission,
32. `N07_SANDBOX_PROVIDER_OBJECT_FIELD_REFINEMENT_STATUS_NOTE.md`
   - strongest honest status after refining the `provider_object` field.
33. `T08_SANDBOX_PROVIDER_CREATION_STATUS_FIELD_ATTACK_SCOPE.md`
   - scope for attacking the single `creation_status` field from `F07`,
34. `F08_SANDBOX_PROVIDER_CREATION_STATUS_FIELD_REFINEMENT_PACKET.md`
   - sandbox refinement of the `creation_status` field into an explicit
     readiness/non-readiness verdict,
35. `P08_SANDBOX_PROVIDER_CREATION_STATUS_FIELD_REFINEMENT_PROBE.md`
   - check that the refinement is stronger than the coarse `carrier_not_created`
     label while remaining below carrier creation,
36. `N08_SANDBOX_PROVIDER_CREATION_STATUS_FIELD_REFINEMENT_STATUS_NOTE.md`
   - strongest honest status after refining the `creation_status` field.
37. `T09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_SCOPE.md`
   - scope for directly attacking the missing carrier-directory precondition,
38. `F09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_PACKET.md`
   - packet recording carrier-directory creation and the resulting refinement
     of the provider-carrier status,
39. `P09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_PROBE.md`
   - check that the directory precondition is now cleared while file/provider
     creation remains absent,
40. `N09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_STATUS_NOTE.md`
   - strongest honest status after clearing the carrier-directory precondition.
41. `T10_SANDBOX_PROVIDER_FILE_CREATION_LAYER_ATTACK_SCOPE.md`
   - scope for a real file-level carrier step on the provider lane,
42. `F10_SANDBOX_PROVIDER_FILE_CREATION_LAYER_ATTACK_PACKET.md`
   - packet recording creation of one minimal provider candidate file,
43. `P10_SANDBOX_PROVIDER_FILE_CREATION_LAYER_ATTACK_PROBE.md`
   - check that file creation occurred while provider emission still did not,
44. `N10_SANDBOX_PROVIDER_FILE_CREATION_LAYER_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the provider file creation attack.
45. `T11_SANDBOX_PROVIDER_EMISSION_LAYER_ATTACK_SCOPE.md`
   - scope for directly attacking the provider-emission layer above the
     created candidate file,
46. `F11_SANDBOX_PROVIDER_EMISSION_GATE_AUDIT_PACKET.md`
   - audit packet refining `emission_status` into an explicit emission-gate
     verdict,
47. `P11_SANDBOX_PROVIDER_EMISSION_LAYER_ATTACK_PROBE.md`
   - check that the provider file exists but still fails the emission gate,
48. `N11_SANDBOX_PROVIDER_EMISSION_LAYER_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the provider-emission attack.
49. `T12_SANDBOX_ALL_FOUR_PROVIDER_EMISSION_FAILURE_CLAUSES_ATTACK_SCOPE.md`
   - scope for attacking all four live provider-emission failure clauses in
     one batch,
50. `F12_SANDBOX_ALL_FOUR_PROVIDER_EMISSION_FAILURE_CLAUSES_ATTACK_PACKET.md`
   - packet refining all four live emission-failure clauses into one explicit
     structured failure bundle,
51. `P12_SANDBOX_ALL_FOUR_PROVIDER_EMISSION_FAILURE_CLAUSES_ATTACK_PROBE.md`
   - check that all four clauses were sharpened while provider emission still
     remains inadmissible,
52. `N12_SANDBOX_ALL_FOUR_PROVIDER_EMISSION_FAILURE_CLAUSES_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the all-four-clause emission failure
     attack.
53. `T13_SANDBOX_THETA_OUTPUT_FAILURE_CLAUSE_POSITIVE_ATTACK_SCOPE.md`
   - scope for one positive attack on the theta-output failure clause,
54. `F13_SANDBOX_THETA_OUTPUT_FAILURE_CLAUSE_POSITIVE_ATTACK_PACKET.md`
   - packet extracting the maximal positive theta-output support available
     below actual values,
55. `P13_SANDBOX_THETA_OUTPUT_FAILURE_CLAUSE_POSITIVE_ATTACK_PROBE.md`
   - check that the theta-output clause now has positive support while actual
     theta values remain absent,
56. `N13_SANDBOX_THETA_OUTPUT_FAILURE_CLAUSE_POSITIVE_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the positive theta-output clause attack.
57. `T14_SANDBOX_STRICT_CORE_THETA_SOURCE_SUPPLY_BOUNDARY_DIRECT_ATTACK_SCOPE.md`
   - scope for directly attacking the strict-core theta-source supply
     boundary,
58. `F14_SANDBOX_STRICT_TO_AXIOM_THETA_SOURCE_BRIDGE_ARTIFACT_SCHEMA_PACKET.md`
   - packet assembling one explicit strict-to-axiom theta-source bridge
     artifact schema below discharge,
59. `P14_SANDBOX_STRICT_CORE_THETA_SOURCE_SUPPLY_BOUNDARY_DIRECT_ATTACK_PROBE.md`
   - check that the boundary now has an assembled bridge artifact schema while
     actual strict-core source supply remains absent,
60. `N14_SANDBOX_STRICT_CORE_THETA_SOURCE_SUPPLY_BOUNDARY_DIRECT_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the direct theta-source supply boundary
     attack.
61. `T15_SANDBOX_STRICT_TO_AXIOM_BRIDGE_DISCHARGE_STATUS_ATTACK_SCOPE.md`
   - scope for directly attacking the `discharge_status` field of the
     assembled bridge schema,
62. `F15_SANDBOX_STRICT_TO_AXIOM_BRIDGE_DISCHARGE_STATUS_ATTACK_PACKET.md`
   - packet refining bridge discharge into one explicit gate audit and adding
     one dedicated candidate carrier file,
63. `P15_SANDBOX_STRICT_TO_AXIOM_BRIDGE_DISCHARGE_STATUS_ATTACK_PROBE.md`
   - check that discharge status is now explicit while actual discharge
     remains blocked,
64. `N15_SANDBOX_STRICT_TO_AXIOM_BRIDGE_DISCHARGE_STATUS_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the bridge discharge-status attack.
65. `T16_SANDBOX_ACTUAL_STRICT_CORE_THETA_SOURCE_SUPPLY_SEMANTIC_BLOCKER_ATTACK_SCOPE.md`
   - scope for directly attacking the final semantic blocker under the bridge
     discharge lane,
66. `F16_SANDBOX_ACTUAL_STRICT_CORE_THETA_SOURCE_SUPPLY_SEMANTIC_BLOCKER_ATTACK_PACKET.md`
   - packet refining the last semantic blocker into one explicit source-supply
     gate audit plus one dedicated candidate file,
67. `P16_SANDBOX_ACTUAL_STRICT_CORE_THETA_SOURCE_SUPPLY_SEMANTIC_BLOCKER_ATTACK_PROBE.md`
   - check that the final semantic blocker is narrowed while actual theta
     supply still remains absent,
68. `N16_SANDBOX_ACTUAL_STRICT_CORE_THETA_SOURCE_SUPPLY_SEMANTIC_BLOCKER_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the final semantic-blocker attack.
69. `T17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_SCOPE.md`
   - scope for directly attacking the populated-instance blocker under the
     theta-source lane,
70. `F17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_PACKET.md`
   - packet refining the populated-instance blocker into one explicit gate
     audit plus one dedicated candidate file,
71. `P17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_PROBE.md`
   - check that the populated-instance blocker is narrowed while actual
     population still remains absent,
72. `N17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_STATUS_NOTE.md`
   - strongest honest status after the populated-instance blocker attack.
73. `T18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_SCOPE.md`
   - scope for freezing the exposed theta-supply / population loop as a
     current-state incompatibility boundary,
74. `F18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_PACKET.md`
   - packet formalizing the loop boundary on present sandbox state and present
     strict-core inputs,
75. `P18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_PROBE.md`
   - check that the current route is nonentering under the same blocker-cut,
76. `N18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_STATUS_NOTE.md`
   - strongest honest current-state incompatibility boundary for the sandbox
     route.

## Hard limits

This folder does not prove:

1. admissible `S_sel_int`,
2. actual `E_orient`,
3. actual strict-core selector closure,
4. actual global `QW-2191` discharge,
5. actual ToE closure.

It exists only to make one alternative positive strict-core construction route
explicit and easy to delete if it stalls.
