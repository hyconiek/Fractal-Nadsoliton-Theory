# SESSION HANDOFF PROMPT AFTER RELEASE 5.4

Use this prompt to continue the FAR work in a new session without losing the
current guardrails or overstating the status.

## Role

Continue work as a strict, no-false-pass FAR coding/research agent inside this
repo.

You are not starting from scratch. You are continuing an already constrained
and partially constructive lane.

## First Read

Before editing `fundamental_action_reconstruction`, read:

1. `AGENTS.md`
2. `fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md`
3. `fundamental_action_reconstruction/K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md`
4. `fundamental_action_reconstruction/F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md`
5. `fundamental_action_reconstruction/F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md`
6. `fundamental_action_reconstruction/S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md`
7. `RELEASE_5_3.md`
8. `RELEASE_5_4.md`

Then inspect:

1. `fundamental_action_reconstruction/T14_SOURCE_TOPOLOGY_SELECTOR_THEOREM_SPEC.md`
2. `fundamental_action_reconstruction/F127_FIRST_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET.md`
3. `fundamental_action_reconstruction/F128_FIRST_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_PACKET.md`
4. `fundamental_action_reconstruction/F133_FIRST_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_PACKET.md`
5. `fundamental_action_reconstruction/F137_FIRST_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_TARGET_PACKET.md`
6. `fundamental_action_reconstruction/F138_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_PACKET.md`
7. `fundamental_action_reconstruction/P226_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_PROBE.md`
8. `fundamental_action_reconstruction/N246_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_THEOREM.md`
9. `fundamental_action_reconstruction/F139_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_PACKET.md`
10. `fundamental_action_reconstruction/P227_CURRENT_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_PROBE.md`
11. `fundamental_action_reconstruction/N247_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_THEOREM.md`
12. `fundamental_action_reconstruction/F140_FIRST_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_PACKET.md`
13. `fundamental_action_reconstruction/P228_CURRENT_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_PROBE.md`
14. `fundamental_action_reconstruction/N248_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_THEOREM.md`
15. `fundamental_action_reconstruction/F141_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_PACKET.md`
16. `fundamental_action_reconstruction/P229_CURRENT_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_PROBE.md`
17. `fundamental_action_reconstruction/N249_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_THEOREM.md`
18. `fundamental_action_reconstruction/F142_FIRST_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_PACKET.md`
19. `fundamental_action_reconstruction/P230_CURRENT_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_PROBE.md`
20. `fundamental_action_reconstruction/N250_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_THEOREM.md`
21. `fundamental_action_reconstruction/F143_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_PACKET.md`
22. `fundamental_action_reconstruction/P231_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_PROBE.md`
23. `fundamental_action_reconstruction/N251_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_THEOREM.md`
24. `fundamental_action_reconstruction/F144_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_PACKET.md`
25. `fundamental_action_reconstruction/P232_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_PROBE.md`
26. `fundamental_action_reconstruction/N252_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_COMPONENTS_PACKAGE_THEOREM.md`
27. `fundamental_action_reconstruction/F145_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_PACKET.md`
28. `fundamental_action_reconstruction/P233_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_PROBE.md`
29. `fundamental_action_reconstruction/N253_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_NONTRIVIALITY_ASSEMBLY_WITNESS_THEOREM.md`
30. `fundamental_action_reconstruction/F146_FIRST_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_PACKET.md`
31. `fundamental_action_reconstruction/P234_CURRENT_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_PROBE.md`
32. `fundamental_action_reconstruction/N254_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_THEOREM.md`
33. `fundamental_action_reconstruction/F147_FIRST_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_PACKET.md`
34. `fundamental_action_reconstruction/P235_CURRENT_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_PROBE.md`
35. `fundamental_action_reconstruction/N255_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_THEOREM.md`
36. `fundamental_action_reconstruction/F148_FIRST_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_PACKET.md`
37. `fundamental_action_reconstruction/P236_CURRENT_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_PROBE.md`
38. `fundamental_action_reconstruction/N256_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_THEOREM.md`
39. `fundamental_action_reconstruction/F149_FIRST_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_PACKET.md`
40. `fundamental_action_reconstruction/P237_CURRENT_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_PROBE.md`
41. `fundamental_action_reconstruction/N257_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_THEOREM.md`
42. `fundamental_action_reconstruction/F150_FIRST_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PACKET.md`
43. `fundamental_action_reconstruction/P238_CURRENT_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PROBE.md`
44. `fundamental_action_reconstruction/N258_CURRENT_FIRST_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM.md`
45. `fundamental_action_reconstruction/P239_CURRENT_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PROMOTION_PROBE.md`
46. `fundamental_action_reconstruction/N259_CURRENT_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PROMOTION_OBSTRUCTION_THEOREM.md`
47. `fundamental_action_reconstruction/P240_CURRENT_T14_DECLARED_SCOPE_COMPLETION_AND_CLOSURE_INCOMPLETENESS_PROBE.md`
48. `fundamental_action_reconstruction/N260_CURRENT_T14_DECLARED_SCOPE_COMPLETION_AND_CLOSURE_INCOMPLETENESS_THEOREM.md`
49. `fundamental_action_reconstruction/T15_LEGACY_TO_STRICT_KERNEL_BRIDGE_THEOREM_SPEC.md`
50. `fundamental_action_reconstruction/F151_FIRST_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_PACKET.md`
51. `fundamental_action_reconstruction/P241_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_PROBE.md`
52. `fundamental_action_reconstruction/N261_CURRENT_FIRST_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_THEOREM.md`

## Hard Guardrails

Keep these rules explicit:

1. Do not silently identify `K_legacy_ont` with `K_strict_gate`.
2. Do not silently transfer legacy physical-role claims onto the strict kernel.
3. Keep the ontology:
   `nadsoliton -> light -> matter -> emergent observer`.
4. Do not promote observer-side asymmetry into primary selector source.
5. Do not claim:
   - strict-core selector closure
   - `QW-2191` discharge
   - global selector closure
   - ToE closure
   unless explicitly proved.
6. Keep all future-route packets explicitly marked as future-route only.
7. No false PASS.

## Current Truthful State

The current repo state already has:

1. a positive constructive preobserver lane:
   `S_preLM_strict_core_source_object_v1 -> E_orient_preLM_v1 -> B_sel_preLM_v1 -> R_sel_preLM_v1 -> O_sel_preLM_v1`
2. a long downstream observer-side chain that remains downstream only
3. `N234` blocking any premature promotion of downstream observer stability into
   global selector closure
4. `T14` as the active source-topology selector route, now containing one
   declared-scope theorem witness below closure
5. the highest-priority post-`T14` frontier remains
   `legacy -> strict bridge or non-bridge`

The latest honest `T14` advance is:

```text
T14_src_selector_declared_scope_actual_witness_v1 :
tau_src_candidate_v1 -> declared_scope_source_topology_selector_theorem_target_v1
```

This is theorem-level packaged in `N258`, and `N259` now blocks any false
promotion of that theorem to current strict-core selector closure or current
global `QW-2191` discharge; `N260` freezes the lane as declared-scope complete
and closure-incomplete on the present export set.

## What Is Proven Right Now

Only this stronger-but-limited statement:

1. `tau_src_candidate_v1` contains one actual source-side scalar nonzero-flow
   component witness
2. `tau_src_candidate_v1` now has one actual source-side
   `source_limit_nonzero_flow_class_v1` witness
3. `tau_src_candidate_v1` contains one actual source-side scalar barrier-sign
   component witness with positive barrier margin
4. the declared core branch now has one actual positive-radius local
   barrier-sign stability witness
5. `tau_src_candidate_v1` now has one actual source-side
   `barrier_protected_sign_class_v1` witness
6. `tau_src_candidate_v1` now has one actual source-side
   `observer_free_scope_tag_v1` witness
7. the three actual witnesses are now bundled into one actual
   source-topology components package
8. that actual package now has one actual assembly witness into
   `Lambda_src_nontriv_target_v1`
9. `tau_src_candidate_v1` now has one actual full source-topology
   nontriviality witness
10. `tau_src_candidate_v1` now has one actual source-side selector witness
    into `Sigma_sel_src_target_v1`
11. `tau_src_candidate_v1` now has one actual source-side basis-independent
    selector-promotion witness into `Sigma_sel_basis_free_target_v1`
12. `tau_src_candidate_v1` now has one actual source-side quotient-safe
    `QW-2191` resolution witness in the declared source-topology scope
13. `tau_src_candidate_v1` now has one actual declared-scope Source Topology
    Selector theorem witness
14. these witnesses remain observer-free in their witness domain
15. the current repo does not justify promoting that theorem to current
    strict-core selector closure or to current global `QW-2191` discharge
16. the current `T14` lane is now frozen as declared-scope complete and
    closure-incomplete on the present export set
17. the current repo also exports one future-only positive bridge-branch target
    `B_legacy_strict_bridge_target_v1 : K_legacy_ont -> K_strict_gate`,
    while keeping the non-bridge branch explicit and leaving actual bridge
    derivation fully open
18. they remain below:
   - current global selector closure
   - global `QW-2191` discharge
   - strict-core selector closure
   - actual legacy-to-strict bridge derivation
   - ToE closure

## What Is Not Yet Proven

Do not claim any of the following:

1. actual global `QW-2191` discharge
2. actual strict-core selector closure
3. actual global selector closure
4. actual legacy-to-strict kernel bridge derivation
5. legacy physical-role transfer onto `K_strict_gate`

## Exact Next Honest Move

Preferred next move:

1. do not further promote the current `T14` lane on the present export set

Only after that:

2. if a stronger claim is still desired, work on the highest-priority
   `bridge or non-bridge` frontier rather than repackage `N258/N259/N260`

## Working Style

1. Use `apply_patch` for edits.
2. Keep commentary short and factual.
3. Prefer real constructive or theorem-spec steps over more meta-audit loops.
4. If a route is only future-route, say so explicitly in the packet status.
5. If a result is only component-level, do not promote it to full-route closure.

## If You Need a One-Line Summary

Current FAR state after Release 5.4:

```text
positive preobserver selector lane exists,
observer remains downstream only,
T14 now exports one declared-scope Source Topology Selector theorem,
and `N260` freezes that lane as declared-scope complete but closure-incomplete;
the repo also exports one future-only positive legacy-to-strict bridge target
while keeping the non-bridge branch open;
the newest honest source-topology state remains below
global `QW-2191` discharge and strict-core selector closure.
```
