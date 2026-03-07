# F4 Legacy Physical Role Transfer Classification Packet

Status: `F4_EXECUTED_LEGACY_PHYSICAL_ROLE_TRANSFER_CLASSIFICATION_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `K1/P47/N50/K2/F2/S2`, the repo already knows:

1. the legacy ontological kernel and the later strict gate kernel are split,
2. no rigorous bridge currently identifies them,
3. the strict gate kernel enters FAR only as a later-pipeline operational
   object unless a bridge is added.

`F4` asks the next honest classification question:

```text
which legacy physical roles may still be carried only by the legacy kernel
package, and which of them are not yet rigorously transferable to the
strict gate kernel on the current repo state?
```

## Result

`F4` establishes the following current-repo-state classification:

1. the old Weinberg-angle role
   `sin^2(theta_W) = alpha_geo / 12` remains a legacy-side
   `heuristic_not_strictly_derived_role`,
2. the old fine-structure role
   `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)` remains a legacy-side
   `partial_model_level_role`,
3. the old gravity-hierarchy role from `beta^N` scaling remains a legacy-side
   `model_consistent_but_not_full_independent_proof_role`,
4. none of those three legacy physical roles is currently rigorously
   transferred onto `K_strict_gate`,
5. the strict gate kernel remains operationally usable, but not as the
   silently inherited ontological carrier of those legacy roles.

## Classification rule

For the current repo state:

```text
legacy physical-role package
  = legacy kernel package plus downgraded legacy claim statuses

strict gate kernel package
  = later-pipeline operational kernel package
```

So the honest current rule is:

```text
do not silently inherit legacy Weinberg-angle, fine-structure, or
gravity-hierarchy roles into K_strict_gate
```

## Why this follows

`F4` is driven by the conjunction of:

1. `QW-2005`, which already downgrades the three legacy claims,
2. `P47/N50`, which already show the absence of a rigorous
   `K_legacy_ont -> K_strict_gate` bridge,
3. `F2`, which already classifies `K_strict_gate` as a later-pipeline
   operational import rather than the restored action-first ontological source
   layer.

## What F4 does not claim

`F4` does not claim:

- that the old legacy formulas are false,
- that the later strict pipeline is false,
- that no future bridge can ever transfer any role,
- that selector closure is solved,
- that `QW-2191` is discharged,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. test whether the repo already exports a rigorous transfer of any of those
   three legacy physical roles onto `K_strict_gate`,
2. or strengthen the current non-transfer frontier theorem-level,
3. without confusing operational strict-kernel outputs with theorem-level
   inheritance of the legacy role package.
