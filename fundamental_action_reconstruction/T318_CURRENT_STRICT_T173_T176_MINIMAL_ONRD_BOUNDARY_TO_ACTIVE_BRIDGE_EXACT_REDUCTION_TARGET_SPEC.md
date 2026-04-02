# T318 Current Strict `T173/T176` Minimal ONRD Boundary To Active Bridge Exact Reduction Target Spec

Status: `TARGET_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_NO_FALSE_PASS`
As of: `2026-03-29`

## Target

Export exactly one target object:

```text
MinimalONRDBoundaryToActiveBridgeExactReductionTarget_v1
```

with intended meaning:

```text
derive one lawful exact reduction,
from MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1
into ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1
```

## Boundary Contract

The target remains honest only if all of the following are preserved:

1. source boundary object remains
   `MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1`,
2. source boundary grade remains
   `candidate_provider_class_seed_only`,
3. target bridge remains
   `ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1`,
4. exact reduction is still missing,
5. no lawful supplier is claimed,
6. no solution is claimed,
7. no strict physical orientation datum is claimed.

## Hard Limits

`T318` does **not** itself export the reduction, supplier, solution, or closure.
