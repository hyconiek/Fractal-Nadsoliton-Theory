# P1095 Current Strict `T173/T176` Post-`F969` Minimal ONRD Boundary To Active Bridge Exact Reduction Target Admission Probe

Status: `P1095_CURRENT_STRICT_T173_T176_POST_F969_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ADMISSION_PROBE_NO_FALSE_PASS`
As of: `2026-03-29`

## Goal

Freeze the next honest move after `P1094/F969`:

```text
admit one exact reduction target
from MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1
into ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1,
without claiming that the reduction already exists.
```

`ONRD := oriented nonreciprocal dephasing`.

## Scope

`P1095` is a target-admission probe only.

It does **not**:

1. export the reduction,
2. export a lawful supplier,
3. export a solution,
4. export a strict physical orientation datum,
5. discharge `T183`, `T176`, `QW-2191`, or ToE closure.

## Acceptance

The probe passes only if it shows all of the following:

1. `F966` still allows one genuinely new inversion-sensitive source-side provider-class route,
2. `F969` exports the minimal ONRD boundary only as `candidate_provider_class_seed_only`,
3. the boundary and the active bridge point to the same frontier context,
4. exact reduction from that boundary to the active bridge is still not exported,
5. old same-lane reentry remains disallowed,
6. and therefore the strongest honest export is one exact reduction target, not a supplier and not a solution.

## Hard Limits

`P1095` does **not** claim:

1. that the boundary already realizes the active bridge,
2. that the reduction already exists,
3. that the route is already lawful as a supplier,
4. that orientation selection is already solved.
