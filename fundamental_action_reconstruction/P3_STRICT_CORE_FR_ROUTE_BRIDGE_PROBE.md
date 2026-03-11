# P3 Strict-Core FR-Route Bridge Probe

Status: `P3_EXECUTED_STRICT_CORE_FR_ROUTE_BRIDGE_COMPUTE_OR_FAIL_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

The current FR/topological route is the only remaining serious internal-source
candidate.

`P3` tests the narrow bridge question:

```text
sigma_int_candidate
  -> residual orientation datum
  -> theta-source
```

The probe is `compute-or-fail`:

- either the current strict-core FR route already exports the bridge objects,
- or it returns an explicit finite missing-object list.

## Result

The current route does **not** reach a strict-core bridge.

The report returns:

```text
NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_FR_ROUTE
```

## Finite missing-object list

The current FR route still lacks:

1. a strict derivation or strict-core source-object upgrade for
   `sigma_int_candidate`,
2. theorem-level gauge-quotient safety for `sigma_int_candidate`,
3. a strict-core equivalence/export map
   `sigma_int_candidate -> residual orientation datum`,
4. a strict-core selector map
   `sigma_int_candidate -> theta`
   or an internal derivation of the `J_ab` selector family,
5. a strict-core actual phase source for `theta_1`, `theta_2`.

## Honest frontier

- `B8` already exposes the FR-route residual blockers.
- `T2_B1` remains active:
  strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent.
- `C35_B1` remains active:
  only an axiom-augmented actual-theta source branch exists.

## What P3 does not claim

`P3` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that the FR route is false,
- that no future strict-core FR bridge can exist,
- that `QW-2191` is discharged.

## Recommended next move

Only two serious routes remain:

1. construct one missing strict-core bridge object and rerun `P3`,
2. or upgrade the negative theorem lane for the FR route itself.
