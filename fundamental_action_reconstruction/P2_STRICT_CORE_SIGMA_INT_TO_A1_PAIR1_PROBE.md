# P2 Strict-Core Sigma-Int To A1(pair1) Probe

Status: `P2_EXECUTED_STRICT_CORE_SIGMA_INT_TO_A1_PAIR1_COMPUTE_OR_FAIL_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N5`, the narrowest executable question is no longer whether the current
strict-core `psi0` route is obstructed in general, but whether the best current
strict-core internal-source candidate can already be pushed all the way to an
operator on `pair1`.

`P2` therefore tests one concrete route:

```text
sigma_int_candidate
  -> residual orientation datum
  -> actual theta_1, theta_2
  -> actual basis pair u_1, u_2 / orientation slice
  -> A_1(pair1)
```

The probe is `compute-or-fail`:

- either the route is sufficiently exported to continue toward `A_1(pair1)`,
- or the report returns an explicit finite missing-object list.

## Inputs

Strict-core route state is read from:

- `B2`
- `B4`
- `B6`
- `B7`
- `T2`
- `C35`
- `C37`
- `C46`
- `C47`
- `C48`
- `C49`
- `QW-2191`

## Result

The current route does **not** reach `A_1(pair1)`.

The report shows:

- `sigma_int_candidate` exists as a candidate object,
- residual `Z2` fit exists only on overlay/control lane,
- strict-core equivalence/export map remains absent,
- strict-core actual `theta_1`, `theta_2` supply remains absent,
- basis-pair export remains only skeleton/conditional-schema level,
- so the route stops before any strict-core operator-level export on `pair1`.

## Current compute-or-fail status

```text
NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE
```

## Finite missing-object list

The probe reduces the blocker to the following current missing objects:

1. a strict-core source object upgrading `sigma_int_candidate` beyond
   candidate-only status,
2. a strict-core equivalence/export map
   `sigma_int_candidate -> residual orientation datum`,
3. a strict-core actual phase source for `theta_1`, `theta_2` on the current
   pair frames,
4. a populated actual basis-pair export `u_1`, `u_2`,
5. a strict-core operator-level bridge from the materialized orientation slice
   to `A_1(pair1)`.

## Honest frontier

- `T2_B1` remains active:
  strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent.
- `C35_B1` remains active:
  strict core does not export actual `theta_1`, `theta_2`; only an
  axiom-augmented source branch exists.
- `C49_B1` remains active:
  the conditional populated-instance schema cannot be instantiated because
  strict core still does not supply actual `theta_1`, `theta_2`.
- `C32_B2` remains active:
  the raw overlap route remains degenerate.

## What P2 does not claim

`P2` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `sigma_int_candidate` is false,
- that no future strict-core route can ever reach `A_1(pair1)`,
- that `QW-2191` is discharged.

## Recommended next move

Only two serious routes remain:

1. add one real strict-core source/bridge object that resolves one item from
   the missing-object list and rerun `P2`,
2. or upgrade the negative lane with a stronger theorem showing why one or more
   of those missing objects cannot arise inside the current strict core.
