# P2 Strict-Core Sigma-Int To A1(pair1) Probe

Status: `P2_EXECUTED_STRICT_CORE_SIGMA_INT_TO_A1_PAIR1_COMPUTE_OR_FAIL_NO_FALSE_PASS`
As of: `2026-03-15`

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

Update (`2026-03-15`): the current route now **does** reach an operator-level
`A_1(pair1)` stage in the declared strict scope.

The report shows:

- `sigma_int_candidate` exists as a candidate object,
- residual `Z2` fit exists only on overlay/control lane,
- a strict-core sigma-int → residual export-map object exists in the declared `R1` scope (`T2`),
- strict-core `theta_1,theta_2` supply exists in the declared `pair1/pair2` scope (`C35`, `C49`),
- an actual populated basis-pair / `R1` inhabitant instance exists in the declared scope (`C48`, `C49`),
- and a strict-core operator object on `V_1 = span{c1,s1}` derived only from the materialized sigma-int `u_1` direction
  is exported (`F456`), providing the minimal downstream operator bridge target:
  `A_1(pair1) := |u_1><u_1|` (projector; residual `Z2` sign gauge invariant).

## Current compute-or-fail status

```text
PASS_COMPUTABLE_STRICT_CORE_SIGMA_INT_TO_A1_PAIR1_OPERATOR_STAGE_DECLARED_SCOPE
```

## Finite missing-object list

No upstream missing objects remain on the declared operator-stage scope.

## Honest frontier

- `T2` remains below any global selector closure and keeps `QW-2191` explicit.
- `F456` exports only a **minimal** strict-core operator object derived from `u_1` (a projector on `V_1`).
  It does not identify this operator with the extension-only `H/O` lane operator `A_1_ext` and does not
  claim any host matching/cancellation.
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

Update (`2026-03-15`): `P2` no longer blocks on the operator bridge itself; the next honest work is:

1. continue under explicit `QW-2191` discipline (no implied global discharge),
2. keep explicit that `A_1(pair1)` here is the minimal projector operator, not an identified `A_1_ext`.
