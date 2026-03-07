# P4 Strict-Core Sigma-Int To Residual Datum Bridge Probe

Status: `P4_EXECUTED_STRICT_CORE_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_PROBE_COMPUTE_OR_FAIL_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P3` and `N6`, the active strict-core choke point is no longer the full
route to `theta` or `A_1(pair1)`.

The narrowest honest question is now:

```text
sigma_int_candidate -> residual orientation datum
```

`P4` tests that bridge in `compute-or-fail` mode:

- either the current strict core already exports a real bridge,
- or the probe returns an explicit finite missing-object list.

## Why this probe is narrower than `P3`

`P3` asked whether the current FR/topological route reaches:

```text
sigma_int_candidate -> residual orientation datum -> theta-source
```

`P4` removes the downstream `theta` question and asks only whether the bridge
to the residual datum itself exists in strict core.

That is the first common upstream blocker for both:

- the FR route,
- the longer `sigma_int -> ... -> A_1(pair1)` route.

## Result

The current route does **not** reach a strict-core residual-datum bridge.

The report returns:

```text
NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE
```

## What is present but still insufficient

The current repo does contain all of the following:

1. `sigma_int_candidate`,
2. candidate-fit of `sigma_int_candidate` to the residual `Z2` slot,
3. semantic target-slot fields in the acceptance carrier,
4. a persisted acceptance-artifact carrier file,
5. an explicit axiom-lane bridge instance.

But these still do **not** amount to a strict-core bridge.

## Finite missing-object list

The current strict-core route still lacks:

1. a strict derivation or source-object upgrade for `sigma_int_candidate`,
2. theorem-level gauge-quotient safety for `sigma_int_candidate`,
3. a strict-core exported target slot for the residual orientation datum,
4. a strict-core equivalence/export map
   `sigma_int_candidate -> residual orientation datum`,
5. selector-track identification beyond overlay-only compatibility.

## Honest frontier

- `C37` already gave candidate internalization only.
- `T2` already gave a conditional theorem spec only.
- `C46` already gave a carrier file only.
- `AX3` already gave a positive bridge instance only on the axiom lane.

So the present route stops at:

```text
candidate-fit + carrier + axiom-lane witness
```

and does not cross into:

```text
strict-core residual datum bridge
```

## What P4 does not claim

`P4` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that the FR/topological idea is false,
- that no future strict-core bridge can exist,
- that `QW-2191` is discharged,
- that the axiom-lane witness can be promoted into strict core.

## Recommended next move

Only two serious routes remain:

1. construct one missing strict-core bridge object for
   `sigma_int_candidate -> residual orientation datum`,
2. or formalize the route-specific negative theorem for this exact bridge lane.
