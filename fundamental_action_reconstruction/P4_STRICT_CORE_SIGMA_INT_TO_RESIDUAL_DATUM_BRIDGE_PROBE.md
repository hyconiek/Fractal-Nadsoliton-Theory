# P4 Strict-Core Sigma-Int To Residual Datum Bridge Probe

Status: `P4_EXECUTED_STRICT_CORE_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_PROBE_COMPUTE_OR_FAIL_NO_FALSE_PASS`
As of: `2026-03-15`

## Goal

After `P3` and `N6`, the active strict-core choke point is no longer the full
route to `theta` or `A_1(pair1)`.

The narrowest honest question is now:

```text
sigma_int_strict_derived_v1 -> residual_orientation_datum_target_slot (R1)
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

Update (`2026-03-15`): after exporting a slot-free strict-core sigma-int → theta-pair source (`F451/N489`) and an
audited `R1` target-slot inhabitant instance constructed from it (`P451`), the strict sigma-int lane is now computable
up to **target-slot population**.

Update (`2026-03-15`): after exporting an explicit post-witness object-support layer above the exported sign-only
export-map object (`F452/N490`), the strict sigma-int lane also discharges the post-`T148` object-support target
(`T130/N395`).

The report returns:

```text
PASS_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE
```

## What is present but still insufficient

The current repo does contain all of the following:

1. strict sigma-int provenance/value (`F307/N418`),
2. theorem-level gauge-quotient safety for that datum (`F308/N419`),
3. a strict-core target-slot export packet (`R1`),
4. an actual strict-core sign-only export-map object into that slot (`F311/N422`),
5. a slot-free strict-core sigma-int → theta-pair supply (`F451/N489`),
6. an audited inhabitant instance populating the `R1` target slot constructed from that theta-pair (`P451`),
7. an explicit axiom-lane bridge instance (outside strict core).

This is sufficient for strict-core **target-slot population**, and the post-`T148` object-support layer above the exported
map object is now exported (`F452/N490`), but it still does not imply selector closure (`QW-2191` remains open).

## Finite missing-object list

Update (`2026-03-15`): the theorem-level discharge of the conditional bridge theorem `T2` is now exported (`N491`).

The current strict-core route still lacks:

1. any strict-core selector closure / internal selector source discharging `QW-2191`.

## Honest frontier

- `C37` already gave candidate internalization only.
- `T2` already gave a conditional theorem spec only.
- `C46` already gave a carrier file only.
- `AX3` already gave a positive bridge instance only on the axiom lane.

So the present route reaches:

```text
target-slot population + post-map object-support (strict lane) + axiom-lane witness
```

but it does not yet export:

```text
any strict-core selector closure (QW-2191 remains open)
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

1. proceed to theorem-level discharge of `T2` without false pass, and/or
2. proceed explicitly on strict-core selector closure work under `QW-2191` discipline (no implied selector closure).
