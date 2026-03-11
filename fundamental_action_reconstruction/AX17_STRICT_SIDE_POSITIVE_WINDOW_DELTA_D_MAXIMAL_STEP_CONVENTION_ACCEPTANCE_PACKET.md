# AX17 Strict-Side Positive-Window Delta_d Maximal-Step Convention Acceptance Packet

Status: `AX17_EXECUTED_STRICT_SIDE_POSITIVE_WINDOW_DELTA_D_MAXIMAL_STEP_CONVENTION_ACCEPTANCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `T119`, the strict sigma-int → theta candidate pipeline admits a typed
positive-window corridor:

```text
delta_d ∈ (0, d_local/11].
```

After `P403/N437`, the repo now explicitly records that different admissible
choices of `delta_d` inside that corridor produce different theta-pair outputs.
Therefore `delta_d` is a real selector slot and is not uniquely fixed by the
corridor definition alone.

`AX17` performs the next honest move required to proceed with a *single*
reproducible strict-side theta-pair representative in a separated extension
scope, without promoting anything into strict core:

```text
accept one explicit convention for choosing delta_d on the strict sigma-int lane:
delta_d := delta_max = d_local/11.
```

## Inputs reused

1. `AX16`
   - strict-side extension scope is already accepted (`strict_extension_only`).
2. `T119`
   - positive-window corridor definition and its bound `delta_max := d_local/11`.
3. `F314/N425`
   - the current strict-input positive-window instantiation already uses the
     maximal-step convention `delta_d := delta_max`.
4. `P403/N437`
   - delta_d nonuniqueness/sensitivity is explicitly recorded (no false pass).

## What is accepted (extension scope only)

`AX17` accepts the following project/theory-level statement:

> In the `strict_extension_only` scope, whenever the strict sigma-int lane uses
> the `T119` positive-window corridor to generate a nad12 carrier for the
> sigma-int → theta pipeline, the corridor step is fixed by convention to:
> `delta_d := delta_max = d_local/11` computed from the strict working kernel tuple.

This is a convention/selector premise, not a strict-core derivation and not a
strict-core uniqueness result.

## Result of AX17

`AX17` establishes:

1. one explicit corridor-step selection rule is now fixed in a separated
   extension scope,
2. the strict sigma-int lane may cite `F314` / `F325` as the *chosen*
   representative theta-pair computation under that convention,
3. strict core remains unchanged.

## Hard limits

`AX17` does not claim:

1. strict-core derivation or uniqueness of `delta_d`,
2. actual strict-core `theta_1`, `theta_2` export,
3. admissible `S_sel_int` or strict-core selector closure,
4. `QW-2191` discharge,
5. ToE closure,
6. legacy-to-strict kernel bridge.

## Product

- one strict-side extension-scope convention packet fixing the `delta_d` choice,
- no false pass.

