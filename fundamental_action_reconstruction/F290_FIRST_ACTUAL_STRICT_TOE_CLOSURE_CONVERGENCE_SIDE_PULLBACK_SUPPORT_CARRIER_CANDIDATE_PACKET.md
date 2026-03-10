# F290 First Actual Strict ToE Closure Convergence-Side Pullback Support Carrier Candidate Packet

Status: `F290_EXECUTED_FIRST_ACTUAL_STRICT_TOE_CLOSURE_CONVERGENCE_SIDE_PULLBACK_SUPPORT_CARRIER_CANDIDATE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

After `T134/F289/N401`, the repo exports one joint coherence witness candidate.
The next honest convergence-side move (still below `N302`) is to package one
explicit **pullback carrier** candidate that:

1. binds the provider projection record (`N399`) and the sigma-int projection
   record (`N400`) on the same pair family,
2. remains noncyclic and observer-free,
3. remains selector-neutral.

`F290` executes exactly that packaging step.

## Inputs reused

1. `N399`
   - provider-object carrier → residual projection candidate (positive-window).
2. `N400`
   - sigma-int → residual projection candidate (positive-window low-delta).
3. `T135`
   - pullback support carrier candidate spec.
4. `N302`
   - boundary below actual bridge/export-map object support remains in force.

## Pullback carrier candidate result

`F290` exports one packaged pullback support carrier candidate:

```text
Mu_strict_convergence_side_pullback_support_carrier_candidate_v1
```

materialized by the persisted artifact instance:

```text
fundamental_action_reconstruction/generated/strict_convergence_side_pullback_support_carrier_candidate_instance.json
```

with the intended reading:

1. the pullback carrier is candidate-only and record-level,
2. it packages two already-exported projection records on the same pair family,
3. it introduces no theta inputs (thetas appear only as candidate outputs),
4. it introduces no populated-instance inputs,
5. it introduces no `K_obs`-indexed selection,
6. it remains selector-neutral,
7. it remains strictly below actual object support (`N302`).

## Status discipline

This packet does **not** claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual export-map object export,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

