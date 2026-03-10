# F287 First Actual Strict ToE Closure Provider-Object Carrier to Residual Bridge/Export-Map Object-Support Projection Candidate (Positive-Window) Packet

Status: `F287_EXECUTED_FIRST_ACTUAL_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_TO_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_PACKET_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T127/F280/N392` exports one explicit provider-carrier → residual projection
candidate, but it remains exposed to the explicit degeneracy frontier named in
`T115`:

```text
(X_i^cand,Y_i^cand)=(0,0) => atan2 undefined
```

`F287` packages a strictly weaker but strictly safer variant:

```text
positive-window corridor projection candidate
```

that enforces a typed `d`-domain restriction guaranteeing:

```text
cos(Theta(d))>0 and sin(Theta(d))>0 on all paths
=> X_i^cand>0 and Y_i^cand>0
=> no atan2 degeneracy
```

This remains:

- noncyclic (no theta inputs; no populated-instance inputs),
- observer-free (no `K_obs` as primary source),
- and candidate-only (no object-support discharge; `N302` remains in force).

## Inputs reused

1. `T126/F279/N391`
   - provider-object carrier orbit-quotient candidate.
2. `T115/F268/N380`
   - strict phasor reduction pair-map-rule candidate form.
3. `R1`
   - residual-datum target-slot export packet (codomain scaffold).
4. `T132`
   - provider positive-window corridor spec (delta-scaled orbit-depth carrier).
5. `N302`
   - boundary below actual object support remains in force.

## Packet result

`F287` exports one packaged projection candidate:

```text
Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_v1
```

materialized by the persisted artifact instance:

```text
fundamental_action_reconstruction/generated/provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_positive_window_instance.json
```

## Status discipline

This packet does **not** claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual export-map object export,
4. actual theta export,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

It claims only:

1. one explicit provider-carrier projection candidate exists,
2. it is noncyclic and observer-free by explicit contract,
3. it is protected against the `atan2` degeneracy frontier by a typed
   positive-window corridor constraint.

