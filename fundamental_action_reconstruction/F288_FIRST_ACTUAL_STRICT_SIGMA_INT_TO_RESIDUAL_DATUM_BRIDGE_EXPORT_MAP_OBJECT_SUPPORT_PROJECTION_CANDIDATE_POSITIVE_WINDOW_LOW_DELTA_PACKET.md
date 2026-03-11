# F288 First Actual Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate (Positive-Window Low-Delta) Packet

Status: `F288_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_LOW_DELTA_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`F272/N384` already export one corridor-protected sigma-int projection
candidate. The strict closure convergence-side evaluation (`T114`) benefits
from also having one additional corridor-protected instance whose output stays
closer to the origin-phase scale (close to `phi`) for coherent cross-lane
comparison.

`F288` packages such a variant, still strictly below:

1. discharge of `N302`,
2. actual bridge/export-map object support above the map object (`N395` remains
   future-only),
3. strict-core theta export / pair population (`N18` loop not broken),
4. admissible `S_sel_int` / selector closure.

## Inputs reused

1. `B4`
   - `sigma_int_candidate ∈ {+1,-1}`.
2. `T115/F268/N380`
   - `M_fractal_light_path_pair_map_rule_candidate_v1`.
3. `R1`
   - residual-datum target-slot export packet (codomain scaffold).
4. `T119`
   - positive-window corridor spec (typed delta-scaled nad12 carrier).
5. `T133`
   - low-delta calibration specialization (still within the corridor).
6. `N302`
   - boundary below actual object support remains in force.

## Packet result

`F288` exports one packaged low-delta positive-window projection candidate:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1
```

materialized by the persisted artifact instance:

```text
fundamental_action_reconstruction/generated/sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_instance.json
```

## Status discipline

This packet does **not** claim:

1. discharge of `T2`,
2. discharge of `N302`,
3. actual bridge/export-map object support above the map object (`N395`),
4. any new bridge/export-map object export beyond `F311/N422`,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

It claims only:

1. one explicit low-delta positive-window candidate projection artifact exists,
2. it is noncyclic and observer-free by explicit contract,
3. it is protected against the `atan2` degeneracy frontier by the same typed
   corridor constraint as `T119`.
