# T133 Current Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate (Positive-Window Low-Delta) Spec

Status: `T133_CURRENT_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_LOW_DELTA_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T119/F272/N384` exports one corridor-protected sigma-int projection candidate:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_v1
```

materialized by one persisted instance choice. That instance is sufficient to
establish:

```text
atan2 nondegeneracy (X_i^cand>0, Y_i^cand>0) on a typed corridor.
```

However, strict closure convergence work (see `T114`) benefits from having one
*calibrated* corridor-protected instance whose output phases remain close to
the strict kernel origin-phase scale (i.e. close to `phi`) so that it can be
coherently compared against other noncyclic projection lanes (e.g. the
provider-object carrier projection lane).

`T133` introduces one additional, still corridor-protected, **low-delta**
positive-window variant as a separate named candidate object.

This remains:

- noncyclic (no theta inputs; no populated-instance inputs),
- observer-free (no `K_obs` as primary source),
- and candidate-only (no discharge of `N302` and no bridge/export-map export).

## Positive-window corridor

Reuse the same corridor definition from `T119`:

```math
\delta_d \in (0, d^{local}/11].
```

`T133` fixes a *low-delta calibration choice*:

```math
\delta_d := 0.05.
```

This satisfies the corridor bound for the current strict working tuple
(`0.05 <= d^{local}/11`).

## Candidate object (low-delta variant)

Export one explicit low-delta positive-window projection candidate:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1
```

with intended type:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1 :
  (sigma_int_candidate, eps)
    -> residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_instance
```

Construction:

1. reuse the delta-scaled nad12 carrier generator from `T119`,
2. fix `delta_d := 0.05`,
3. reduce via `M_fractal_light_path_pair_map_rule_candidate_v1` (`T115`),
4. emit a residual-datum target-slot **candidate** population record below `R1`.

## Hard limits

`T133` must not claim:

1. discharge of `T2`,
2. discharge of `N302`,
3. satisfaction of the `N301` bridge/export-map object target,
4. actual bridge/export-map object support,
5. actual bridge/export map export,
6. admissible `S_sel_int`,
7. strict-core selector closure or `QW-2191` discharge,
8. ToE closure.

