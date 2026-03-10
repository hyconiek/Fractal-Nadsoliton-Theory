# N400 Current First Actual Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate (Positive-Window Low-Delta) Theorem

Status: `N400_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_LOW_DELTA_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Export one additional corridor-protected sigma-int projection candidate whose
typed `d`-domain restriction guarantees `atan2` nondegeneracy and whose chosen
instance keeps `theta_i^{cand}` closer to the origin-phase scale (useful for
cross-lane convergence-side coherence comparison), while preserving the `N302`
boundary.

## Theorem-level conclusion

From `T133/P374/F288`, the current repo exports one actual packaged low-delta
positive-window projection candidate:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_v1
```

materialized as the persisted artifact instance:

```text
fundamental_action_reconstruction/generated/sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_low_delta_instance.json
```

with the following exact meaning:

1. the projection is noncyclic (no theta inputs; no populated-instance inputs),
2. the projection is observer-free (no `K_obs` primary source),
3. the projection is pair-indexed on `[pair1,pair2]`,
4. the projection remains candidate-level only,
5. the projection is corridor-protected:
   - its carrier `E_pair` satisfies `d_{i,k} ∈ [0,d^{local}]`,
   - hence each term has `cos(Theta(d_{i,k}))>0` and `sin(Theta(d_{i,k}))>0`,
6. therefore for each pair slot `i`:
   - `X_i^{cand}>0` and `Y_i^{cand}>0`,
   - so `atan2(Y_i^{cand},X_i^{cand})` is well-defined.

## Operational evaluation (scope-limited)

For the persisted instance choice:

```text
sigma_int_candidate = +1
eps = 1/2
delta_d = 0.05
```

the projection produces:

```text
pair1: (X_cand,Y_cand,theta_cand) ≈ (0.8819272510, 0.1864891620, 0.2083866709)
pair2: (X_cand,Y_cand,theta_cand) ≈ (0.8700903452, 0.1880349424, 0.2128365652)
```

This is an operational check only and does not imply uniqueness.

## What N400 does not prove

`N400` does not prove:

1. discharge of `T2`,
2. satisfaction of the strict-core bridge/export-map object target `N301`,
3. discharge of `N302`,
4. actual residual bridge/export-map object support,
5. actual bridge/export map export,
6. actual `theta_1`, `theta_2`,
7. actual pair population,
8. admissible `S_sel_int`,
9. strict-core selector closure,
10. `QW-2191` discharge,
11. ToE closure.

