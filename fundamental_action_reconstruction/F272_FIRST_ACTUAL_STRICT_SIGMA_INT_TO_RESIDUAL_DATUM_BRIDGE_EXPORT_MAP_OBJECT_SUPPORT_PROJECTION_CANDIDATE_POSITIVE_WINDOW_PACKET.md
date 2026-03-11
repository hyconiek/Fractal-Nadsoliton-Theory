# F272 First Actual Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate (Positive-Window) Packet

Status: `F272_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T118/F271/N383` exports one explicit **candidate** object-to-map support
projection artifact, but it remains exposed to the explicit degeneracy frontier
named in `T115`:

```text
(X_i^cand,Y_i^cand)=(0,0) => atan2 undefined
```

`F272` packages a strictly weaker but strictly safer variant:

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
- and candidate-only (no bridge/export-map discharge).

## Inputs reused

1. `B4`
   - `sigma_int_candidate ∈ {+1,-1}` (candidate object; hybrid-supported).
2. `F307/N418`
   - `sigma_int_strict_derived_v1 ∈ {+1,-1}` (strict-side source upgrade; explicit premise provenance).
3. `T115/F268/N380`
   - `M_fractal_light_path_pair_map_rule_candidate_v1`.
4. `R1`
   - residual-datum target-slot export packet (codomain scaffold).
5. `T119`
   - positive-window corridor spec for a delta-scaled nad12 carrier.

## Packet result

`F272` exports one packaged projection candidate:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_v1
```

materialized by the persisted artifact instance:

```text
fundamental_action_reconstruction/generated/sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance.json
```

## Status discipline

This packet does **not** claim:

1. discharge of `T2`,
2. discharge of `N302`,
3. satisfaction of the `N301` bridge/export-map object target,
4. actual bridge/export-map object support,
5. actual `theta_1`, `theta_2`,
6. actual pair population,
7. admissible `S_sel_int`,
8. strict-core selector closure or `QW-2191` discharge,
9. ToE closure.

It claims only:

1. one explicit candidate projection artifact exists,
2. it is noncyclic and observer-free by explicit contract,
3. it is protected against the `atan2` degeneracy frontier by a typed
   positive-window corridor constraint.

No identification theorem between `sigma_int_candidate` (`B4`) and
`sigma_int_strict_derived_v1` (`F307/N418`) is used or implied here; both are
admissible instantiations of the abstract `Z2` input `sigma_int_input ∈ {+1,-1}`.
