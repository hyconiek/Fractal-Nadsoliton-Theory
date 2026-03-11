# F314 First Actual Strict Sigma-Int Strict-Derived v1 Theta-Candidate Population Record (Positive-Window) Packet

Status: `F314_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_STRICT_DERIVED_V1_THETA_CANDIDATE_POPULATION_RECORD_POSITIVE_WINDOW_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `F312/N423`, the repo already exports a strict-input instantiation
artifact of the sigma-int → residual-datum **projection-candidate** pipeline,
evaluated at:

```text
sigma_int_strict_derived_v1 = -1
```

However, that instantiation does not enforce the positive-window corridor
introduced in `T119/F272` for avoiding the explicit `atan2` degeneracy frontier
named in `T115`.

This packet performs the next honest move:

```text
instantiate the already-exported positive-window projection-candidate pipeline
at the strict sigma-int datum sigma_int_strict_derived_v1 = -1,
derive a maximal admissible delta_d from the strict kernel tuple,
and persist one positive-window candidate theta-pair population record
for the R1 target-slot scaffold.
```

Hard limits are preserved:

- noncyclic (no theta input; no populated-instance input),
- observer-free (no `K_obs` primary source),
- candidate-only output (no actual theta export; no selector closure; `QW-2191` stays open).

## Inputs reused (strict-admissible)

1. `F307/N418`
   - strict sigma-int source upgrade:
     `sigma_int_strict_derived_v1 = -1` (premise-based provenance, no hybrid reuse).
2. `T115/F268`
   - strict fractal-light pair-map-rule **candidate** form
     `M_fractal_light_path_pair_map_rule_candidate_v1` (atan2-reduction).
3. `T119/F272`
   - positive-window corridor projection-candidate form
     `Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_v1`.
4. `R1`
   - residual orientation datum target-slot export scaffold
     `residual_orientation_datum_target_slot` (unpopulated in strict core).

## Derived positive-window step parameter (delta-scaled nad12 corridor)

Using the strict working kernel tuple:

```text
omega = 0.18575
phi   = 0.16250
beta  = 1.0
eta   = 1.8
```

define (as in `T119`):

```text
delta_barrier := (pi/2) - |phi|
eps_local     := delta_barrier / 2
d_local       := eps_local / omega
delta_max     := d_local / 11
```

and instantiate the positive-window corridor with the maximal admissible step:

```text
delta_d := delta_max.
```

This guarantees:

```text
all paths satisfy d in [0, d_local]
=> cos(Theta(d))>0 and sin(Theta(d))>0
=> X_i^cand>0 and Y_i^cand>0
=> atan2 well-defined (no degeneracy).
```

## Packet result (strict-input instantiation; still candidate-level output)

`F314` exports one strict-input instantiation artifact for the already exported
positive-window projection-candidate object (still **candidate-only** output):

```text
fundamental_action_reconstruction/generated/
  sigma_int_strict_derived_v1_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance.json
```

This artifact:

1. uses `sigma_int_strict_derived_v1 = -1` as the sigma-int input value,
2. fixes one explicit parameter choice `eps = 1/2` (candidate-only; not strict-derived),
3. derives `delta_d` from the strict kernel tuple via the `T119` corridor,
4. generates a finite `E_pair` carrier by the positive-window generator form,
5. reduces `E_pair` to `(theta_1^cand, theta_2^cand)` via the `T115` candidate
   phasor rule,
6. records a candidate population record for the `R1` target-slot scaffold.

### Recorded candidate outputs (for this one instance only)

```text
theta_1^cand ≈ 0.36273330535417875
theta_2^cand ≈ 0.33287066305007130
```

This is only an instantiated candidate record. It is not an actual theta-source export.

## Status discipline

This packet does **not** claim:

1. actual strict-core `theta_1`, `theta_2`,
2. actual populated basis-pair instance,
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

