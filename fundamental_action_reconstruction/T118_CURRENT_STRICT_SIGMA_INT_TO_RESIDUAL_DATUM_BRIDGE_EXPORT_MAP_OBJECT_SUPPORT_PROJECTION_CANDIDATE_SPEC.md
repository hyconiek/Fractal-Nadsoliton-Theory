# T118 Current Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate Spec

Status: `T118_CURRENT_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`N302` isolates a precise missing layer on the residual-datum / `sigma_int_candidate`
route:

```text
missing: actual bridge/export-map object support witness
missing: actual object-to-map support projection
```

`N381/N382` show that the repo can already build, at candidate level:

1. a finite pair-indexed carrier field `E_pair`,
2. a noncyclic reduction `E_pair -> (theta_1^cand, theta_2^cand)`,
3. a sigma-int-driven generator `sigma_int_candidate -> E_pair`.

The next honest move is therefore not a closure relabel.
It is to attempt a **typed projection** from the internal datum side into the
bridge/export-map object-support lane, while staying explicitly below:

- discharge of `T2`,
- discharge of `QW-2191`,
- actual theta export,
- actual bridge/export-map object support,
- ToE closure.

`T118` specifies one admissible **candidate** projection interface intended to
attack the `N302` missing object-support projection layer without cycles.

## Proposed candidate projection object

Export one explicit candidate projection:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_v1
```

with intended type:

```text
Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_v1 :
  (sigma_int_candidate, eps)
    -> residual_datum_bridge_export_map_object_support_projection_candidate_instance
```

## Construction (composition; noncyclic)

Inputs:

1. `sigma_int_candidate ∈ {+1,-1}` (`B4`),
2. `eps ∈ [0,1]` (candidate parameter).

Steps:

1. Generate `E_pair` without theta inputs:

   ```text
   E_pair := G_sigma_int_to_E_pair_generator_candidate_v1(sigma_int_candidate, eps)
   ```

2. Reduce `E_pair` through the strict kernel coupling channel:

   ```text
   (theta_1^cand, theta_2^cand) := M_fractal_light_path_pair_map_rule_candidate_v1(E_pair)
   ```

3. Package a residual-datum target-slot *candidate population record* below
   `R1` (no claim of actual population):

   ```text
   residual_orientation_datum_target_slot_candidate_population_record :=
     (
       theta_1^cand, theta_2^cand,
       u_1^cand = cos(theta_1^cand)c_1 + sin(theta_1^cand)s_1,
       u_2^cand = cos(theta_2^cand)c_2 + sin(theta_2^cand)s_2,
       S_orient^cand = span{u_1^cand,u_2^cand}
     )
   ```

4. Emit the object-support projection instance as a persisted artifact
   containing:
   - inputs `(sigma_int_candidate, eps)`,
   - generated `E_pair`,
   - computed phasor components and `theta_i^cand`,
   - the target-slot candidate population record,
   - explicit `current_absence` / `forbidden_claims`.

Noncyclic contract:

1. no `theta` is used as input,
2. no populated basis-pair instance is used as input.

Observer-free contract:

1. no `K_obs`-indexed selection is used as a primary source.

## What this does and does not accomplish

This projection is admissible only as:

```text
candidate-level object-to-map support projection attempt
```

It must not be promoted into:

1. a strict-core equivalence theorem,
2. an actual bridge/export-map object support witness,
3. an admissible `S_sel_int`,
4. strict-core selector closure,
5. `QW-2191` discharge,
6. ToE closure.

