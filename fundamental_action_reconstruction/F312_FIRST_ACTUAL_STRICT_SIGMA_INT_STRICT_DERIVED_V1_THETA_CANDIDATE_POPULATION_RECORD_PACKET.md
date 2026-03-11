# F312 First Actual Strict Sigma-Int Strict-Derived v1 Theta-Candidate Population Record Packet

Status: `F312_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_STRICT_DERIVED_V1_THETA_CANDIDATE_POPULATION_RECORD_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `P390`, the repo has already discharged:

1. `T149` (strict sigma-int source upgrade) via `F307/N418`, and
2. `T148` (actual strict residual-datum sigma-int export-map object) via
   `F311/N422` (residual `Z2` population only; no theta export).

The next honest missing layer remains explicit:

```text
strict-core theta-source export and/or residual target-slot population (absent)
```

This packet does **not** pretend to export actual `theta_1`, `theta_2`.

It executes a narrower, audit-safe move:

```text
instantiate the already-exported strict-side projection-candidate pipeline
at the now-exported strict sigma-int datum sigma_int_strict_derived_v1 = -1,
and persist one candidate theta-pair population record for the R1 target-slot scaffold
without theta inputs and without populated-instance inputs.
```

This keeps the strict discipline:

- noncyclic (no theta input; no populated-instance input),
- observer-free (no `K_obs` primary source),
- no implied selector closure; `QW-2191` remains open.

## Inputs reused (strict-admissible)

1. `F307/N418`
   - strict sigma-int source upgrade:
     `sigma_int_strict_derived_v1 ∈ {+1,-1}` with value `-1` (premise-based).
2. `T115/F268`
   - strict fractal-light pair-map-rule **candidate** form
     `M_fractal_light_path_pair_map_rule_candidate_v1` (atan2-reduction).
3. `T117/F270`
   - strict sigma-int driven `E_pair` generator **candidate**
     `G_sigma_int_to_E_pair_generator_candidate_v1`.
4. `T118/F271`
   - strict sigma-int → residual-datum projection **candidate**
     `Pi_sigma_int_to_residual_datum_bridge_export_map_object_support_projection_candidate_v1`.
5. `R1`
   - residual orientation datum target-slot export scaffold
     `residual_orientation_datum_target_slot` (unpopulated in strict core).

## Packet result (strict-input instantiation; still candidate-level output)

`F312` exports one strict-input instantiation artifact for the already exported
projection-candidate object (still **candidate-only** output):

```text
fundamental_action_reconstruction/generated/
  sigma_int_strict_derived_v1_to_residual_datum_bridge_export_map_object_support_projection_candidate_instance.json
```

This artifact:

1. uses `sigma_int_strict_derived_v1 = -1` as the sigma-int input value,
2. uses the exported eps value object
   `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2` (`F317/N428`),
3. generates a finite `E_pair` carrier by the `T117` candidate generator,
4. reduces `E_pair` to `(theta_1^cand, theta_2^cand)` via the `T115` candidate
   phasor rule,
5. records a candidate population record for the `R1` target-slot scaffold.

### Recorded candidate outputs (for this one instance only)

```text
theta_1^cand ≈ 0.4625142242896771
theta_2^cand ≈ 0.3554808343611710
```

This is only an instantiated candidate record. It is not an actual theta-source export.

## Status discipline

This packet does **not** claim:

1. actual strict-core `theta_1`, `theta_2`,
2. actual populated basis-pair instance,
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.
