# F325 First Actual Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Candidate Packet

Status: `F325_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_CANDIDATE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After the strict sigma-int provenance + gauge-quotient safety package (`F307/N418`, `F308/N419`)
and after the strict-provenance upgrades for:

- eps amplitude (`F317/N428`: `eps=1/2`),
- nad12 sign mask (`F324/N435`: `b_{i,k}` sourced from `chi_FR_strict_v1`),

the strict lane still lacked an explicit internal **theta-supply / selector ingredient**
capable of cutting the continuous `O(2)` family from `QW-2191` without importing the
axiom-augmented selector (`QW-2192/2193`) as strict evidence.

This packet executes the narrow `T157` move at the **candidate ingredient** level:

```text
export one strict-side candidate selector ingredient:
  sigma_int_strict_derived_v1 -> (theta_1^cand, theta_2^cand),
and record the induced O(2)-cut witness on the QW-2190 mode scaffold.
```

## Inputs reused (strict-admissible; explicit provenance)

1. `F307/N418`
   - strict-side sigma-int source upgrade: `sigma_int_strict_derived_v1 = -1` (premise-based; no hybrid reuse).
2. `F317/N428`
   - eps amplitude value object: `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2` (premise-based).
3. `F324/N435`
   - strict-provenance nad12 sign-mask value object: `b_sigma_int_E_pair_sign_mask_strict_provenance_v1`.
4. `T119`
   - positive-window corridor formula guaranteeing `atan2` nondegeneracy.
5. `F314`
   - existing strict-input positive-window instantiation artifact already recording the same numeric
     theta-candidate pair inside the sigma-int → residual-datum projection pipeline:
     `generated/sigma_int_strict_derived_v1_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance.json`.
6. `material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2049_spectral_micro_stagec_intersection_gate.json`
   - strict working kernel tuple `(omega,phi,beta,eta)` used by the corridor.
7. `QW-2190` / `QW-2191`
   - deterministic real Fourier mode scaffold and its continuous `O(2)` nonuniqueness obstruction.

## Packet product

`F325` exports one persisted strict-side **candidate** selector ingredient artifact:

```text
fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json
```

with exact meaning:

1. define the positive-window corridor step `delta_d` from the strict working kernel tuple (as in `T119`),
2. generate a finite nad12 carrier `E_pair` from:
   - `sigma_int_strict_derived_v1`,
   - `eps_sigma_int_E_pair_amplitude_strict_provenance_v1`,
   - `b_sigma_int_E_pair_sign_mask_strict_provenance_v1`,
3. compute phasor sums `(X_i,Y_i)` and define:
   - `theta_i^cand := atan2(Y_i, X_i)`,
4. record `(theta_1^cand, theta_2^cand)` and the induced basis vectors `u_1^cand,u_2^cand` on the `QW-2190` mode scaffold,
5. record the explicit `O(2)`-cut argument in the declared scope (one selected representative point).

### Exported values (current strict tuple; sigma-int = -1)

From the persisted artifact:

```text
theta_1^cand ≈ 0.3627333053541785
theta_2^cand ≈ 0.33287066305007096
```

These numeric values match (up to rounding) the `theta_cand` entries already persisted by `F314` in:

```text
fundamental_action_reconstruction/generated/
  sigma_int_strict_derived_v1_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance.json
```

The novelty of `F325` is therefore **semantic packaging**, not new computation:
we re-export the theta pair as a dedicated strict-side selector-ingredient object and attach an explicit
`O(2)`-cut witness argument (with induced `u_1^cand,u_2^cand`).

## Status discipline

This packet remains below actual theta export:
it exports only a **candidate ingredient** object with explicit premise-based provenance recorded in the artifact.

This packet does **not** claim:

1. kernel-alone discharge of `QW-2191`,
2. actual strict-core `theta_1`, `theta_2` export,
3. admissible `S_sel_int` or strict-core selector closure,
4. ToE closure.
