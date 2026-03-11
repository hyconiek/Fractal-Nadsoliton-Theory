# N436 Current First Actual Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Candidate Theorem

Status: `N436_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_CANDIDATE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Package, theorem-level, the strongest honest current statement answering the
`T157` target question:

```text
does the repo export a strict-side **candidate** selector ingredient
sigma_int_strict_derived_v1 -> (theta_1^cand, theta_2^cand)
strong enough to cut the continuous O(2) family from QW-2191 in the declared scope?
```

## Theorem-level conclusion

From `F325`, the current repo exports one strict-side **candidate** selector ingredient object:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1
```

persisted as:

```text
fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json
```

This is a dedicated selector-ingredient packaging of the theta-candidate record already present in the
strict-input positive-window instantiation artifact exported by `F314`
(`generated/sigma_int_strict_derived_v1_to_residual_datum_bridge_export_map_object_support_projection_candidate_positive_window_instance.json`),
augmented with an explicit `O(2)`-cut witness argument and induced `u_i^cand` vectors.

with the following exact meaning:

1. the input sigma-int datum is the exported strict-side source-upgrade value
   `sigma_int_strict_derived_v1 = -1` (`F307/N418`) with explicit premise-based provenance,
2. eps amplitude and the nad12 sign mask are imported only through the exported
   strict-provenance value objects (`F317/N428`, `F324/N435`),
3. the construction uses a positive-window corridor (T119-style) derived from the strict
   working kernel tuple (`QW-2049`) to guarantee `atan2` is well-defined (no degeneracy) **for the
   chosen corridor step `delta_d` recorded in the exported artifact**,
4. the output exports explicit candidate values `(theta_1^cand, theta_2^cand)` and the induced
   candidate basis vectors `u_1^cand,u_2^cand` on the `QW-2190` deterministic mode scaffold,
5. therefore a single declared representative point in each degenerate `O(2)` rotation orbit
   from `QW-2191` is selected in the declared scope **for that chosen `delta_d`**: the continuous
   family is cut only after the extra corridor-step choice (up to the residual sign convention).

No false pass: `T119` admits `delta_d ∈ (0, delta_max]`, and `P403/N437` record explicit theta-pair
dependence on admissible `delta_d` choices. Therefore `N436` remains a *candidate selector ingredient
existence* statement, not a strict-core uniqueness claim.

## What N436 proves

`N436` proves only:

1. existence of one exported strict-side candidate selector ingredient object
   `sigma_int_strict_derived_v1 -> (theta_1^cand, theta_2^cand)` (`F325`),
2. explicit nondegeneracy certification via the positive-window corridor (no `atan2(0,0)`),
3. an explicit `O(2)`-cut witness statement in the declared `QW-2190/QW-2191` scope.

## What N436 does not prove

`N436` does not prove:

1. kernel-alone discharge of `QW-2191`,
2. actual strict-core `theta_1`, `theta_2` export,
3. admissible `S_sel_int` or strict-core selector closure,
4. discharge of post-`T148` object-support targets (`N302/N395/T130`),
5. ToE closure.

## Status discipline

All premise-based inputs remain explicitly marked as premise-based in the exported artifact.
No “hybrid FR reuse” is used as strict evidence.
