# T157 Current Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Target Spec

Status: `T157_CURRENT_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After the strict sigma-int source upgrade (`F307/N418`), theorem-level gauge-quotient safety (`F308/N419`),
and the strict sigma-int → residual target-slot export-map object (`F311/N422`),
the post-`T148` strict bottleneck remains explicit (`N302/N395`):

```text
strict-side theta-supply / selector ingredient is still missing,
and QW-2191 O(2) nonuniqueness remains open in strict core.
```

`T157` targets one narrow missing object required to *cut* the continuous `O(2)` family from `QW-2191`
without importing the axiom-augmented selector family (`QW-2192/2193`) as strict evidence:

```text
an explicit strict-side **candidate** (observer-free, noncyclic) map:
  sigma_int_strict_derived_v1  ->  (theta_1^cand, theta_2^cand)
which selects one point in the O(2) family (basis choice) on the QW-2190 mode scaffold.
```

This is the “strict sigma_int → theta selector map” item from the missing-object list of `P3`
(`strict_core_sigma_int_to_theta_map_or_internal_derivation_of_Jab_selector_family`),
specialized to the strict sigma-int lane.

**Note (non-novelty discipline):** `F314` already records one concrete `(theta_1^cand,theta_2^cand)` pair
inside the strict-input positive-window instantiation artifact of the sigma-int → residual-datum
projection pipeline. The missing item was not the existence of a numeric record, but the lack of a
*dedicated selector-ingredient object* with an explicit `O(2)`-cut witness argument.

## Intended construction class (strict-side, no hybrid reuse)

`T157` admits only a construction that is:

1. **Noncyclic**: no `theta` inputs and no populated-instance inputs (respects sandbox `N18`).
2. **Observer-free**: no `K_obs` primary indexing/choice.
3. **Strict-provenance explicit**:
   - sigma-int strict-side source upgrade (`N418`) is the only sigma-int input,
   - eps amplitude uses the exported strict-provenance value object (`N428`: `1/2`),
   - the nad12 sign-mask uses the exported strict-provenance value object (`N435`),
   - the positive-window corridor uses the strict working tuple and the corridor formula (`T119`),
   - no “hybrid FR reuse” is used as evidence.
4. **O(2)-cut certified**: the output must be strong enough to select a unique representative in each
   degenerate mode-pair rotation orbit from `QW-2191` (at least up to the unavoidable residual `Z2` sign).

## Target object

If the construction is achieved, export one strict-side **candidate** selector ingredient object:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1
```

with intended meaning:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1 :
  sigma_int_strict_derived_v1  ->  (theta_1^cand, theta_2^cand)

where (theta_1^cand,theta_2^cand) selects one basis representative
in the QW-2191 O(2) assignment family on the QW-2190 mode scaffold.
```

The object may additionally export the induced basis vectors:

```text
u_1^cand := cos(theta_1^cand)c_1 + sin(theta_1^cand)s_1
u_2^cand := cos(theta_2^cand)c_2 + sin(theta_2^cand)s_2
```

as explicit witness that the continuous `O(2)` orbit is cut to a single declared point.

## Acceptance tests

`T157` is discharged only at the **candidate ingredient** level if all of the following are satisfied:

1. **Typed output:** explicit exported values `(theta_1^cand, theta_2^cand)` exist in a declared strict domain.
2. **Well-definedness:** `atan2` degeneracy is ruled out by an exported positive-window corridor
   (e.g. `T119`-style corridor) so that both components are computable.
3. **Strict provenance:** the construction cites `N418/N419/N428/N435` explicitly and keeps all
   premise-based inputs marked as such (no silent conventions).
4. **O(2)-cut proof:** the export includes a proof/argument that the constructed datum fixes a unique
   representative in the `O(2)` family from `QW-2191` (for the same mode pairs used by `QW-2190`).
5. **No false pass:** no claim of:
   - strict-core selector closure,
   - `QW-2191` axiom-free discharge,
   - ToE closure,
   is smuggled in by wording.

## Hard limits

`T157` must not claim:

1. that kernel alone discharges `QW-2191` (it does not),
2. that the export is an axiom-free global uniqueness theorem,
3. admissible `S_sel_int` or strict-core selector closure,
4. ToE closure.

`T157` targets only a **new strict-side selector ingredient** that cuts the `O(2)` family in a declared scope.
