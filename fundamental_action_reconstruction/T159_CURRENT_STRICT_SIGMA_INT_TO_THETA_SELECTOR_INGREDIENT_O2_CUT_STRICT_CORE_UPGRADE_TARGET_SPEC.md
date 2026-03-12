# T159 Current Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Strict-Core Upgrade Target Spec

Status: `T159_CURRENT_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_STRICT_CORE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

`T157` targeted a **candidate-ingredient** that cuts the `QW-2191` continuous `O(2)` family
in a declared scope.

On the current repo state:

1. `T157` is discharged at the candidate-ingredient level (`F325/N436/P400`),
2. the strict lane also exports explicit strict-provenance value objects for:
   - eps amplitude (`F317/N428`: `eps := 1/2`, premise-based),
   - delta_d step (`F328/N440`: `delta_d := delta_max`, premise-based),
3. `N437` and `N441` prove the computed theta-pair depends on admissible
   corridor-step and amplitude choices (`delta_d` and `eps`), i.e. both remain
   real selector slots on the current strict sigma-int → theta lane,
4. but `QW-2191` strict-core uniqueness remains open and strict-core theta export remains absent.

`T159` names the next missing object sharply:

```text
upgrade from a delta_d/eps-parameterized candidate representative
to a strict-core theta-supply / selector ingredient whose O(2)-cut is canonical
in strict core (no remaining hidden selector slots).
```

This is narrower than ToE closure. It targets only the missing **strict-core upgrade**
needed if one wants to stop relying on “one chosen representative” language for `QW-2191`.

## Target object

If achieved, export one strict-core upgraded selector ingredient object:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_v1
```

with intended meaning:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_v1 :
  sigma_int_strict_derived_v1  ->  (theta_1, theta_2)
```

such that:

1. `(theta_1, theta_2)` are **strict-core** theta values (not merely a candidate record),
2. the induced basis representatives `u_1,u_2` cut the `QW-2191` `O(2)` family to a canonical
   representative (up to the unavoidable residual `Z2` sign convention),
3. no hidden selector slot remains (no silent `delta_d` choice, no silent `eps` choice).

## Acceptance tests (what would count as strict-core upgrade)

An **actual discharge** of this target must at minimum provide:

1. **Typed strict output:** explicit exported values `(theta_1, theta_2)` and induced
   basis vectors on the declared `QW-2190` scaffold (or an equivalent declared strict scaffold).
2. **Selector-slot closure:** the construction must eliminate the currently exposed strict-lane
   selector slots, by one of the following (must be explicit which):
   - **(A) slot elimination:** prove the theta output is invariant under all admissible
     corridor-step choices `delta_d ∈ (0, delta_max]` and under all admissible eps choices
     allowed by the generator contract (no hidden dependence), or
   - **(B) strict derivation:** export strict-derived (not premise-only) provenance chains
     that uniquely derive the needed `delta_d` and `eps` values on a declared strict domain,
     and prove that no additional free selector slot remains.
3. **Noncyclic contract:** no `theta` inputs and no populated basis-pair instance as input
   (respects sandbox `N18`).
4. **Observer-free contract:** no `K_obs`-indexed selection as a primary source.
5. **QW-2191 compatibility:** the theorem-level statement must explicitly acknowledge `QW-2191`:
   - kernel alone does not fix the basis,
   - therefore the upgrade must exhibit a strict-core internal selector source (not a silent convention),
   and must prove that this source canonically cuts the `O(2)` family in the declared scope.
6. **No false pass discipline:** the discharge must not imply:
   - global strict-core selector closure (`S_sel_int`) unless separately proved,
   - global `QW-2191` discharge beyond the declared scope,
   - ToE closure.

## Hard limits

`T159` must not claim:

1. that `F325/N436` already satisfy strict-core upgrade (they are candidate-ingredient only),
2. that premise-based eps/delta_d value objects (`F317/F328`) constitute strict derivation,
3. global selector closure,
4. ToE closure.
