# T159 Current Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Strict-Core Upgrade Target Spec

Status: `T159_CURRENT_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_STRICT_CORE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

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
4. update (`2026-03-15`): a slot-free strict-core sigma-int → theta-pair source is exported, discharging `T162` and
   satisfying the `T159` strict-core upgrade target in the declared `R1` scope (`F451`, packaged by `N489`), while
   keeping global `QW-2191` kernel-alone uniqueness open.

`T159` names the strict-core upgrade object sharply:

```text
upgrade from a delta_d/eps-parameterized candidate representative
to a strict-core theta-supply / selector ingredient whose O(2)-cut is canonical
in strict core (no remaining hidden selector slots).
```

Update (`2026-03-15`): this target is now satisfied in the declared `R1` scope by the slot-free construction-class
change route `T162` (`F451/N489`). The slot-elimination and strict-derived slot-selection routes remain open only for
the older parameterized corridor class.

## Current state note (no false pass)

On the current repo state, the slot-elimination route (acceptance test 2A) is already blocked:

1. `N437` proves theta-pair dependence on admissible `delta_d`,
2. `N441` proves theta-pair dependence on admissible `eps`,
3. therefore the current exported sigma-int → theta candidate pipeline cannot satisfy slot-elimination
   by invariance (packaged as `N443`).

So for the older delta_d/eps-parameterized corridor class, the slot-elimination upgrade route is blocked; the only
remaining strict continuation would be strict-derived slot selection (acceptance test 2B), or a construction-class
change (2C).

Update (`2026-03-15`): the construction-class change route is now discharged (`T162`) and satisfies `T159` in the
declared `R1` scope (`F451`, packaged by `N489`).

## Extension-lane continuation note (explicitly non-strict)

Even though the strict-core upgrade target `T159` is now satisfied via the slot-free `T162` route, the repo still
records the following continuation for the older parameterized sigma-int → theta class: if one insists on a single
reproducible representative of that parameterized corridor class **today**, one must work in an explicit separated
scope and declare the selector choices as premises.

The repo now records exactly that continuation:

1. `AX21` freezes both exposed selector slots in `strict_extension_only` scope:
   - `eps := 1/2`,
   - `delta_d := delta_max := d_local/11`.
2. `AX22` packages a publication-ready strict-extension summary of this lane.

This does **not** discharge `T160/T161` and does not upgrade the old parameterized corridor class into strict core.
It only makes the extension-lane representative explicit and reproducible while keeping strict-core claims unchanged.
`T159` is satisfied separately via the slot-free `T162` route (`F451/N489`).

On the current repo state, the strict-derived slot-selection ingredients required by (2B) are
not exported yet and remain future-only targets for the old parameterized class:

1. eps strict-derived selection target: `T160` (status packaged by `P410/N444`),
2. delta_d strict-derived selection target: `T161` (status packaged by `P411/N445`).

## Target object

Target object (strict-core theta supply):

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1
```

Exported meaning:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1 :
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
     The strict-derived slot-selection targets are named as `T160` (eps) and `T161` (delta_d).
   - **(C) construction-class change:** export a genuinely new strict sigma-int → theta construction
     class in which the `eps` / `delta_d` slot families do not exist at all, and prove that no
     replacement hidden selector slot has been introduced.
     The slot-free construction-class target is named as `T162`.
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
