# T162 Current Strict Sigma-Int → Theta Slot-Free Construction-Class Target Spec

Status: `T162_CURRENT_STRICT_SIGMA_INT_TO_THETA_SLOT_FREE_CONSTRUCTION_CLASS_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

On the current strict sigma-int → theta lane, the exported candidate construction class contains two
exposed selector slots:

1. `eps ∈ [0,1]` (generator amplitude; `T117`), and
2. `delta_d ∈ (0, delta_max]` (positive-window corridor step; `T119`).

The current repo also exports theorem-level sensitivity results:

- theta depends on admissible `eps` (`P407/N441`),
- theta depends on admissible `delta_d` (`P403/N437`),

so the invariance-based slot-elimination route is closed negatively for the current class (`N443`).

The strict-derived slot-selection route is named (`T160`, `T161`) but not discharged (`P410/P411`, `N444/N445`),
and the common proposed “Final Stroke” derivation slogans (charge parity balance / maximum information packing)
are closed negatively as strict-derived sources (`N446`, `N447`, packaged as `N448`).

Therefore, on the current repo state, the only remaining honest strict-core upgrade route for the sigma-int
selector ingredient is:

```text
change the construction class so the eps / delta_d slots do not exist.
```

`T162` names that missing ingredient class sharply, as a future-only target with explicit acceptance tests.

## Scope

`T162` is scoped only to the strict sigma-int → theta selector-ingredient frontier on `QW-2191`.

It does **not** decide:

1. discharge of post-`T148` object-support targets (`N302/N395/T130`),
2. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge beyond the declared scope,
3. ToE closure.

## Target object

If achieved, export one strict-core selector ingredient object in a **slot-free construction class**:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1
```

with intended meaning:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1 :
  sigma_int_strict_derived_v1  ->  (theta_1, theta_2)
```

such that:

1. `(theta_1, theta_2)` are strict-core theta values (not merely a candidate record),
2. the induced basis representatives canonically cut the `QW-2191` `O(2)` family (up to residual `Z2` sign),
3. the construction contains **no free selector slots** of the `eps` / `delta_d` kind:
   - no generator amplitude parameter family,
   - no corridor-step discretization family.

Equivalently: all inputs are strict-core objects or strict-derived values already exported on the declared
strict domain, and the construction does not require choosing a representative from a parameter family.

## Acceptance tests (what would count as discharge)

An **actual discharge** of this target must at minimum provide:

1. **Typed strict output:** explicit exported values `(theta_1, theta_2)` on a declared strict scaffold.
2. **Slot-free construction:** the construction must not quantify over a family parameter:
   - no `eps ∈ [0,1]` slot,
   - no `delta_d ∈ (0, delta_max]` slot,
   - no replacement hidden slot of the same kind.
3. **Noncyclic contract:** no theta inputs and no populated basis-pair instance as input (`N18`).
4. **Observer-free contract:** no `K_obs`-indexed selection as a primary source.
5. **QW-2191 compatibility:** the discharge must explicitly identify the strict-core internal selector source
   that canonically cuts the `O(2)` family (not a silent convention), and state its scope.
6. **No false pass discipline:** the discharge must not imply:
   - global strict-core selector closure (`S_sel_int`),
   - global `QW-2191` discharge beyond the declared scope,
   - ToE closure,
   unless separately proved.

## Hard limits

`T162` must not claim:

1. that the current exported sigma-int → theta class (`T117/T119` family) is already slot-free,
2. that premise-based slot point choices (`F317/F328`) constitute slot elimination,
3. strict-core theta export on the current repo state,
4. ToE closure.

