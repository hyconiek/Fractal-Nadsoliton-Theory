# T161 Current Strict Sigma-Int Positive-Window Delta_d Step Strict-Derived Selection Target Spec

Status: `T161_CURRENT_STRICT_SIGMA_INT_POSITIVE_WINDOW_DELTA_D_STEP_STRICT_DERIVED_SELECTION_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After the positive-window corridor spec (`T119`) the strict sigma-int → theta
candidate lane uses a typed corridor step parameter:

```text
delta_d ∈ (0, delta_max]   where   delta_max := d_local/11.
```

After the delta_d provenance target (`T158`) and its discharge (`F328/N440`),
the strict sigma-int lane exports one dedicated delta_d value object:

```text
delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max
```

with explicit provenance classification:

```text
strict_source_upgraded (explicit premise; not strict-derived).
```

However, `P403/N437` prove that the computed theta-pair depends on admissible
delta_d choices inside the corridor. Therefore delta_d remains a real selector
slot and cannot be silently upgraded to strict-core by citing only one
premise-based corridor-saturation choice.

`T159` names a strict-core upgrade target for an `O(2)`-cut selector ingredient.
After `N443`, the invariance-based slot-elimination route is blocked on the
current exports. The only remaining honest route for slot removal in strict
core is:

```text
export a genuinely strict-derived (not premise-only) delta_d selection law/value object,
or change the construction class so the slot does not exist.
```

`T161` names the missing **strict-derived delta_d selection** ingredient sharply,
as a future-only target with explicit acceptance tests.

## Extension-lane continuation note (explicitly non-strict)

On the current repo state, the strict-core delta_d selection target remains open:
no strict-derived delta_d law/value object is exported (`N447`, `P408`).

If one insists on proceeding *today* with a single reproducible sigma-int → theta
representative without false pass, the repo explicitly separates that move into
`strict_extension_only` scope:

1. `AX21` freezes `delta_d := delta_max` as a declared corridor-saturation premise
   (premise-based; not strict-derived).
2. `AX22` packages a publication-ready strict-extension summary of that lane.

This does **not** discharge `T161`. It only records the explicit extension-scope
premise needed for a reproducible representative while keeping strict-core delta_d
selection open.

## Scope

`T161` is scoped only to the **delta_d step selection** used by the strict
positive-window sigma-int → theta candidate constructions.

It does **not** decide:

1. strict-core `theta_1`, `theta_2` export,
2. discharge of post-`T148` object-support targets (`N302/N395/T130`),
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

## Target object

If the repo cannot yet export a strict-derived delta_d selection ingredient,
export one explicit future-only target object:

```text
Delta_sigma_int_positive_window_delta_d_step_strict_derived_selection_target_v1
```

with intended meaning:

```text
export one dedicated delta_d value object with STRICT_DERIVED provenance
(not premise-only), together with an explicit strict derivation/selection chain
on a declared strict domain, so the strict sigma-int → theta construction no
longer depends on a free delta_d selector slot in strict core
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_sigma_int_positive_window_delta_d_step_strict_derived_selection_target_v1`
must at minimum provide:

1. **Dedicated delta_d value object:** an exported value object
   `delta_d_sigma_int_positive_window_step_strict_derived_v1` with contract:
   - `0 < delta_d_sigma_int_positive_window_step_strict_derived_v1 <= delta_max`,
   - where `delta_max := d_local/11` is the corridor bound from `T119`.
2. **Strict-derived provenance (no premise-only):** the object must be explicitly
   classified as `strict_derived` with an explicit derivation/selection chain on a
   declared strict domain. It must not be only a `strict_source_upgraded` premise
   (e.g. `delta_d := delta_max`) and must not be left as an implicit convention.
3. **Slot-closure discipline:** the derivation/selection chain must not smuggle
   in any free selector slot. In particular:
   - it must not take `theta_{1,2}` outputs or any populated basis-pair instance as input (`N18`),
   - if it depends on `eps`, then `eps` must itself be supplied by a strict-derived
     value object (not a premise-only amplitude convention), and the dependence must be explicit.
4. **Observer-free contract:** no `K_obs`-indexed selection as a primary source.
5. **Compatibility with `T159`:** the discharge must state how this strict-derived delta_d object
   removes the delta_d selector slot required by the `T159` strict-core upgrade acceptance test (2B).
6. **No false pass discipline:** the discharge must not imply:
   - strict-core theta export,
   - discharge of `N302/N395/T130`,
   - admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
   unless separately proved.

## Relation to existing exports (current-state discipline)

On the current repo state:

1. delta_d is exported only as strict provenance (premise-based):
   `delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max` (`F328/N440`),
   which does **not** discharge `T161`,
2. `delta_max` is currently only a corridor bound (`T119`), and no strict “maximum information packing”
   objective exists on this lane selecting `delta_d = delta_max`, packaged as `N447`,
3. no strict exported theorem currently derives the corridor-saturation choice from
   any typed strict “information saturation” or uniqueness objective (`P408` audit).

## Hard limits

`T161` must not claim:

1. that delta_d is already strict-derived on the current strict sigma-int lane,
2. that premise-based delta_d provenance (`F328/N440`) constitutes strict derivation,
3. strict-core theta export or strict-core selector closure,
4. ToE closure.
