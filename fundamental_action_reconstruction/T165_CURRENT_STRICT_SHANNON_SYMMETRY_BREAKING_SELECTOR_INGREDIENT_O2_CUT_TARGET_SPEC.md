# T165 Current Strict Shannon Symmetry‑Breaking Selector Ingredient (O(2)‑Cut) Target Spec

Status: `T165_CURRENT_STRICT_SHANNON_SYMMETRY_BREAKING_SELECTOR_INGREDIENT_O2_CUT_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `QW-2191`, the strict program explicitly exports a continuous `O(2)` non‑uniqueness family for the
kernel‑mode assignment problem: kernel‑alone does not canonically pick a representative.

After `N420` and `N418`, the strict program also exports two strict‑side *source* value objects that are
plausible ingredients for an **internal** symmetry‑breaking mechanism:

1. `alpha_geo_strict_derived_v1 := 4 ln 2` (strict Shannon equipartition witness; `F309/N420`),
2. `sigma_int_strict_derived_v1 := -1` (strict sigma‑int source upgrade on a declared domain; `F306..F308/N418..N419`).

`S2` records the current strategic intention:

```text
selector-axiom discharge via a Shannon symmetry-breaking premise
```

But the repo still exports **no** strict‑core internal selector source (`N124`), and the existing
Shannon‑weighted objects on the nad12‑sigma route are packaged only as *refinement candidates*
(`F319..F322`, `N430..N433`).

`T165` makes the missing object explicit:

```text
export an actual strict internal symmetry-breaking / selector ingredient whose content is Shannon-typed
and whose theorem-level effect is to canonically cut the QW-2191 O(2) family in a declared scope.
```

This is a *strict-core* target spec. It is not an acceptance packet and it does not introduce any new
non‑strict premise.

## Scope

`T165` is scoped only to the missing **selector ingredient** for the `QW-2191` `O(2)` ambiguity.

It does **not** decide:

1. strict-core theta export for the sigma‑int → theta lane (`T159`),
2. strict-derived slot selection for `eps` / `delta_d` (`T160/T161`),
3. any slot‑free sigma‑int → theta construction class (`T162`),
4. admissible `S_sel_int`, strict-core selector closure, or global `QW-2191` discharge,
5. ToE closure.

## Target object

If achieved, export one strict-core selector ingredient object:

```text
S_shannon_symmetry_breaking_selector_ingredient_o2_cut_v1
```

with intended meaning (typed at minimum at the level of declared strict inputs and outputs):

```text
S_shannon_symmetry_breaking_selector_ingredient_o2_cut_v1 :
  (alpha_geo_strict_derived_v1, sigma_int_strict_derived_v1, declared_mode_scaffold_v1, ...) -> o2_cut_datum_v1
```

where `o2_cut_datum_v1` is a typed object sufficient to reduce the continuous `O(2)` family from `QW-2191`
to a canonical representative (up to an explicitly declared unavoidable residual `Z2` convention).

The exact carrier of `o2_cut_datum_v1` may be chosen by the discharge (must be explicit), e.g.:

1. an angle parameter `theta_fix_v1 ∈ [0,2π)` with theorem-level uniqueness modulo the declared residual symmetry,
2. an explicit basis‑pair selector operator on the declared `QW-2190` mode scaffold selecting one representative,
3. a canonical-fixing datum reducing the admissible symmetry group from `O(2)` to a discrete residual subgroup.

## Acceptance tests (what would count as discharge)

An **actual discharge** of `S_shannon_symmetry_breaking_selector_ingredient_o2_cut_v1` must at minimum provide:

1. **Typed strict inputs:** explicitly reuse strict-source value objects at least:
   - `alpha_geo_strict_derived_v1` (`N420`),
   - `sigma_int_strict_derived_v1` (`N418`),
   and name any additional strict inputs (e.g. declared scaffold carriers such as `QW-2190`).
2. **Exported Shannon-typed objective/selector primitive:** export a **typed** strict functional/objective
   `J_shannon_v1` on the relevant candidate family whose optimization (or extremality) performs the symmetry breaking.
   - “Shannon symmetry breaking” must be represented as an explicit object, not as rhetoric.
3. **No hidden selector slot:** prove that `J_shannon_v1` (and its domain) do not silently reintroduce
   an untracked continuous knob that is equivalent to the original `O(2)` freedom.
4. **Nontriviality:** prove `J_shannon_v1` is not constant on the `O(2)` family in the declared scope.
5. **Uniqueness / O(2)-cut theorem:** prove that the induced selection is unique in the declared scope:
   - either a unique extremizer exists (mod the explicitly declared residual symmetry), or
   - the remaining ambiguity is reduced to an explicitly declared discrete residual subgroup (tracked as part of the result).
6. **Compatibility with `QW-2191`:** the discharge must explicitly state how the exported selector ingredient
   interacts with the `QW-2191` hypothesis (“kernel alone is insufficient”) and why the new ingredient is not a silent
   `QW-2192` reuse.
7. **Noncyclic + observer-free:** no `theta` inputs, no populated basis-pair inputs (`N18` discipline), and no `K_obs`
   indexing as a uniqueness source.
8. **No false pass discipline:** the discharge must not imply:
   - strict-core selector closure (`S_sel_int`),
   - global `QW-2191` discharge,
   - ToE closure,
   unless separately proved.

## Hard limits

`T165` must not claim:

1. that the strict-source Shannon-weighted refinement candidates (`F319..F322`) already constitute an internal
   selector ingredient,
2. that `alpha_geo_strict_derived_v1 := 4 ln 2` by itself supplies an `O(2)` cut,
3. that any “time-arrow selector / ontological symmetry breaking” rhetoric is already an exported strict ingredient
   (`P399`),
4. strict-core selector closure, global `QW-2191` discharge, or ToE closure.

