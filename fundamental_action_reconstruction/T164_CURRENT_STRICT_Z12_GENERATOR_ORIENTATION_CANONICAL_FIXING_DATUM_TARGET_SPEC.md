# T164 Current Strict `Z_12` Generator/Orientation Canonical-Fixing Datum Target Spec

Status: `T164_CURRENT_STRICT_Z12_GENERATOR_ORIENTATION_CANONICAL_FIXING_DATUM_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

`T163` names the missing ingredient required to treat a `Z_12 -> Phase_12` embedding as *slotless*:

- either export a **canonical-fixing datum** (acceptance route `T163(4a)`), or
- export only a quotient-safe embedding and prove the downstream numeric object is invariant under the full
  `Aut(Z_12)` ambiguity (`T163(4b)`).

On the current repo state, the strict boundary theorem `N461` shows that demanding **pure**
`Aut(Z_12)`-invariance in the “no hidden generator choice” sense collapses nontrivial phase information to the
parity-only sector `{±1}` and thus cannot supply nontrivial angles.

Therefore any future attempt to obtain a **nontrivial** canonical phase embedding on the typed `AX20/T162`
lane must (at minimum) introduce an additional strict-core datum that breaks the `Aut(Z_12)` ambiguity in a
tracked way.

Additionally, `N462` makes explicit that there is no `Aut(Z_12)`-invariant way to select a generator/orientation
from the typed `Z_12` + `Aut(Z_12)` structure alone, so any nontrivial fixing datum must genuinely add (and
track) symmetry breaking rather than pretending the choice is “from air”.

`T164` names that missing datum sharply as a future-only strict-core target object with explicit acceptance
tests, so later work cannot silently treat “choose generator 1” or “choose an orientation” as already
canonical.

## Scope

`T164` is scoped only to the **generator/orientation fixing** sub-ingredient required for `T163(4a)`.

It does **not** decide:

1. discharge of `T163` (canonical/quotient-safe phase embedding),
2. any Berry/holonomy construction or gauge discipline,
3. any strict sigma-int → theta export (`T162`),
4. any `O(2)`-cut selector ingredient or `QW-2191` discharge,
5. ToE closure.

## Target object

Export one future-only strict target object:

```text
Kappa_Z12_generator_orientation_canonical_fixing_datum_target_v1
```

with intended meaning:

```text
export a strict-core internal datum D_fix_v1
that canonically fixes (in a declared scope) a preferred generator / orientation of Z_12_v1,
so that a nontrivial phase embedding Z_12_v1 -> Phase_12_v1 can be treated as slotless under T163(4a)
without smuggling in an untracked generator-choice selector.
```

## Acceptance tests (what would count as discharge)

An **actual discharge** of `Kappa_Z12_generator_orientation_canonical_fixing_datum_target_v1` must at minimum
provide:

1. **Typed prerequisite carriers:** reuse (or re-export) the typed group objects:
   - `Z_12_v1` (from `F329/N450`),
   - `Aut_Z12_v1` as the relevant discrete ambiguity group acting on embeddings (from `F331/N453`).
2. **Exported fixing datum:** an exported strict-core datum object `D_fix_v1` with:
   - explicit declared domain of validity,
   - explicit provenance classification (strict-derived or explicit strict-side premise),
   - and an explicit statement of which transformations it is invariant under.
3. **Induced generator/orientation selection:** an exported strict object that represents the “fixed”
   generator/orientation, e.g. one of the following (must be typed and explicit which is used):
   - a generator element `g_fix_v1 ∈ Z_12_v1` with proof that `⟨g_fix_v1⟩ = Z_12_v1`, or
   - an oriented 12-cycle successor map `suc_fix_v1 : I_12_v1 -> I_12_v1` with proof it is a 12-cycle and
     that it is induced by a generator.
4. **No hidden slot:** a theorem-level statement that the selection in (3) is not merely a relabeling choice:
   - the datum `D_fix_v1` must *reduce* the `Aut(Z_12)` ambiguity in a declared sense (e.g. by restricting to
     automorphisms preserving `D_fix_v1`), and
   - the selected `g_fix_v1` / `suc_fix_v1` must be invariant under that reduced admissible symmetry.
5. **Compatibility statement with `T163(4a)`:** an explicit note explaining how this datum would discharge the
   “canonical-fixing datum” branch of `T163` and what residual ambiguity (if any) remains.
6. **Noncyclic + observer-free:** no theta inputs, no populated basis-pair inputs (`N18`), and no `K_obs`
   indexing as a uniqueness source.
7. **No false pass discipline:** the discharge must not imply:
   - discharge of `T162`,
   - strict-core theta export,
   - strict-core selector closure or `QW-2191` discharge,
   unless separately proved.

## Hard limits

`T164` must not claim:

1. that a strict canonical-fixing datum is already exported,
2. that “choose generator 1” or “choose orientation” is canonical without satisfying the above acceptance
   tests,
3. that `Aut(Z_12)` ambiguity is already eliminated in strict core,
4. strict theta export or ToE closure.
