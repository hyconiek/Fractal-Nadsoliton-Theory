# T163 Current Strict AX20 Typed `Z_12` Phase-Embedding Canonicity (Slotless) Target Spec

Status: `T163_CURRENT_STRICT_AX20_TYPED_Z12_PHASE_EMBEDDING_CANONICITY_SLOTLESS_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `F329/N450`, the repo exports a typed cyclic group carrier `Z_12_v1` and a typed regular action on the
12-slot scaffold.

`P415/P416` keep explicit that the next missing typed step for any serious “slotless projector” (`AX20` / `T162`)
attempt is not merely “a `Z_12` exists”, but:

```text
a canonical (slotless) phase embedding of the 12-cycle,
with no hidden offset/scale/generator-choice selector.
```

`T163` names that missing ingredient sharply as a future-only strict-core target object, with explicit
acceptance tests, so later work cannot silently treat “pick 2π/12 phases” as already canonical.

This target is **below** `T162` (it is a prerequisite sub-ingredient for any `Z_12`-phase-based slot-free
construction class).

## Current state note (canonicity is still blocked)

On the current repo state, the repo exports strict boundary theorems clarifying why “just use `2π/12` phases”
cannot be treated as canonical:

1. a canonical oriented 12-cycle successor map is not exported for the abstract phase carrier (`N456`), and
2. purely Aut(`Z_12`)-invariant phase-embedding/holonomy claims collapse to the parity-only sector `{±1}` and
   therefore yield only trivial angles `{0, π}` (`N461`),
3. there is no Aut(`Z_12`)-invariant way to canonically pick a generator/orientation from the typed
   `Z_12` + `Aut(Z_12)` structure alone (`N462`).

So an actual discharge of `T163` must explicitly add a new canonical-fixing datum (acceptance route 4a; named
as `T164`) or
prove a quotient-safe embedding whose *downstream* use is invariant under the full `Aut(Z_12)` ambiguity
(acceptance route 4b), without smuggling in an untracked generator/orientation choice.

## Scope

`T163` is scoped only to phase-embedding canonicity for the typed `Z_12` carrier.

It does **not** decide:

1. any density-operator rigidity forcing eigenvalues `1/2`,
2. any Berry/holonomy construction and gauge discipline,
3. any strict sigma-int → theta export,
4. any `O(2)`-cut selector ingredient or `QW-2191` discharge,
5. ToE closure.

## Target object

Export one future-only strict target object:

```text
Phi_Z12_phase_embedding_canonical_or_quotient_safe_slotless_target_v1
```

with intended meaning:

```text
Phi_Z12_phase_embedding_canonical_or_quotient_safe_slotless_target_v1 :
  Z_12_v1  ->  Phase_12_v1
```

where `Phase_12_v1` is a typed 12-phase carrier (e.g. a discrete `U(1)`-like object of 12th roots),
and the embedding is *canonical* in the strict sense below.

## Acceptance tests (what would count as discharge)

An **actual discharge** must at minimum provide:

1. **Typed carrier present:** an exported strict `Z_12_v1` carrier object and action (already present via
   `F329/N450`; may be reused).
2. **Typed phase carrier:** an exported strict finite phase carrier object `Phase_12_v1` of cardinality 12,
   with explicit group law / composition rule.
3. **Typed embedding:** an exported strict map object:
   ```text
   emb_v1 : Z_12_v1 -> Phase_12_v1
   ```
   satisfying the typed homomorphism law.
4. **No hidden generator/offset/scale slot:** the construction must **not** rely on an untracked selector such
   as:
   - choosing a preferred generator of `Z_12_v1`,
   - choosing an orientation of the 12-cycle,
   - choosing an arbitrary phase offset.

   This must be resolved in one of two acceptable strict ways:
   - **(4a) Canonical-fixing datum:** export a strict internal datum (see `T164`) that canonically fixes the generator/orientation
     (and prove the fixing is invariant under the admissible symmetries), **or**
   - **(4b) Quotient-safe embedding:** export the embedding only up to the relevant automorphism/gauge group and
     prove that the downstream quantities used in later steps (Berry/holonomy/theta) are invariant under that
     equivalence (so the hidden slot is removed by quotient invariance).
5. **Noncyclic + observer-free:** no theta inputs, no populated instances, and no `K_obs`-indexed selection.
6. **No false pass discipline:** no implied discharge of `T162`, no implied `O(2)` cut, no implied `QW-2191`
   discharge.

## Hard limits

`T163` must not claim:

1. that such a canonical/quotient-safe phase embedding is already exported,
2. that “`2π/12`” or “12th roots” are canonical without satisfying (4),
3. strict theta export or strict-core selector closure,
4. `QW-2191` discharge,
5. ToE closure.
