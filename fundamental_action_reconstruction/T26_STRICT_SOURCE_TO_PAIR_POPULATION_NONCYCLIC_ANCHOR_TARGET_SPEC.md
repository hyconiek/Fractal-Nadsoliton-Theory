# T26 Strict Source-To-Pair-Population Noncyclic Anchor Target Spec

Status: `T26_CURRENT_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N275`, the exact closure frontier still contains:

1. one genuine strict-side selector ingredient,
2. one basis-independent quotient-safe promotion/discharge layer,
3. one actual non-strict declared-scope discharge ingredient if that lane is
   pursued.

After `N283`, the official strict-side lane is additionally frozen on the same
`strict_extension_only` ladder for the remaining four `F34` clauses.

After sandbox `N18`, the strongest available strict-core attempt also freezes
at one explicit dependency loop:

```text
theta supply waits for populated basis-pair instance
populated basis-pair instance waits for theta inputs
```

So the next honest strict-side solution candidate is no longer:

```text
one more same-lane extension lift
or
one more same-loop sandbox positive move
```

It is narrower:

```text
one source-side, observer-free, kernel-split-safe, pair-indexed,
noncyclic anchor target
which could in principle break both current blocker cuts
without reactivating K_legacy_ont
and without promoting K_obs to primary selector source
```

`T26` specifies only that target.

It does not discharge:

1. an actual strict-core ingredient,
2. actual theta values,
3. actual populated basis-pair instance,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.

## Reused support

`T26` is built only from already exported current-state material:

1. `N256`
   - source-side basis-independent promotion is actual,
2. `N257`
   - source-side quotient-safe `QW-2191` resolution is actual,
3. `N258`
   - one actual declared-scope source-topology selector theorem is exported,
4. `N275`
   - one genuine strict-side selector ingredient is still missing,
5. `N282`
   - present closure stack is frozen at current-state incompatibility boundary,
6. `N283`
   - the remaining four strict-side clauses are nonentering on the same
     official extension lane,
7. `H15/N13`
   - `K_obs` remains a distinct hypothesis-extension lane and is not
     identified with selector-facing strict-core source,
8. `N234`
   - downstream observer stability may not be promoted to selector closure,
9. `N18`
   - the sandbox strict-core route is nonentering under the same loop and
     needs a noncyclic anchor from outside that loop.

## Formal target

Freeze one future-only strict-side target:

```text
V_strict_src_pair_population_anchor_target_v1 :
  tau_src_candidate_v1 -> strict_src_pair_population_noncyclic_anchor_target_v1
```

with component packet:

```text
Gamma_strict_src_pair_population_anchor_components_target_v1 :=
(
  strict_source_orientation_seed_target_v1,
  pair_indexed_population_anchor_target_v1,
  strict_internal_orientation_provider_target_v1
)
```

## Meaning of the three target components

### 1. `strict_source_orientation_seed_target_v1`

This component is intended only as:

```text
one future source-side seed
derived from already exported source-topology selector support
but still below actual theta population
```

It must remain:

1. upstream of observer,
2. independent of `K_obs`,
3. below actual `E_orient`.

### 2. `pair_indexed_population_anchor_target_v1`

This component is intended only as:

```text
one future pair-indexed anchor
for at least the minimal designated pair family [pair1, pair2]
which supplies a noncyclic entry point
for populated basis-pair construction
```

It must remain:

1. pair-indexed rather than only class-level,
2. external to the `theta-supply <-> populated-instance` loop,
3. below actual populated basis-pair discharge.

### 3. `strict_internal_orientation_provider_target_v1`

This component is intended only as:

```text
one future downstream provider route
from the noncyclic pair-population anchor
toward actual strict internal orientation data
```

It is still weaker than:

1. actual `E_orient`,
2. downstream `B_sel`,
3. downstream `R_sel`,
4. downstream `O_sel`,
5. admissible `S_sel_int`.

## Intended route shape

If this target were ever later discharged, the strict-side route would have to
take the shape:

```text
tau_src_candidate_v1
  -> strict_source_orientation_seed_target_v1
  -> pair_indexed_population_anchor_target_v1
  -> strict_internal_orientation_provider_target_v1
```

The target is scientifically honest only if all of the following remain
explicit:

1. the route is source-side and observer-free,
2. the route does not reactivate `K_legacy_ont` as live constructive kernel,
3. the route does not use `K_obs` as primary selector source,
4. the route is future-only,
5. the route is still below actual strict-core admission and closure.

## Why this target is the correct proposal

This target is proposed because it is the smallest route that could in
principle answer both current blockers at once:

1. `N283`
   - by adding a genuinely new strict-core ingredient class instead of one
     more same-lane extension lift,
2. `N18`
   - by providing one noncyclic anchor from outside the exposed
     `theta-supply <-> populated-instance` loop.

It is therefore a solution target, not yet a solved solution.

## Hard limits

`T26` does not claim:

1. actual anchor discharge,
2. actual theta supply,
3. actual populated basis-pair instance,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. global selector closure,
8. global `QW-2191` discharge,
9. ToE closure,
10. any theorem erasing the `bridge/non-bridge` frontier.
