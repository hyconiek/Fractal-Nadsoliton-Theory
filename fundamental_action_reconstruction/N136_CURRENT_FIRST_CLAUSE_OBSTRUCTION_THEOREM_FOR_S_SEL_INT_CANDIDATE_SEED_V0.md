# N136 Current First Clause Obstruction Theorem For S_sel_int_candidate_seed_v0

Status: `N136_EXECUTABLE_THEOREM_PACKET_READY`
As of: `2026-03-08`

## Theorem statement

Within the current repo state, `S_sel_int_candidate_seed_v0` does not yet
discharge the first admissibility clause

```text
genuinely_new_strict_core_source_object_required
```

for future admissible `S_sel_int`.

## Scope

This theorem is scoped only to:

- current repo state,
- `S_sel_int_candidate_seed_v0`,
- the first admissibility clause.

## Proof basis

1. `F36` defines `S_sel_int_candidate_seed_v0` as a packaged pair built from
   current route ingredients.
2. `F34` requires a genuinely new strict-core source object for admissible
   `S_sel_int`.
3. `P125` checks whether the candidate seed has crossed that threshold.

## Hard limits

This theorem does not claim:

- that no future genuinely new source object can be constructed,
- that later admissibility clauses are resolved,
- that admissible `S_sel_int` is impossible forever,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
