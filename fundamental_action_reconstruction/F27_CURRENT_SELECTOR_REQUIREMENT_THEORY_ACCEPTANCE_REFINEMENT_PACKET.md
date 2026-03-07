# F27 Current Selector Requirement Theory-Acceptance Refinement Packet

Status: `F27_EXECUTED_CURRENT_SELECTOR_REQUIREMENT_THEORY_ACCEPTANCE_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P108/N118`, one higher-order project gap remains explicit:

```text
explicit_theory_level_acceptance_of_selector_or_symmetry_breaking_requirement_if_no_internal_source_is_derived
```

`F27` asks the next honest refinement question:

```text
what are the two narrowest ways that this remaining theory-level gap could in
principle be discharged on the current repo state?
```

## Result

`F27` establishes the following refinement:

1. one possible discharge would be an explicit theory-level `acceptance`
   verdict adopting the selector/symmetry-breaking requirement into the
   axiom-augmented theory scope if no internal source is derived,
2. the other possible discharge would be an explicit theory-level `deferral`
   verdict keeping the selector/symmetry-breaking requirement as an active
   boundary while the project continues to seek an internal selector source,
3. the current repo exports neither decision verdict yet.

## Real reduction after F27

So the missing object is no longer:

```text
one vague theory-level acceptance gap after QW-2191
```

It is now:

```text
two narrower decision branches
```

namely:

1. `explicit_theory_level_acceptance_verdict_adopting_the_selector_or_symmetry_breaking_requirement_into_axiom_augmented_theory_scope_if_no_internal_source_is_derived`
2. `explicit_theory_level_deferral_verdict_keeping_the_selector_or_symmetry_breaking_requirement_as_an_active_boundary_while_internal_source_search_continues`

## Why this follows

The split is forced by the current repo state:

1. `P108/N118` already support the requirement conclusion for `QW-2191`,
2. `S2` already records that formalizing the selector requirement is one of
   the highest-level project priorities,
3. but the repo still exports no explicit theory-level decision verdict,
4. therefore the narrowest honest refinement is exactly:
   acceptance branch vs deferral branch.

## What F27 does not claim

`F27` does not claim:

- that the selector/symmetry-breaking requirement is already accepted into the
  theory,
- that the selector/symmetry-breaking requirement is already formally deferred
  by the theory,
- that strict-core selector closure is achieved,
- `QW-2191` discharge,
- a legacy-to-strict kernel bridge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. test whether the current repo already exports either the acceptance branch
   or the deferral branch,
2. or formalize theorem-level that both decision branches remain absent.
