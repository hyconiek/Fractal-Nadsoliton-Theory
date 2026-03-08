# N125 Current Selector Requirement Theory-Level Acceptance Theorem

Status: `N125_DISCHARGED_CURRENT_SELECTOR_REQUIREMENT_THEORY_LEVEL_ACCEPTANCE_THEOREM_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P114`, the strongest honest theorem-level question is:

```text
what is the strongest updated-repo theorem one may now make about the
selector/symmetry-breaking requirement decision?
```

## Statement

Consider the updated repo state containing all of the following:

1. `N118`:
   selector/symmetry-breaking requirement is supported after `QW-2191`,
2. `N124`:
   no strict-core internal selector source derivation discharge is exported,
3. `AX15`:
   an explicit theory-level acceptance verdict is now present in
   `axiom_augmented_only` scope,
4. `P114`:
   the updated repo now exports an explicit theory-level selector-requirement
   decision verdict.

The theorem is:

> On the updated repo state, the selector/symmetry-breaking requirement is now
> accepted at theory level in axiom-augmented scope, while strict core remains
> unchanged.
>
> Equivalently: the project no longer leaves the selector requirement merely as
> a supported boundary. It now explicitly adopts that requirement outside
> current strict core unless and until a genuinely new strict-core internal
> selector source is later derived.

## Result

`N125` discharges:

- a theorem-level decision result for selector/symmetry-breaking requirement,
- a clean separation between `accepted in axiom-augmented scope` and
  `not derived in strict core`,
- a project-level resolution of the selector-decision indecision.

## Hard limits

`N125` does not discharge:

- strict-core selector closure,
- `QW-2191` discharge,
- legacy-to-strict bridge,
- ToE closure,
- a proof that no future strict-core source can exist.

## Recommended next move

The correct next move is now:

1. either add a genuinely new strict-core internal selector source object in a
   future route,
2. or work from the explicitly accepted axiom-augmented selector scope rather
   than pretending the decision is still open.
