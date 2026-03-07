# N119 Current Selector Requirement Theory-Level Decision Boundary Theorem

Status: `N119_DISCHARGED_CURRENT_SELECTOR_REQUIREMENT_THEORY_LEVEL_DECISION_BOUNDARY_THEOREM_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P109`, the strongest honest theorem-level question is:

```text
what is the strongest current-repo-state statement one may now make about the
theory-level decision status of the selector/symmetry-breaking requirement?
```

## Statement

Consider the current repo state containing all of the following:

1. `P108/N118`:
   the selector/symmetry-breaking requirement conclusion is supported for the
   `QW-2191` frontier,
2. `S2`:
   the requirement question is already elevated to top-level project priority,
3. `F27/P109`:
   the open theory-level decision has been refined into acceptance vs deferral,
   and neither branch is currently exported.

The theorem is:

> On the current repo state, the repo supports the selector/symmetry-breaking
> requirement conclusion for `QW-2191`, but exports neither an explicit
> theory-level acceptance verdict nor an explicit theory-level deferral verdict.
>
> Therefore the project has crossed the requirement-support boundary, but has
> not yet crossed the theory-level decision boundary.

Equivalently:

> The current repo does justify saying that kernel alone is insufficient and
> that an explicit selector/symmetry-breaking premise is the active boundary,
> but it does **not** yet justify saying that the theory has formally adopted
> that premise.

## Result

`N119` discharges:

- a current-repo-state theorem-level decision-boundary result for the
  selector/symmetry-breaking requirement after `QW-2191`,
- a theorem-level warning against citing `N118` as if the theory had already
  accepted the selector requirement,
- a clean handoff to the remaining decision branches:
  explicit acceptance verdict vs explicit deferral verdict.

## Hard limits

`N119` does not discharge:

- an explicit theory-level acceptance verdict,
- an explicit theory-level deferral verdict,
- strict-core selector closure,
- `QW-2191` discharge,
- a legacy-to-strict kernel bridge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. either formalize explicit theory-level acceptance of the selector or
   symmetry-breaking requirement,
2. or formalize explicit theory-level deferral while internal-source search
   continues,
3. without pretending that the decision is already made.
