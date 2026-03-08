# N128 Last Positive Branch Reduced To One Initial Source Seed Target Theorem

Status: `N128_DISCHARGED_LAST_POSITIVE_BRANCH_REDUCED_TO_ONE_INITIAL_SOURCE_SEED_TARGET_THEOREM_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P117`, the strongest honest constructive question is:

```text
what is the strongest current-repo-state theorem one may now make about the
first future construction move?
```

## Statement

Consider the current repo state containing all of the following:

1. `N127`:
   the last positive branch is reduced to one full future chain,
2. `F31`:
   the first forced seed target is `S_sel_int -> E_orient`,
3. `P117`:
   the current repo now reduces the last positive branch to that one initial
   seed target.

The theorem is:

> On the current repo state, the first future construction move is reduced to
> one explicit initial source-seed target:
>
> `S_sel_int -> E_orient`
>
> while the downstream remainder
> `B_sel -> R_sel -> O_sel`
> stays explicitly open for later.

## Result

`N128` discharges:

- a theorem-level reduction of the first future move to one explicit initial
  seed target,
- a clean separation between what must be built first and what remains
  downstream.

## Hard limits

`N128` does not discharge:

- existence of the seed,
- downstream bridge/reduction/operator completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
