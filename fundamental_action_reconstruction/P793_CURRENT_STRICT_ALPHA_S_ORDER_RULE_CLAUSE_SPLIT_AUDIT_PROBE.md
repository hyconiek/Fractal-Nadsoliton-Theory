# P793 Current Strict Alpha_s Order-Rule Clause-Split Audit Probe

Status: `P793_CURRENT_STRICT_ALPHA_S_SOURCE_AUTHORITY_CANDIDATE_SUPPORTED_NORMALIZATION_BOUNDARY_BLOCKED`
As of: `2026-03-19`

## Goal

After `P792/F792`, the next honest question is:

```text
which clause of the current probe-local order rule
is actually the sharp blocker to export?
```

## Scope

`P793` does not export the order rule.
It only audits the two leading clauses from `P792` separately:

1. `source_authority`,
2. `normalization_boundary`.

## Main Checks

1. test whether `source_authority` already has strong strict-side support from
   the current Release-7 observable hierarchy,
2. test whether `normalization_boundary` already has any exported strict-side
   rule beyond the probe-local ranking used in `P792`,
3. keep geometric-mean language on the non-strict calibration side unless a
   new strict theorem exists,
4. identify which clause is now the sharp blocker to export-grade authority.

## Result

`P793` finds an asymmetric outcome:

1. `source_authority` is **candidate-supported** on the current repo state,
2. `normalization_boundary` remains **blocked / nonexport**.

So the blocker narrows again:

```text
not "the whole order rule is equally missing",
but "the normalization-boundary clause is now the sharp blocker"
```

## Hard Limit

`P793` does not claim that `source_authority_rule_ref` is already exported.
It only says that the current repo evidence supports it much more strongly
than the still-missing `normalization_boundary_rule_ref`.
