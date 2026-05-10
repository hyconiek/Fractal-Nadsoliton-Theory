# P815 Current Strict Alpha_s Domain Member Grade Handling Clause Export Admitted Selection Or Preference Rule Schema Still Blocked

Status: `P815_CURRENT_STRICT_ALPHA_S_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED`
As of: `2026-03-19`

## Goal

After `F814`, the next honest question is:

```text
now that the F813 support domain is explicit
and the schema-only target is frozen,
can the current repo already support
one exact domain_member_grade_handling_clause,
or does even that clause remain implicit?
```

## Scope

`P815` does not export a selection/preference schema.
It only splits the `F814` blocker into:

1. `domain_member_grade_handling_clause`,
2. `selection_or_preference_rule_schema`.

## Main Checks

1. confirm `F814` already freezes the exact schema-only target and requires
   both `domain_member_grade_handling_clause` and
   `selection_or_preference_rule_schema`,
2. confirm `F813` already exports the exact finite support domain with explicit
   member grades,
3. confirm `F321` keeps the current Shannon support at candidate-only grade,
4. confirm `P755/P756` keep the lawful future entry support at future-only
   target grade below actual realization,
5. decide whether one exact grade-handling clause may already be exported from
   those current facts without selecting one source,
6. keep explicit that the selection/preference rule schema itself remains
   unexported.

## Result

`P815` returns a split verdict:

```text
domain_member_grade_handling_clause export is admitted,
but selection_or_preference_rule_schema remains blocked
```

Therefore the blocker narrows again:

```text
not "how are the two support grades handled",
but "what exact schema acts on that already grade-disciplined domain
to choose or prefer one source on the F813 route"
```

## Hard Limit

`P815` does not claim:

1. that any source has already been selected,
2. that any source-binding relation is already exported,
3. that any adapter action schema is already exported,
4. that `F321` is promoted beyond candidate-only grade,
5. that `T209/P755` are promoted beyond future-only target grade,
6. alpha_s boundary export readiness,
7. QCD closure,
8. ToE closure.
