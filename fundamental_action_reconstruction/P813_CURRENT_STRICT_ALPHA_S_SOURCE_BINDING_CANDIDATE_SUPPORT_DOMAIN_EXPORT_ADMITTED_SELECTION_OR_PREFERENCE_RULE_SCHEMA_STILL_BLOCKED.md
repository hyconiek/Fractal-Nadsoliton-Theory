# P813 Current Strict Alpha_s Source-Binding Candidate Support Domain Export Admitted Selection Or Preference Rule Schema Still Blocked

Status: `P813_CURRENT_STRICT_ALPHA_S_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED`
As of: `2026-03-19`

## Goal

After `F812`, the next honest question is:

```text
can the current repo already export
one exact finite candidate support domain
for the source-binding rule target,
or does even that domain still remain implicit?
```

## Scope

`P813` does not export a selection rule.
It only splits the `F812` blocker into:

1. `candidate_source_support_domain_ref`,
2. `selection_or_preference_rule_schema`.

## Main Checks

1. confirm `F812` already freezes the exact rule target and requires both
   `candidate_source_support_domain_ref` and `selection_or_preference_rule_schema`,
2. confirm the current source-side support pair is already explicit on current
   exports,
3. confirm the future-only support target remains future-only and not actual
   entry,
4. decide whether that exact support pair may already be exported as one finite
   support-domain object,
5. keep explicit that the selection/preference rule schema itself remains
   unexported.

## Result

`P813` returns a split verdict:

```text
candidate support domain export is admitted,
but selection_or_preference_rule_schema remains blocked
```

Therefore the blocker narrows again:

```text
not "what is the support domain",
but "what exact rule schema acts on that support domain
to choose or prefer one source for the F811 route"
```

## Hard Limit

`P813` does not claim:

1. that any source has already been selected,
2. that any source-binding relation is already exported,
3. that any adapter action schema is already exported,
4. that `F321` is promoted beyond candidate-only grade,
5. that `T209` is promoted beyond future-only entry-target grade,
6. alpha_s boundary export readiness,
7. QCD closure,
8. ToE closure.
