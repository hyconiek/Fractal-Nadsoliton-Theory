# P167 Current Additive Preobserver Source Object Nonreduction Probe

Status: `P167_EXECUTED_CURRENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_NONREDUCTION_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Test one narrow current-state claim:

```text
does the repo now export one explicit nonzero nonreduction witness
showing that S_preLM_additive_candidate_v1 is not equal to the same-basis
packaged F75 target realization?
```

## Probe rule

The probe may return a positive result only if all of the following are true:

1. `F76` still exports `S_preLM_additive_candidate_v1`,
2. `F79` exports one explicit same-basis packaged realization of the `F75`
   target at `d_* = 1`,
3. the difference vector is explicitly nonzero,
4. the conclusion is stated only as `nonreduction`, not as clause satisfaction.

## Allowed conclusion

If the probe passes, the only allowed conclusion is:

```text
the repo exports one explicit nonzero nonreduction witness between
S_preLM_additive_candidate_v1 and the same-basis packaged realization of the
F75 target
```

No stronger conclusion is allowed.
