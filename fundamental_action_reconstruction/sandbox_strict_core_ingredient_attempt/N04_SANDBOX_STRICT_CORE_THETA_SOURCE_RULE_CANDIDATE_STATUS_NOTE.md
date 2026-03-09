# N04 Sandbox Strict-Core Theta-Source Rule Candidate Status Note

Status: `N04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_STATUS_NOTE_NO_FALSE_PASS`
As of: `2026-03-09`

## Strongest honest conclusion

The sandbox attack succeeds in one narrow sense:

```text
the sandbox now contains one conditional strict-core theta-source rule
candidate above the non-placeholder skeleton attempt
```

This is stronger than `N03` because the strict-core-only theta route is no
longer only a named skeleton.
It now also contains one explicit conditional rule:

```text
if strict core ever supplies an actual populated pair instance (u_1,u_2),
then theta_1, theta_2 are serialized by the local atan2 phase formulas
```

## What did not happen

This attack did **not** derive:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. a packet-ready strict-core minimal source skeleton,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure.

## Why this attack is still useful

It sharpens the honest frontier:

1. the sandbox now separates three levels cleanly:
   - source sketch,
   - non-placeholder skeleton attempt,
   - conditional rule candidate,
2. the exact failure point remains the same blocker from `C50`:
   strict core still does not supply the populated basis-pair instance needed
   to fire the rule.

## Honest next move

If the sandbox is continued, the next clean move is:

1. either try one incompatibility boundary showing that the conditional rule
   candidate cannot become packet-ready on current strict-core inputs,
2. or attack the missing populated basis-pair instance layer directly.
