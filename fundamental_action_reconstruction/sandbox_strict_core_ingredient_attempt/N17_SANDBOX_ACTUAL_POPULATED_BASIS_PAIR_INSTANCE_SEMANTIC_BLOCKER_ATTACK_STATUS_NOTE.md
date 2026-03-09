# N17 Sandbox Actual Populated Basis-Pair Instance Semantic Blocker Attack Status Note

Status: `N17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_STATUS_NOTE_NO_FALSE_PASS`
As of: `2026-03-09`

## Strongest honest conclusion

The sandbox attack succeeds in one narrow but real sense:

```text
the actual populated basis-pair instance blocker is no longer a flat absence
claim; it now has one explicit semantic gate audit and one dedicated
populated-instance candidate file below actual population
```

This is stronger than `N16` because the blocker is no longer represented only
as:

```text
actual_populated_basis_pair_instance_present := false
```

It now has a finer internal split:

1. minimal basis-pair export skeleton is present,
2. conditional populated-instance schema is present,
3. dedicated populated-instance candidate file is present,
4. theta input slots are written,
5. but actual theta input values remain absent,
6. therefore actual populated instance remains absent.

## Loop exposure

This attack makes the circular frontier explicit rather than hidden:

1. source supply waits for actual populated instance,
2. actual populated instance waits for actual theta inputs.

## What did not happen

This attack did **not** derive:

1. actual populated basis-pair instance,
2. actual `u_1`, `u_2`,
3. actual `theta_1`, `theta_2`,
4. actual strict-core theta source supply,
5. actual bridge discharge,
6. actual provider emission,
7. actual `E_orient`,
8. admissible `S_sel_int`,
9. strict-core selector closure.

## Why this attack is still useful

It sharpens the frontier to its honest fixed point:

1. both sides of the dependency are now explicit,
2. almost no hidden semantic slack remains,
3. the next honest move is no longer another blind positive lift.

## Honest next move

If the sandbox is continued, the clean next move is:

1. either write the incompatibility boundary for this exposed dependency loop,
2. or leave the sandbox and return to a different noncyclic blocker-cut.
