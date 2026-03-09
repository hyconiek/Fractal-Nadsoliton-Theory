# T17 Sandbox Actual Populated Basis-Pair Instance Semantic Blocker Attack Scope

Status: `T17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F16/P16/N16`, the last live negative point under the theta-source lane
is explicitly:

```text
actual_populated_basis_pair_instance_present := false
```

The next direct question is:

```text
can this populated-instance blocker be attacked directly and narrowed without
pretending that an actual populated instance already exists?
```

## Repo-consistent support reused

Current repo and sandbox state already give:

1. `C48`
   - minimal basis-pair export skeleton is present,
2. `C49`
   - conditional populated-instance schema is present,
3. `F05`
   - the populated-instance lane already has an explicit artifact schema,
4. `F16`
   - the theta-source lane is now narrowed exactly down to this blocker.

## Intended move

`T17` does not try to turn:

```text
actual_populated_basis_pair_instance_present := false
```

into `true`.

It attacks that blocker in a narrower way:

1. add one dedicated populated-instance candidate file,
2. make the theta-input slots explicit inside that file,
3. isolate the live negative point to:
   - actual theta input values are still absent,
   - therefore population still does not happen.

## Hard limits

`T17` must not claim:

1. actual populated basis-pair instance,
2. actual `u_1`, `u_2`,
3. actual `theta_1`, `theta_2`,
4. actual strict-core theta source supply,
5. actual bridge discharge,
6. actual provider emission,
7. actual `E_orient`,
8. admissible `S_sel_int`,
9. strict-core selector closure,
10. ToE closure.
