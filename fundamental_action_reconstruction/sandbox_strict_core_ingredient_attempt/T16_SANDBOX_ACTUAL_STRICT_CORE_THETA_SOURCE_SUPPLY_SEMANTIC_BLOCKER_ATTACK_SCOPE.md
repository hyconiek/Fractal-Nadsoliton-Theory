# T16 Sandbox Actual Strict-Core Theta-Source Supply Semantic Blocker Attack Scope

Status: `T16_SANDBOX_ACTUAL_STRICT_CORE_THETA_SOURCE_SUPPLY_SEMANTIC_BLOCKER_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F15/P15/N15`, the last remaining semantic blocker under the sandbox
bridge-discharge lane is explicitly:

```text
actual_strict_core_theta_source_supply_present := no
```

The next direct question is:

```text
can this last semantic blocker be attacked directly and narrowed without
pretending that actual theta supply already exists?
```

## Repo-consistent support reused

Current repo and sandbox state already give:

1. `C33`
   - local phase formula class is present,
2. `C34`
   - local representative class is present,
3. `C49`
   - conditional populated-instance schema is present,
4. `F04`
   - conditional phase-serialization rule candidate is explicit,
5. `F15`
   - the bridge-discharge gate is narrowed to the single remaining semantic
     blocker above.

## Intended move

`T16` does not try to turn:

```text
actual_strict_core_theta_source_supply_present := no
```

into `yes`.

It attacks that blocker in a narrower way:

1. expose all positive semantic support already present below supply,
2. add one dedicated source-supply candidate file,
3. isolate the live negative point to:
   - no actual populated basis-pair instance,
   - therefore still no actual theta values supplied.

## Hard limits

`T16` must not claim:

1. actual strict-core theta source supply,
2. actual `theta_1`, `theta_2`,
3. actual populated `u_1`, `u_2`,
4. actual bridge discharge,
5. actual provider emission,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. ToE closure.
