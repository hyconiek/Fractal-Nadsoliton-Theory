# P04 Sandbox Strict-Core Theta-Source Rule Candidate Probe

Status: `P04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the new rule candidate is honestly stronger than `F03` while
still consistent with `C49/C50`.

## What is checked

`P04` checks whether the new candidate:

1. stays strict-core-only,
2. upgrades the sandbox from a bare skeleton attempt to a conditional rule
   form,
3. avoids importing the axiom-augmented theta branch,
4. avoids claiming actual populated instances or actual exported phases,
5. remains below a packet-ready strict-core minimal source skeleton.

## Result matrix

### Conditional strict-core theta-source rule candidate

Current verdict after `F04`:

```text
YES
```

Reason:

1. `F04` explicitly states the rule that populated `u_1,u_2` would serialize
   phases through the `C33` formulas,
2. the rule is anchored to the `C49` conditional populated-instance schema,
3. therefore the sandbox is sharper than a skeleton attempt alone.

### Import of axiom-augmented theta branch

Current verdict after `F04`:

```text
NO
```

Reason:

1. the rule uses only strict-core-side formulas and the conditional schema,
2. it does not use the non-strict `theta_1^*`, `theta_2^*` branch.

### Actual populated strict-core basis-pair instance

Current verdict after `F04`:

```text
NO
```

Reason:

1. `C49/C50` still leave population blocked,
2. `F04` preserves that blocker explicitly.

### Actual strict-core `theta_1`, `theta_2` export

Current verdict after `F04`:

```text
NO
```

Reason:

1. `C50` still denies a packet-ready strict-core source skeleton,
2. `F04` does not override that denial.

### Packet-ready strict-core minimal source skeleton

Current verdict after `F04`:

```text
NO
```

Reason:

1. the rule candidate is conditional on a populated instance that strict core
   does not export,
2. therefore it stays below the strong sourcehood denied by `C50`.

## Hard limits

`P04` does not establish:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. actual strict-core selector closure,
7. actual ToE closure.
