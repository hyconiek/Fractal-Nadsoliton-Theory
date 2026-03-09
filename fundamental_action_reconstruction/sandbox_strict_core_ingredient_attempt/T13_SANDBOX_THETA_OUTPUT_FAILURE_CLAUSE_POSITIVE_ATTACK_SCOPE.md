# T13 Sandbox Theta-Output Failure Clause Positive Attack Scope

Status: `T13_SANDBOX_THETA_OUTPUT_FAILURE_CLAUSE_POSITIVE_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F12/P12/N12`, all four provider-emission failure clauses are explicitly
structured, but none has yet received a positive partial lift.

The most honest next clause to attack positively is the theta-output clause,
because current repo state already exports:

1. `C48`
   - a minimal basis-pair export skeleton,
2. `C49`
   - a packet-ready conditional populated-instance schema,
3. `C50`
   - an exact boundary showing that the still-missing piece is the strict-core
     source supply of actual `theta_1`, `theta_2`.

## Direct question

Can the theta-output failure clause be lifted from:

```text
pure absence
```

to:

```text
positive support present below actual values
```

without claiming actual theta-output emission?

## Intended move

`T13` will not turn:

```text
actual_theta_output_values_present := no
```

into `yes`.

It will attack the same clause positively in a narrower way:

1. lock the output slot names,
2. register that the minimal basis-pair export skeleton is already present,
3. register that the conditional populated-instance schema is already present,
4. isolate the remaining negative point to:
   - no strict-core theta source skeleton,
   - therefore no actual theta values.

## Hard limits

`T13` must not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual provider emission,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
