# P17 Sandbox Actual Populated Basis-Pair Instance Semantic Blocker Attack Probe

Status: `P17_SANDBOX_ACTUAL_POPULATED_BASIS_PAIR_INSTANCE_SEMANTIC_BLOCKER_ATTACK_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the populated-instance blocker now has one explicit gate audit
while actual population still remains absent.

## What is checked

`P17` checks whether:

1. minimal basis-pair export skeleton is present,
2. conditional populated-instance schema is present,
3. dedicated populated-instance candidate file is present,
4. theta input slots are written,
5. actual theta input values remain absent,
6. actual populated instance remains absent,
7. the blocker is therefore narrowed below actual theta inputs.

## Result matrix

### Minimal basis-pair export skeleton present

Current verdict after `F17`:

```text
YES
```

### Conditional populated-instance schema present

Current verdict after `F17`:

```text
YES
```

### Dedicated populated-instance candidate file present

Current verdict after `F17`:

```text
YES
```

### Theta input slots written

Current verdict after `F17`:

```text
YES
```

### Actual theta input values present

Current verdict after `F17`:

```text
NO
```

### Actual populated basis-pair instance present

Current verdict after `F17`:

```text
NO
```

### Populated-instance blocker narrowed below actual theta inputs

Current verdict after `F17`:

```text
YES
```

## Hard limits

`P17` does not establish:

1. actual populated basis-pair instance,
2. actual `u_1`, `u_2`,
3. actual `theta_1`, `theta_2`,
4. actual strict-core theta source supply,
5. actual bridge discharge,
6. actual provider emission,
7. actual `E_orient`,
8. admissible `S_sel_int`,
9. actual strict-core selector closure,
10. actual ToE closure.
