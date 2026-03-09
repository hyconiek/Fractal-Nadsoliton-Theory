# P09 Sandbox Carrier-Directory Precondition Attack Probe

Status: `P09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Check whether the carrier-directory precondition is now genuinely cleared
while file/provider creation remains absent.

## What is checked

`P09` checks whether the new step:

1. clears the directory-level blocker identified by `F08`,
2. remains below provider-file creation,
3. remains below provider emission,
4. remains below actual `theta_1`, `theta_2` export.

## Result matrix

### Carrier-directory precondition cleared

Current verdict after `F09`:

```text
YES
```

Reason:

1. the sandbox `generated/` directory now exists,
2. the exact blocker from `F08` no longer applies.

### Provider file created

Current verdict after `F09`:

```text
NO
```

Reason:

1. only the directory was created,
2. no candidate JSON provider file exists yet.

### Provider emission performed

Current verdict after `F09`:

```text
NO
```

Reason:

1. no provider file was created,
2. no emitted provider instance exists.

### Actual strict-core `theta_1`, `theta_2` export

Current verdict after `F09`:

```text
NO
```

Reason:

1. clearing the carrier-directory precondition is a carrier-side move only,
2. it does not alter the missing provider logic itself.

## Hard limits

`P09` does not establish:

1. actual provider file creation,
2. actual provider emission,
3. actual `theta_1`, `theta_2`,
4. actual populated `u_1`, `u_2`,
5. actual internal orientation datum,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. actual strict-core selector closure,
9. actual ToE closure.
