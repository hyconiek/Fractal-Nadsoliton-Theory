# T09 Sandbox Carrier-Directory Precondition Attack Scope

Status: `T09_SANDBOX_CARRIER_DIRECTORY_PRECONDITION_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F08/P08/N08`, the sandbox already knows the sharp blocker on the
provider-carrier side:

```text
carrier directory absent
```

The next narrow move is:

```text
clear exactly that one precondition by creating the sandbox carrier
directory, without creating any provider file and without emitting any
provider instance
```

## Support reused

1. `F07`
   - pair-indexed provider-carrier candidate,
2. `F08`
   - refined `creation_status` showing that the missing directory is the live
     blocker,
3. `C45`
   - non-destructive additive carrier steps are methodology-sensitive,
4. current sandbox isolation rule
   - the whole subtree remains removable.

## Question

Is it honest to execute the following narrow additive move:

```text
create
  fundamental_action_reconstruction/sandbox_strict_core_ingredient_attempt/generated/
```

while preserving:

1. no provider file created,
2. no provider emitted,
3. no actual `theta_1`, `theta_2`,
4. no strict-core closure claim?

## Hard limits

`T09` must not claim:

1. actual provider carrier creation at file level,
2. actual provider emission,
3. actual `theta_1`, `theta_2`,
4. actual populated `u_1`, `u_2`,
5. actual internal orientation datum,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. ToE closure.
