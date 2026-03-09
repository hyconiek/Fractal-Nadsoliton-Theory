# T08 Sandbox Provider Creation-Status Field Attack Scope

Status: `T08_SANDBOX_PROVIDER_CREATION_STATUS_FIELD_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F07/P07/N07`, the sandbox already has one pair-indexed
provider-carrier candidate, but its `creation_status` field is still too
coarse:

```text
creation_status := carrier_not_created
```

The next narrow move is:

```text
refine only that single field into one explicit readiness/non-readiness
verdict tied to the actual sandbox filesystem state
```

## Support reused

1. `F07`
   - pair-indexed provider-carrier candidate,
2. `C43`
   - filename/path convention discipline,
3. `C44`
   - minimal template content discipline,
4. `C45`
   - non-destructive creation admission requires an existing carrier lane,
5. current sandbox filesystem state
   - `sandbox_strict_core_ingredient_attempt/generated/` is absent.

## Question

Is it honest to refine `creation_status` into a structured verdict saying:

1. filename/path convention is specified,
2. minimal template content is specified,
3. carrier directory is currently absent,
4. therefore non-destructive carrier creation is not yet admissible on the
   present sandbox state,
5. and the file is still not created?

## Hard limits

`T08` must not claim:

1. actual provider carrier creation,
2. actual provider emission,
3. actual `theta_1`, `theta_2`,
4. actual populated `u_1`, `u_2`,
5. actual internal orientation datum,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. ToE closure.
