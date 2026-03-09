# T07 Sandbox Provider Object Field Refinement Scope

Status: `T07_SANDBOX_PROVIDER_OBJECT_FIELD_REFINEMENT_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F06/P06/N06`, the sandbox already has one upstream provider artifact
schema, but its `provider_object` field is still semantically coarse:

```text
strict_core_theta_input_provider_for_pair1_pair2
```

The next narrow move is:

```text
refine only that single field into one pair-indexed provider-carrier
candidate, without claiming provider emission
```

## Support reused

1. `F06`
   - upstream provider artifact schema,
2. `C42`
   - persisted carrier absent vs schema present distinction,
3. `C43`
   - filename/path convention discipline,
4. `C44`
   - minimal template content discipline.

## Question

Is it honest to refine `provider_object` into one explicit carrier candidate
of the form:

```text
pair-indexed provider carrier candidate
```

with:

1. pair scope,
2. intended outputs,
3. filename/path convention,
4. minimal template content,
5. creation/emission status,

while still keeping:

1. no created carrier file,
2. no emitted provider instance,
3. no actual `theta_1`, `theta_2`.

## Hard limits

`T07` must not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual provider emission,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
