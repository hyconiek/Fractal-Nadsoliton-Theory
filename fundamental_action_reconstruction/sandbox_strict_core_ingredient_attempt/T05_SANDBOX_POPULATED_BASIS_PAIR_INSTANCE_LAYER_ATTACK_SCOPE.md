# T05 Sandbox Populated Basis-Pair Instance Layer Attack Scope

Status: `T05_SANDBOX_POPULATED_BASIS_PAIR_INSTANCE_LAYER_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F04/P04/N04`, the sandbox already contains one conditional
strict-core theta-source rule candidate.

The next direct question is:

```text
can the sandbox attack the missing populated basis-pair instance layer
itself, rather than only the theta-source side feeding it?
```

This step does **not** try to produce an actual populated basis-pair instance.
It only tries to expose that missing layer as one explicit packet-ready
artifact schema.

## Support reused

1. `C48`
   - minimal actual basis-pair export skeleton,
2. `C49`
   - conditional populated-instance schema,
3. `C50`
   - strict-core minimal source skeleton absent,
4. `C40`
   - field-list present vs assembled artifact absent pattern,
5. `C41`
   - packet-ready artifact schema discipline,
6. `F04`
   - conditional theta-source rule candidate.

## Question

Is the following narrow move honest?

```text
write one packet-ready artifact schema for the missing populated basis-pair
instance layer, while preserving:
  - no actual theta_1/theta_2 supply,
  - no actual populated u_1/u_2 instance,
  - no actual export.
```

## Hard limits

`T05` must not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.

It must also avoid importing the axiom-augmented branch.
