# T06 Sandbox Strict-Core Theta/Input Provider Attack Scope

Status: `T06_SANDBOX_STRICT_CORE_THETA_INPUT_PROVIDER_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F05/P05/N05`, the sandbox already has one explicit downstream contract
for the missing populated basis-pair instance layer.

The next direct question is:

```text
can the sandbox write one explicit strict-core-only artifact schema for the
missing upstream theta/input provider that would have to feed F05?
```

This step does **not** try to produce actual `theta_1`, `theta_2`.
It only tries to make the missing provider layer explicit as one packet-ready
artifact schema.

## Support reused

1. `C35`
   - strict-core actual phase source absent; axiom-augmented branch present,
2. `C50`
   - no packet-ready strict-core minimal source skeleton,
3. `C52`
   - field-list present vs assembled bridge artifact absent pattern,
4. `C53`
   - artifact schema discipline for source-side bridge-like frontiers,
5. `F05`
   - downstream populated basis-pair instance layer artifact schema.

## Question

Is the following narrow move honest?

```text
write one packet-ready artifact schema for the missing strict-core theta/input
provider, with F05 as its downstream consumer contract,
while preserving:
  - no actual theta_1/theta_2 export,
  - no actual provider instance,
  - no identification with the axiom-augmented fallback lane
```

## Hard limits

`T06` must not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.

It must also avoid treating the axiom-augmented branch as a strict-core
provider.
