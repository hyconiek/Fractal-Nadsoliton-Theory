# T14 Sandbox Strict-Core Theta-Source Supply Boundary Direct Attack Scope

Status: `T14_SANDBOX_STRICT_CORE_THETA_SOURCE_SUPPLY_BOUNDARY_DIRECT_ATTACK_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F13/P13/N13`, the theta-output clause has positive support below actual
values, but the one remaining negative point inside that clause is exact:

```text
strict-core theta-source supply boundary
```

Current repo support says:

1. `C50`
   - strict-core minimal source skeleton is absent,
2. `C51`
   - fallback axiom-augmented lane is present, but no strict-to-axiom
     bridge-spec packet exists,
3. `C52`
   - the minimal bridge field list is already present, but no assembled bridge
     artifact exists.

## Direct question

Can the sandbox attack this boundary directly by assembling one explicit
strict-to-axiom theta-source bridge artifact schema while still refusing to
claim:

1. actual strict-core theta supply,
2. actual bridge discharge,
3. actual theta-output emission?

## Intended move

`T14` does not try to turn:

```text
actual_strict_core_theta_source_supply_present := no
```

into `yes`.

It attacks the boundary directly in a narrower way:

1. cite the exact strict blocker,
2. cite the exact non-strict fallback lane,
3. cite the current bridge class,
4. assemble one explicit schema artifact joining those fields,
5. keep the whole result below source supply and below emission.

## Hard limits

`T14` must not claim:

1. fallback lane = strict-core source,
2. control-route overlay = strict bridge discharge,
3. actual `theta_1`, `theta_2`,
4. actual provider emission,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
