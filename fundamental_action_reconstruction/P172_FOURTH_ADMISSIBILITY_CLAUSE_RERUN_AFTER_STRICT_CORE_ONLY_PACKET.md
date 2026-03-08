# P172 Fourth Admissibility Clause Rerun After Strict-Core-Only Packet

Status: `P172_EXECUTED_FOURTH_ADMISSIBILITY_CLAUSE_RERUN_AFTER_STRICT_CORE_ONLY_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Rerun only the fourth admissibility clause:

```text
strict_core_only
```

for `S_preLM_strict_core_source_object_v1`.

## Positive rule

The rerun may return a positive result only if:

1. the first three clauses were already discharged,
2. the exported object is explicitly marked `strict_core_only = true`,
3. it is not laundered through the axiom lane,
4. observer remains excluded from the exported carrier.

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
