# P171 Third Admissibility Clause Rerun After Source-Seed-Only Packet

Status: `P171_EXECUTED_THIRD_ADMISSIBILITY_CLAUSE_RERUN_AFTER_SOURCE_SEED_ONLY_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Rerun only the third admissibility clause:

```text
source_seed_only
```

for `S_preLM_strict_core_source_object_v1`.

## Positive rule

The rerun may return a positive result only if:

1. the first two clauses were already discharged,
2. the object is explicitly marked `source_seed_only = true`,
3. `E_orient` is not exported,
4. `B_sel`, `R_sel`, and `O_sel` are not exported.

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
