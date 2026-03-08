# P174 Sixth Admissibility Clause Rerun After Selector-Acceptance-Independence Packet

Status: `P174_EXECUTED_SIXTH_ADMISSIBILITY_CLAUSE_RERUN_AFTER_SELECTOR_ACCEPTANCE_INDEPENDENCE_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Rerun only the sixth admissibility clause:

```text
selector_acceptance_independent
```

for `S_preLM_strict_core_source_object_v1`.

## Positive rule

The rerun may return a positive result only if:

1. the first five clauses were already discharged,
2. the exported object does not use the axiom lane,
3. it remains strict-core-only,
4. theory-level selector acceptance exists only as `axiom_augmented_only`,
5. that theory-level acceptance leaves strict core unchanged.

## Hard limits

Even a positive result still does **not** imply:

- strict-core selector closure,
- `QW-2191` discharge,
- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- ToE closure.
