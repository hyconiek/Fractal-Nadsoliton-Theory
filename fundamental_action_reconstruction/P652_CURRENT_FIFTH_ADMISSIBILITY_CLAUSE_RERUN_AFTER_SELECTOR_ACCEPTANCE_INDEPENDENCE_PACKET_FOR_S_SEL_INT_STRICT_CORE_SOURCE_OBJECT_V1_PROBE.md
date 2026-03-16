# P652 Current Fifth Admissibility Clause Rerun After Selector‑Acceptance‑Independence Packet (S_sel_int strict‑core source object v1)

Status: `P652_EXECUTABLE_CURRENT_FIFTH_ADMISSIBILITY_CLAUSE_RERUN_AFTER_SELECTOR_ACCEPTANCE_INDEPENDENCE_PACKET_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Rerun only the fifth admissibility clause:

```text
selector_acceptance_outside_strict_core_may_not_count_as_source_construction
```

for `S_sel_int_strict_core_source_object_v1`, after the explicit selector‑acceptance‑independence
packet `F652`, and assuming the first four clauses were already discharged
(`N540`, `N541`, `N542`, `N543`).

## Positive rule

The rerun may return a positive result only if:

1. the first four clauses were already discharged,
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

