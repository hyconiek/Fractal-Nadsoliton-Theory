# P651 Current Fourth Admissibility Clause Rerun After Strict‑Core‑Only Packet (S_sel_int strict‑core source object v1)

Status: `P651_EXECUTABLE_CURRENT_FOURTH_ADMISSIBILITY_CLAUSE_RERUN_AFTER_STRICT_CORE_ONLY_PACKET_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Rerun only the fourth admissibility clause:

```text
strict_core_only_required
```

for `S_sel_int_strict_core_source_object_v1`, after the explicit strict‑core‑only
packet `F651`, and assuming the first three clauses were already discharged
(`N540`, `N541`, `N542`).

## Positive rule

The rerun may return a positive result only if:

1. the first three clauses were already discharged,
2. the exported object is explicitly marked `strict_core_only = true` (`F651`),
3. it is not laundered through the axiom lane (`F651`),
4. it remains upstream of the observer (`F651`).

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

