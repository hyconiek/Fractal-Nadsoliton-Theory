# P650 Current Third Admissibility Clause Rerun After Source‑Seed‑Only Packet (S_sel_int strict‑core source object v1)

Status: `P650_EXECUTABLE_CURRENT_THIRD_ADMISSIBILITY_CLAUSE_RERUN_AFTER_SOURCE_SEED_ONLY_PACKET_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Rerun only the third admissibility clause:

```text
source_seed_only_not_counted_as_E_orient_or_bridge
```

for `S_sel_int_strict_core_source_object_v1`, after the explicit source‑seed‑only
packet `F650`, and assuming the first two clauses were already discharged
(`N540`, `N541`).

## Positive rule

The rerun may return a positive result only if:

1. the first two clauses were already discharged,
2. the object is explicitly marked `source_seed_only = true` (`F650`),
3. `E_orient` is not exported (`F650`),
4. `B_sel`, `R_sel`, and `O_sel` are not exported (`F650`).

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

