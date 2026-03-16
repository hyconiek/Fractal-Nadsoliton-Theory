# P649 Current Second Admissibility Clause Rerun After Carrier Typing Probe (S_sel_int strict‑core source object v1)

Status: `P649_EXECUTABLE_CURRENT_SECOND_ADMISSIBILITY_CLAUSE_RERUN_AFTER_CARRIER_TYPING_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Rerun only the second admissibility clause:

```text
carrier_typed_enough_for_later_E_orient_export_required
```

for the currently exported strict‑core source object:

```text
S_sel_int_strict_core_source_object_v1
```

after the explicit carrier typing packet `F649`, and assuming the first clause
was already discharged by `N540`.

## Positive rule

The rerun may return a positive result only if:

1. the first clause was already discharged (`N540`),
2. the exported object has an explicit typed carrier interface (`F649`),
3. the future `E_orient` target frame is explicit (`F649`),
4. the designated sine axis has nonzero support (`F649` state support),
5. `E_orient` itself is **not** falsely claimed as already exported.

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

