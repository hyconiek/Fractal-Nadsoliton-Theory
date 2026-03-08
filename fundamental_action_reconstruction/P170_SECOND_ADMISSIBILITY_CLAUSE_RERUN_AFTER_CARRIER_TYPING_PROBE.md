# P170 Second Admissibility Clause Rerun After Carrier Typing Probe

Status: `P170_EXECUTED_SECOND_ADMISSIBILITY_CLAUSE_RERUN_AFTER_CARRIER_TYPING_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Rerun only the second admissibility clause:

```text
carrier_typed_enough_for_later_export
```

for `S_preLM_strict_core_source_object_v1`.

## Positive rule

The rerun may return a positive result only if:

1. the first clause was already discharged,
2. the exported object has explicit typed carrier decomposition,
3. topological, light, and matter slots are explicit,
4. the light slot has nonzero support,
5. a future `E_orient` target frame is explicit,
6. `E_orient` itself is not falsely claimed as already exported.

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
