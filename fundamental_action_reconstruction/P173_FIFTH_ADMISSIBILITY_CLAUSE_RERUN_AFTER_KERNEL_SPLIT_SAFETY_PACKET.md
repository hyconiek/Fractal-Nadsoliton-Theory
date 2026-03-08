# P173 Fifth Admissibility Clause Rerun After Kernel-Split-Safety Packet

Status: `P173_EXECUTED_FIFTH_ADMISSIBILITY_CLAUSE_RERUN_AFTER_KERNEL_SPLIT_SAFETY_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Rerun only the fifth admissibility clause:

```text
non_substitutive_wrt_kernel_split
```

for `S_preLM_strict_core_source_object_v1`.

## Positive rule

The rerun may return a positive result only if:

1. the first four clauses were already discharged,
2. the exported object is explicitly `kernel_split_safe = true`,
3. it does not import an external selector,
4. the admissibility-upgrade guardrail remains `kernel_split_safe = true`.

## Hard limits

Even a positive result still does **not** imply:

- a legacy-to-strict kernel bridge,
- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
