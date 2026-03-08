# P169 First Admissibility Clause Rerun After Strict-Core Export Probe

Status: `P169_EXECUTED_FIRST_ADMISSIBILITY_CLAUSE_RERUN_AFTER_STRICT_CORE_EXPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Rerun only the first admissibility clause after the explicit export:

```text
does S_preLM_strict_core_source_object_v1 now satisfy
genuinely_new_strict_core_source_object_required?
```

## Positive rule

The rerun may return a positive result only if all of the following hold:

1. a new strict-core exported object identity is present,
2. that object is marked `constructed_source_object_exported = true`,
3. it is strict-core-only and upstream of observer,
4. it carries no external selector import,
5. it does not reuse the blocked artifact families
   (`psi0`, FR, `sigma_int` as source object, overlay-fit, axiom-lane object),
6. a nonzero nonreduction witness against the same-basis `F75` packaging is
   present.

## Hard limits

Even a positive result still does **not** imply:

- full admissibility of `S_sel_int`,
- admissible `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
