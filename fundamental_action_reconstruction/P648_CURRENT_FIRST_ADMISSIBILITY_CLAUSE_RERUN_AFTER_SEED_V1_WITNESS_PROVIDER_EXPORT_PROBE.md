# P648 Current First Admissibility Clause Rerun Probe (Seed‑v1 witness provider)

Status: `P648_EXECUTABLE_CURRENT_FIRST_ADMISSIBILITY_CLAUSE_RERUN_AFTER_SEED_V1_WITNESS_PROVIDER_EXPORT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `F647/N539`, the next honest admissibility move is to rerun the **first**
clause of the frozen minimal admissibility contract `F34` for `S_sel_int`:

```text
genuinely_new_strict_core_source_object_required
```

This probe checks whether the newly exported strict-core constructed source
object satisfies the first clause **in the narrow export/property sense**,
without implying any later-clause discharge.

## Hard limits

`P648` does not claim admissible `S_sel_int`, admissible `E_orient`, selector
closure, `QW‑2191` discharge, or ToE closure.

