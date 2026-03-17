# P672 Current Strict T172 Readiness Probe (global selector closure + QW‑2191 discipline)

Status: `P672_EXECUTABLE_CURRENT_STRICT_T172_READINESS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit whether the **current repo state** is past the “global infrastructure export” phase and is now honestly at the `T172`
frontier:

```text
global strict selector closure + the corresponding QW-2191 uniqueness discipline
```

This probe is **readiness-only**: it does **not** export any closure object and it does **not** claim any `QW-2191` discharge.

## Positive rule

The probe may return “prereqs present” only if:

1. `T170` is discharged (`F469`): global selector atlas + transition objects on `C_v1`,
2. global selector state objects exist on `C_v1`:
   - projective (`F470`),
   - directed in an explicit premise-tracked scope (`F473` + `F474`, via `T164/T171`),
3. the seed‑v1 strict-core internal selector-source lane has been promoted globally on `C_v1` (`N550–N553`),
4. those global promotions remain explicitly **projector/section-level only** (respect `N512` boundary).

## What the probe intentionally does not “PASS”

Even if prereqs are present, this probe does **not** claim:

- a global closure object exists (`SelectorClosure_global_C_v1_*`),
- strict-core selector closure,
- global `QW-2191` discharge,
- ToE closure.

