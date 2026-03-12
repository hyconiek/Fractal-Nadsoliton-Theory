# P411 Current Strict Sigma-Int Positive-Window Delta_d Step Strict-Derived Selection Target Probe

Status: `P411_EXECUTED_CURRENT_STRICT_SIGMA_INT_POSITIVE_WINDOW_DELTA_D_STEP_STRICT_DERIVED_SELECTION_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Probe whether the current repo already exports any **actual** strict-derived
(not premise-only) delta_d-selection ingredient discharging `T161`, or whether
the strongest honest current state remains:

```text
delta_d is fixed only by explicit strict-side premise (strict provenance),
and delta_d remains a real selector slot on the strict sigma-int → theta lane.
```

## Probe table

| Check | Verdict | Evidence |
|---|---|---|
| corridor admits free `delta_d ∈ (0, delta_max]` | YES | `T119` |
| theta-pair depends on admissible `delta_d` choice | YES | `P403/N437` |
| dedicated strict-provenance delta_d value object exported | YES | `F328/N440` export `delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max` as `strict_source_upgraded` (premise) |
| dedicated strict-derived delta_d value object exported | NO | no export of `delta_d_sigma_int_positive_window_step_strict_derived_v1`; `P408` audit: no strict-derived uniqueness/objective theorem selects `delta_d = delta_max` |

## Exact verdict

The strongest honest current verdict is:

```text
strict-derived delta_d selection ingredient (T161): NOT EXPORTED
delta_d remains a real selector slot on the strict sigma-int → theta lane (N437).
```

Therefore any strict-core upgrade attempt requiring delta_d-slot closure
(e.g. `T159`) must remain explicitly conditional on either:

1. an explicit delta_d selector premise (current state), or
2. a future strict-derived delta_d selection law (not yet exported).

No strict-core theta export, strict-core selector closure, `QW-2191` discharge,
or ToE closure is implied.

