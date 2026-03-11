# P405 Current Strict Sigma-Int Positive-Window Delta_d Step Strict-Provenance Target Probe

Status: `P405_EXECUTED_CURRENT_STRICT_SIGMA_INT_POSITIVE_WINDOW_DELTA_D_STEP_STRICT_PROVENANCE_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo already exports an **actual** dedicated delta_d
step value object with explicit strict provenance/selection method, or whether
the strongest honest state remains only:

```text
future-only target naming for the missing strict-provenance delta_d step ingredient (T158/F327).
```

This probe exists to prevent false pass by silently treating a corridor-step
choice inside a candidate instantiation artifact as if it were strict-derived.

## Probe table

| Check | Verdict | Evidence |
|---|---|---|
| corridor admits free `delta_d ∈ (0, delta_max]` | YES | `T119` |
| theta-pair depends on admissible `delta_d` choice | YES | `P403/N437` |
| extension-scope convention fixing `delta_d := delta_max` exists | YES (extension only) | `AX17` (`strict_extension_only`) |
| dedicated strict-provenance delta_d value object exported | NO | no exported artifact `generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json` (or equivalent) exists on current repo state |
| delta_d choice recorded inside instantiation artifacts | YES (embedded only) | `F314` / `F325` artifacts record `delta_d` as a chosen corridor step, but do not export it as a dedicated value object |

## Exact verdict

The strongest honest current verdict is:

```text
strict-provenance dedicated delta_d step value object: NOT EXPORTED
delta_d selection remains either:
  (a) an embedded explicit choice inside candidate artifacts, or
  (b) an extension-scope convention (AX17),
therefore T158 remains a live missing-ingredient target for strict-core provenance hygiene.
```

No strict-core theta export, object-support discharge, selector closure, `QW-2191`
discharge, or ToE closure is implied.

