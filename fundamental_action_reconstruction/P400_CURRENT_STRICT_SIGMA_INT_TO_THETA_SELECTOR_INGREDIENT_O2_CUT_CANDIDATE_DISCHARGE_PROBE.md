# P400 Current Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Candidate Discharge Probe

Status: `P400_EXECUTED_CURRENT_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_CANDIDATE_DISCHARGE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo discharges the `T157` target:

```text
export a strict-side sigma_int_strict_derived_v1 -> (theta_1^cand,theta_2^cand)
candidate selector ingredient cutting the QW-2191 O(2) family (declared scope),
without importing QW-2192/2193 as strict evidence.
```

## Probe table (T157 acceptance tests)

| Acceptance test (T157) | Verdict | Evidence |
|---|---|---|
| typed output `(theta_1^cand, theta_2^cand)` exported | YES | `F325/N436` export `ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1` persisted as `generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json` (numeric theta pair consistent with the existing `F314` positive-window instantiation record) |
| `atan2` well-defined (no degeneracy) | YES | positive-window corridor recorded in the artifact (`T119`-style) with explicit numeric `delta_d` and `X_i>0,Y_i>0` |
| strict provenance explicit (no silent conventions) | YES | artifact cites `N418/N428/N435` and records premise-based provenance explicitly |
| O(2)-cut argument included | YES | artifact includes explicit `o2_cut_argument` referencing `QW-2190/QW-2191` and the induced `u_1,u_2` |
| no false pass (no kernel-alone uniqueness / no closure smuggling) | YES | `F325/N436` keep `QW-2191` as obstruction and forbid selector closure / ToE closure claims |

## Exact verdict

On the current repo state:

```text
T157: DISCHARGED at candidate-ingredient level (F325/N436),
with explicit premise-based provenance and an explicit O(2)-cut witness
in the declared QW-2190/QW-2191 scope.
```
