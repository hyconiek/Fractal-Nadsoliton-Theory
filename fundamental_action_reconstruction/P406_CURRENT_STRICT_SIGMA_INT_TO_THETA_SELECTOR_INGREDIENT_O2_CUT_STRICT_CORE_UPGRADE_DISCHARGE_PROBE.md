# P406 Current Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Strict-Core Upgrade Discharge Probe

Status: `P406_EXECUTED_CURRENT_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_STRICT_CORE_UPGRADE_DISCHARGE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

`P400` confirmed discharge of `T157` at the **candidate selector-ingredient** level.

`T159` names the next missing upgrade target: a strict-core selector ingredient that canonically
cuts the `QW-2191` `O(2)` family without leaving hidden selector slots (`delta_d`, `eps`, …).

`P406` probes whether the current repo already discharges that stricter `T159` target.

## Probe table (T159 acceptance tests)

| Acceptance test (T159) | Verdict | Evidence |
|---|---|---|
| strict-core (non-candidate) theta output `(theta_1,theta_2)` exported | NO | strict lane exports only theta **candidate** records (`F312/N423`, `F314/N425`) and one **candidate** selector ingredient (`F325/N436/P400`); `N1/C50` remain: no strict-core theta source |
| hidden selector slot `delta_d` eliminated or strict-derived | NO | `P403/N437` record theta dependence on admissible `delta_d`; `F328/N440` exports `delta_d := delta_max` only as `strict_source_upgraded` (premise), not strict-derived |
| hidden selector slot `eps` eliminated or strict-derived | NO | `F317/N428` exports `eps := 1/2` only as `strict_source_upgraded` (premise), not strict-derived; `T117` still treats eps as a free generator parameter |
| `O(2)` family canonically cut in strict core (no “chosen representative” caveat) | NO | `F325/N436` cut the `O(2)` family only after the explicit corridor-step choice; `QW-2191` still blocks kernel-alone uniqueness and no strict-core internal selector source has been exported beyond candidate-ingredient level |
| noncyclic and observer-free contracts satisfied | YES (for the candidate layer only) | `F325/N436` are noncyclic and observer-free, but they remain candidate-ingredient level and do not satisfy the strict-core upgrade requirements of `T159` |

## Exact verdict

On the current repo state:

```text
T157: discharged (candidate selector ingredient)  (F325/N436/P400)
T159: NOT discharged (strict-core upgrade target remains open)
QW-2191: remains open in strict core
```

## Next honest step

The next honest move is not rhetorical “time-arrow selector” promotion (`P399` already audits that).

It must be one genuinely new strict-core ingredient that satisfies `T159`:

1. eliminate the exposed selector slots (`delta_d`, `eps`) by invariance or strict derivation, and
2. upgrade the sigma-int → theta selector from candidate-ingredient into a strict-core canonicalization
   source capable of resolving the `QW-2191` `O(2)` family in declared scope,

or proceed explicitly on a separated non-strict scope (axiom-augmented / strict-extension) without
claiming strict-core internalization.

