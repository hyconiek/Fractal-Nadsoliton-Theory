# P409 Current Strict Sigma-Int → Theta Selector Ingredient (O(2)-Cut) Strict-Core Upgrade Incompatibility Boundary Probe

Status: `P409_EXECUTED_CURRENT_STRICT_SIGMA_INT_TO_THETA_SELECTOR_INGREDIENT_O2_CUT_STRICT_CORE_UPGRADE_INCOMPATIBILITY_BOUNDARY_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Probe whether the current repo may already treat the strict sigma-int selector-ingredient lane as having
crossed the strict-core upgrade boundary:

```text
candidate O(2)-cut (T157)  ->  strict-core canonical O(2)-cut (T159)
```

without introducing any hidden selector slot (`eps`, `delta_d`, …) and without implying any discharge of
`QW-2191`.

## Probe table

| Check | Verdict | Evidence / meaning |
|---|---|---|
| candidate selector ingredient exists (T157 level) | YES | `F325/N436/P400` export `ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1` (cuts O(2) only for the recorded `delta_d` and premise-based inputs) |
| strict-core (non-candidate) theta output exists | NO | `N1/C50` remain: no strict-core internal theta source; `F312/F314` are candidate-only records |
| `delta_d` slot eliminated or strict-derived | NO | corridor admits `delta_d ∈ (0,delta_max]` (`T119`); dependence is proved (`P403/N437`); `delta_d := delta_max` is premise-based (`F328/N440`) |
| `eps` slot eliminated or strict-derived | NO | generator admits `eps ∈ [0,1]` (`T117`); dependence is proved (`P407/N441`); `eps := 1/2` is premise-based (`F317/N428`) |
| “charge parity split ⇒ eps=1/2” exported as strict theorem | NO | `P408` audit: no strict exported law derives eps=1/2; current status is premise-based value object only |
| “information saturation ⇒ delta_max unique” exported as strict theorem | NO | `P408` audit: no strict exported objective/uniqueness theorem selects `delta_d = delta_max`; `delta_d` remains a selector slot (`N437`) |
| strict-core canonical O(2)-cut upgrade discharged (T159) | NO | `P406` discharge probe: fails due to exposed selector slots and lack of strict-core theta output |

## Probe result

The strongest honest current repo reading is:

```text
T157: discharged at candidate-ingredient level
T159: not discharged; strict-core upgrade is blocked by exposed selector slots eps and delta_d
```

So the strict sigma-int lane may be used only in candidate / premise-scoped form unless a genuinely new
strict-core selector source eliminating those slots is exported.

