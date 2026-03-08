# N240 Current First Source Topology Observer-Free Scope Subtarget Theorem

Status: `N240_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Claim

Given the already exported packet `F132` and the current-repo verification
probe `P220`, the current repo state supports exactly the following scoped
theorem:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET
BELOW_ACTUAL_OBSERVER_FREE_SCOPE_DISCHARGE
```

and nothing stronger.

## Inputs

This theorem depends only on:

1. `F132_FIRST_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_PACKET`,
2. `P220_CURRENT_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_PROBE`.

## Theorem content

The theorem states:

1. the repo exports one explicit future-only subtarget
   `Omega_src_observer_free_scope_target_v1 : tau_src_candidate_v1 -> observer_free_scope_tag_v1`,
2. the exported subtarget remains source-side and observer-free in the witness
   domain,
3. the result remains strictly below:
   - actual observer-free scope discharge,
   - actual barrier-protected sign discharge,
   - actual nonzero-flow discharge,
   - full source-topology nontriviality,
   - selector promotion,
   - quotient-safe `QW-2191` resolution,
   - current selector closure.

## Why this does not overclaim

`N240` does not identify:

1. observer-free scope with actual source-topology nontriviality,
2. observer-free scope with selector promotion,
3. observer-free scope with a basis-independent witness,
4. observer-free scope with `QW-2191` discharge,
5. observer-free scope with current selector closure.

It is only the third scoped subtarget theorem below `Lambda_src_nontriv_target_v1`.

## Result

`N240` discharges only the following theorem-level statement:

```text
N240_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_THEOREM_NO_FALSE_PASS
```

## Hard limits

`N240` does not discharge:

1. actual observer-free scope of `tau_src_candidate_v1`,
2. actual barrier-protected sign,
3. actual nonzero-flow,
4. full source-topology nontriviality,
5. basis-independent selector promotion,
6. quotient-safe `QW-2191` resolution,
7. current selector closure,
8. current global `QW-2191` discharge,
9. ToE closure.
