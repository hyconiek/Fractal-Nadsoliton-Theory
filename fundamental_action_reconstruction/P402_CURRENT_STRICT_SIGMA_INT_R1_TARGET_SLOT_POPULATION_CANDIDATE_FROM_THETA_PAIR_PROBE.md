# P402 Current Strict Sigma-Int R1 Target-Slot Population Candidate (from Theta-Pair) Probe

Status: `P402_EXECUTED_CURRENT_STRICT_SIGMA_INT_R1_TARGET_SLOT_POPULATION_CANDIDATE_FROM_THETA_PAIR_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `R1/F311`, the strict sigma-int residual lane exports:

1. the residual-orientation target slot `R1` (unpopulated; requires `theta_1,theta_2`), and
2. a strict-core sigma-int → target-slot export-map object (`F311/N422`),
   but **sign-only** (no theta supply; no population).

After `F325/N436` and the `P401` selected-point witness, the repo also exports
one strict-side **candidate** theta-pair selector ingredient:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1
```

persisted as:

```text
fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json
```

This probe checks one narrow question:

```text
can we construct (and audit) a candidate inhabitant instance of the R1 target slot
S_orient_cand(theta_1,theta_2)=span{u_1(theta_1),u_2(theta_2)}
using the exported theta-pair selector ingredient,
without claiming strict-core theta export, object-support, selector closure, or QW-2191 discharge?
```

## Implementation artifact

Executed by:

```text
python3 fundamental_action_reconstruction/p402_current_strict_sigma_int_r1_target_slot_population_candidate_from_theta_pair_probe.py
```

Outputs:

```text
fundamental_action_reconstruction/generated/r1_residual_orientation_datum_target_slot_population_candidate_from_sigma_int_theta_pair_v1.json
fundamental_action_reconstruction/generated/p402_current_strict_sigma_int_r1_target_slot_population_candidate_from_theta_pair_probe.json
fundamental_action_reconstruction/generated/p402_current_strict_sigma_int_r1_target_slot_population_candidate_from_theta_pair_probe_summary.json
```

## Exact verdict

On the current repo state:

```text
P402: PASS (candidate inhabitant instance constructed and audited).
```

Meaning:

1. the R1 target-slot object class admits a concrete candidate inhabitant instance
   formed from `(theta_1^cand,theta_2^cand)` and the induced `(u_1^cand,u_2^cand)` from `F325`,
2. basic internal audits (orthonormality, non-collinearity, scaffold consistency) pass,
3. no strict-core promotion is implied:
   this does not upgrade `F311` into a populated strict-core export-map object and does not discharge `N302` nor `QW-2191`.

## What P402 does not prove

`P402` does not prove:

1. actual strict-core theta export (`theta_1/theta_2` remain candidate-only),
2. actual strict-core population of the R1 target slot by a strict-core export-map object,
3. actual bridge/export-map object-support above the map object (`N302/N395` remain open),
4. admissible `S_sel_int` / strict-core selector closure,
5. `QW-2191` discharge,
6. ToE closure.

No false pass: this is a candidate inhabitant construction + audit only.

