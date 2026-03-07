# CF1 Compute-Or-Fail a1 b1 d1

## Goal

Run one executable test that either:

- computes `a_1, b_1, d_1` from the current repository state, or
- returns a machine-level failure explaining why this is not computable today.

## Sources Used

- `fundamental_action_reconstruction/generated/o2_exported_composite_a1_ext_instance.json`
- `fundamental_action_reconstruction/generated/o3_a1_ext_coefficient_evaluation_rule.json`
- `material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2165_l13_exhaustive_canonical_eom_gate.json`

## Executable Criterion

The run is considered successful only if the persisted `A_1_ext(pair1)` instance
contains actual numeric or directly evaluable entries for:

- `a_1`
- `b_1`
- `d_1`

Otherwise the run must terminate with:

- `NOT_COMPUTABLE_FROM_CURRENT_REPO_STATE`

## Honest Result

At the current repository state:

- `O2` provides a persisted `A_1_ext(pair1)` object,
- `O3` provides the readout rule,
- but the entries of `A_1_ext(pair1)` remain symbolic placeholders,
- and `QW-2165` exports symbolic canonical EoM plus kernel-mixing terms, not an actual selector-sector `2x2` block.

So the executable outcome is:

- `NOT_COMPUTABLE_FROM_CURRENT_REPO_STATE`

## Meaning

This is not a theorem-level closure result.
It is an execution-level result saying that, today, the repository does not yet
contain enough operator-level populated data to produce `a_1, b_1, d_1`.
