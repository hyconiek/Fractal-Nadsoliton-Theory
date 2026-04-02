# T306 Current Strict `T173/T176` Source-Side Input-Leg Attempt Exact Lower Supplier-Boundary Target AR Attempt Exact Further Lower Boundary Target AR Attempt Spec

Status: `T306_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_SPEC_NO_FALSE_PASS`
As of: `2026-03-23`

## Goal

After `P1057/N890`, the exact `T305` target still remains future-only on the
current repo state.

The honest next question is therefore:

```text
what exact first actual-realization attempt instance should now be treated as
the active primary move on that exact T305 target?
```

## First actual-realization attempt instance

Freeze one exact first attempt instance:

```text
W_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_attempt_v1
```

with the following scoped attempt shape:

```text
exact_further_lower_boundary_target_actual_realization_attempt_v1(
  source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_v1,
  W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1,
  Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1
)
```

## Frozen Attempt Discipline

- `attempt_is_over_exact_t305_target := yes`
- `attempt_uses_t224_route_local_actual_attempt_as_best_current_support_not_as_already_lawful_supplier := yes`
- `attempt_must_not_reenter_exhausted_pair12_entry_point_same_lane_as_primary_descent := yes`
- `attempt_keeps_noncyclic_provider_shift_anchor_explicit := yes`
- `attempt_must_remain_below_actual_source_side_input_leg_export := yes`
- `attempt_must_remain_below_actual_bridge_output_schema_export := yes`
- `attempt_must_remain_below_actual_full_C_v1_transported_section_lift := yes`
- `attempt_must_remain_below_actual_T176_discharge := yes`

## Hard Limits

`T306` does **not** claim:

1. actual export of the exact further lower-boundary,
2. actual export of the exact lower supplier-boundary,
3. actual export of the exact `source_side_input_leg`,
4. actual bridge-output-schema export,
5. actual full `C_v1` transported-section lift,
6. actual `T176` discharge,
7. `QW-2191` discharge,
8. ToE closure.
