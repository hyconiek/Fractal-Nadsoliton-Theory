# P16 Existing Kernel Feedback Legacy Chart-Reduced Operator Export Probe

Status: `P16_EXECUTED_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P15`, the first remaining blocker became:

```text
explicit_coefficient_filled_legacy_chart_reduced_operator_object_on_the_chosen_current_pair_chart_pair1_or_equivalent_actual_target
```

`P16` tests that object directly in `compute-or-fail` form.

The acceptance gate is:

- either the repo now exports a coefficient-filled legacy chart-reduced
  operator on `pair1`,
- or the missing structure is sharpened into a finite upstream blocker-set.

## Result

The route still does **not** compute.

The report returns:

```text
NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_ROUTE
```

## What is present but still insufficient

The repo now contains all of the following:

1. a host-scope existing-kernel-feedback carrier from `R8`,
2. a typed host-to-control pushforward from `R9`,
3. a formal control pullback formula
   `M_control = T_control^T H_PsiPsi T_control` from `C15`,
4. a chosen current-pair chart reduction `Pi_pair1 : M_control -> V_1` from
   `R10`.

But this still does **not** amount to a coefficient-filled legacy operator on
`pair1`, because:

1. the legacy host is still not identified with a concrete Psi-sector block,
2. no executed coefficient-filled Psi-sector block export is present,
3. no coefficient-filled control pullback `M_control` and therefore no
   coefficient-filled chart-reduced `pair1` block is exported.

## Sharpened decomposition of the first `P15` blocker

`P16` reduces the first `P15` blocker to three current missing objects:

1. host-to-concrete Psi-sector quadratic block identification for the
   existing-kernel-feedback host operator,
2. an executed and persisted coefficient-filled Psi-sector block export
   supporting `H_PsiPsi` evaluation,
3. a coefficient-filled control pullback `M_control` and its chart-reduced
   `pair1` block export.

## Honest frontier

`P16` shows that the route fails before any legacy-side matrix on `pair1`
exists. The obstruction is now localized upstream of `pair1` to:

1. host-to-Psi block identification,
2. coefficient-filled Psi-block export,
3. coefficient-filled control pullback and chart-reduced export.

## What `P16` does not claim

`P16` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `M_control` is already coefficient-filled,
- that a concrete legacy `2x2` matrix on `pair1` exists,
- that existing kernel feedback has been identified with the computed `P10`
  block,
- that `QW-2191` is discharged,
- that ToE is closed.

## Recommended next move

Only two honest routes remain:

1. export a concrete existing-feedback Psi-sector block and evaluate the
   control pullback to `pair1`,
2. or keep the legacy operator export route negative and do not claim any
   coefficient-filled legacy-side matrix on `pair1`.
