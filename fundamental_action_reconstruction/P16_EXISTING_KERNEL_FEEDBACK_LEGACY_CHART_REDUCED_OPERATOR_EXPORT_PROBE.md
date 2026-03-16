# P16 Existing Kernel Feedback Legacy Chart-Reduced Operator Export Probe

Status: `P16_EXECUTED_EXISTING_KERNEL_FEEDBACK_LEGACY_CHART_REDUCED_OPERATOR_EXPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-16`

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
2. therefore the declared-control artifacts (canonical `H_PsiPsi` and declared
   `M_control`) remain below a strict existing-feedback promotion.

Update (2026-03-16): in declared scope the repo now **does** export:

- an explicit coefficient-filled canonical Psi x Psi block `H_PsiPsi` on full
  declared transport support (`R12`),
- an explicit coefficient-filled declared control pullback
  `M_control = T_control^T H_PsiPsi T_control` (`P476`).

This removes two of the three original `P16` sub-blockers, but it still does
not supply a host-to-canonical matching witness identifying `QW-2186` with the
canonical carrier (`C10_B1`).

## Sharpened decomposition of the first `P15` blocker

`P16` now reduces the first `P15` blocker to one remaining missing object:

1. host-to-concrete Psi-sector quadratic block identification for the
   existing-kernel-feedback host operator (matching `QW-2186` to the canonical
   carrier).

## Honest frontier

`P16` shows that the route fails before any legacy-side matrix on `pair1`
exists. The obstruction is now localized upstream of `pair1` to:

1. host-to-Psi block identification.

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

1. export a host-to-canonical Psi-block matching witness identifying `QW-2186`
   with the canonical carrier (or equivalently with its declared control
   pullback),
2. or keep the legacy operator export route negative and do not claim any
   coefficient-filled legacy-side matrix on `pair1`.
