# P15 Existing Kernel Feedback Intertwiner Equality Witness Probe

Status: `P15_EXECUTED_EXISTING_KERNEL_FEEDBACK_INTERTWINER_EQUALITY_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R10`, `P14` reduced the factorization route to a single nominal blocker:

```text
intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_object_with_the_computed_current_pair_H3_block
```

`P15` tests that remaining blocker directly in `compute-or-fail` form.

The acceptance gate is:

- either the repo already exports the required witness,
- or the missing structure is sharpened again into a smaller machine-readable
  blocker-set.

## Result

The route still does **not** compute.

The report returns:

```text
NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_INTERTWINER_EQUALITY_WITNESS_ROUTE
```

## What is present but still insufficient

The repo now contains all of the following:

1. a chosen current-pair chart reduction `Pi_pair1 : M_control -> V_1`,
2. a computed current-pair `H3` block from `P10`,
3. a provenance-valid extension-lane `Route A` witness from `H18`,
4. an exported composite `A_1_ext` instance from `O2`.

But this still does **not** amount to an intertwiner/equality witness, because:

1. the legacy side still does not export a coefficient-filled chart-reduced
   operator object on that chosen chart,
2. the strongest extension-lane composite witness remains coefficient-unresolved,
3. therefore no operator-level equality or intertwiner can be honestly written
   between legacy feedback and the computed `P10` block.

## Sharpened decomposition of the last `P14` blocker

`P15` reduces the single nominal `P14` blocker to two current missing objects:

1. an explicit coefficient-filled legacy chart-reduced operator object on the
   chosen current-pair chart `pair1` or an equivalent actual target,
2. an intertwiner/equality witness identifying that legacy object with the
   computed current-pair `H3` block.

## Honest frontier

`P15` shows that the current route still fails before operator identification.
The remaining obstruction is no longer phrased as one vague witness token, but
as:

1. missing legacy chart-reduced operator export,
2. missing equality/intertwiner witness after that export.

## What `P15` does not claim

`P15` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that the legacy route already exports a matrix on `pair1`,
- that the extension-lane composite witness is evaluated,
- that existing kernel feedback is identified with the computed `H3` block,
- that `QW-2191` is discharged,
- that ToE is closed.

## Recommended next move

Only two honest routes remain:

1. export the coefficient-filled legacy chart-reduced operator object on
   `pair1` first,
2. or keep the factorization route negative and do not claim any
   intertwiner/equality witness.
