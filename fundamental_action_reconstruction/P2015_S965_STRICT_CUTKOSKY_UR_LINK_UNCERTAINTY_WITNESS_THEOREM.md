# P2015 S965 Strict Cutkosky UR-Link Uncertainty Witness Theorem

Status: `P2015_EXECUTED_UR_LINK_UNCERTAINTY_BOUNDED_SCAN_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

Export a machine-checkable strict-lane witness that extends the Cutkosky proxy
chain with a bounded uncertainty table and residue-positivity stability scan for
`graviton -> gauge_gauge`, without claiming full all-state theorem closure.

## Construction

Using `P1997` grid outputs `(ImM, CutSum, Delta_opt)` on the strict lane,
`P2015` applies an admissible local perturbation family:

- multiplicative dressing uncertainty `epsilon_z1 ∈ [-0.03, 0.03]`,
- channel uncertainty `epsilon_ch ∈ [-0.05, 0.05]`.

For each energy sample `s` from the exported grid, it computes:

1. bounded interval of `Delta_opt` under the scan,
2. center/std of that uncertainty band,
3. positivity flags for residue proxy transport on strict-safe poles.

## Result

`P2015` exports:

- bounded uncertainty table,
- positivity stability flags under the local scan,
- gatekeeper checks and explicit no-false-pass guard.

This is a strict-lane UR-link **diagnostic strengthening**, not a theorem that
full unitarity closure is achieved.

## Next honest step

Replace the multiplicative proxy perturbation with explicit loop-derived
channel amplitudes and re-run uncertainty transport on a wider RG-safe atlas.
