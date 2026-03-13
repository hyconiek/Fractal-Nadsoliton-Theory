# P432 Current Strict Canonical Local-Diagonal Mode‑2 Defect Value-Instantiation Repo-Scan Probe

Status: `P432_EXECUTED_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_VALUE_INSTANTIATION_REPO_SCAN_PROBE_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

After `N466/N467/P426`, the strict diagonal/local accelerator question on `pair1` is:

```text
canonical FIN D_local_residual cuts O(2) on pair1  <=>  F2(d) ≠ 0.
```

`N472/P431` already close one ambiguity:

```text
at the level of the exported canonical coefficient class (R15), F2(d) is underdetermined.
```

`P432` performs one additional hygiene check demanded by “no false pass” discipline and by the user workflow:

```text
does the repo already contain any exported numeric/value instantiation
of the canonical diagonal/local residual profile (or its coefficient symbols)
that would decide F2(d) anyway?
```

This probe is **not** a physics admissibility claim.
It is a repo-state scan for any already-exported value assignments sufficient to evaluate the `P426` expression.

## Strict-admissible evidence reused

1. `T166`
   - the current strict decision target: strict-derived decision of `F2(d)` for canonical FIN `D_local_residual`.
2. `P426`
   - the explicit six-class expression for `Re(F2)` and `Im(F2)` in terms of
     `Sigma_psi0_psi6, ..., Sigma_psi5_psi11`.
3. `N472/P431`
   - theorem+probe stating underdetermination at the level of the exported coefficient class.

## Scan method

`P432` is executed by the repo-local script:

- `fundamental_action_reconstruction/p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe.py`

The scan searches the repo (as text) for **numeric value assignments** to any of the following families:

1. direct diagonal/profile values (e.g. `d_k` profiles or opposite-pair sums),
2. the six opposite-pair sum symbols:
   `Sigma_psi0_psi6, ..., Sigma_psi5_psi11`,
3. the canonical diagonal coefficient symbols inside `R15`:
   `m2_psi*`, `g4_psi*`, `g6_psi*`, `gY*`, `vpsi*`, `vphi`.

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe.json`
- `fundamental_action_reconstruction/generated/p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe_summary.json`

## Result (no false pass)

On the current repo state, the scan finds **no decision‑ready exported numeric instantiation** of:

1. the canonical diagonal/local residual profile `d=(d_0,...,d_{11})`,
2. the six opposite-pair sums `Sigma_psi0_psi6, ..., Sigma_psi5_psi11`,
3. the canonical diagonal coefficient symbols sufficient to evaluate `P426`.

Therefore the strict diagonal accelerator decision target remains open:

```text
T166 is not decidable from already-exported numeric instantiations on current repo state.
```

The scan does detect *some* numeric assignments in toy/probe artifacts (e.g. witness specializations such as `vphi=0`)
but those are not an exported canonical diagonal profile instantiation and do not decide `F2(d)` for the canonical FIN
`D_local_residual`.

So the practical verdict remains consistent with `N472/P431`:
the diagonal/local lane cannot be promoted into strict core until additional strict-derived structure exports actual
values/constraints deciding `F2(d)` for the canonical FIN `D_local_residual`.
