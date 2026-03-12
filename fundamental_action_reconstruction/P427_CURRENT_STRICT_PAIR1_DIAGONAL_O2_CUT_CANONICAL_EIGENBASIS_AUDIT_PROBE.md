# P427 Current Strict `pair1` Diagonal `O(2)`-Cut Canonical Eigenbasis Audit Probe

Status: `P427_EXECUTED_CURRENT_STRICT_PAIR1_DIAGONAL_O2_CUT_CANONICAL_EIGENBASIS_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After:

1. the diagonal `pair1` `O(2)`-cut criterion (`N466`, audited on toy profiles in `P425`), and
2. the canonical `n=12` six-class reduction of the diagonal mode-2 defect (`N467/P426`),

the remaining strict missing step is not another slogan but one explicit reconstruction formula:

```text
if F2(d) ≠ 0, what is the canonical axis/eigenbasis on pair1 induced by the diagonal sector?
```

`P427` audits the explicit reconstruction formulas packaged as `N468` on concrete toy diagonal profiles
(constant, mode-1, mode-2 cosine/sine, mixed).

This probe does **not** claim any strict-derived diagonal coefficient instantiation for FIN.

## Inputs reused

1. `P425`
   - toy diagonal profiles on the `n=12` carrier.
2. `N466`
   - relation between the `pair1` block anisotropy and `F2(d)`.
3. `N468`
   - canonical eigenbasis reconstruction formula from `F2(d)` (conditional on `F2(d)≠0`).

## Result

`P427` exports a persisted artifact computing, for each toy diagonal profile:

1. the complex defect `F2(d)` and its magnitude `|F2(d)|`,
2. the canonical angle
   $$
   \theta_*=\frac12\,\operatorname{atan2}(\mathrm{Im}F_2,\mathrm{Re}F_2),
   $$
   (only defined when `F2(d)≠0`),
3. the rotated `pair1` block in the basis $(c_1',s_1')$ induced by $\theta_*$,
4. the eigenvalues predicted by `N468`:
   $$
   \lambda_\pm = \frac1n\sum_k d_k \pm \frac1n|F_2(d)|,
   $$
   and the absolute errors vs the direct eigenvalues of the `pair1` block.

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p427_current_strict_pair1_diagonal_o2_cut_canonical_eigenbasis_audit_probe.json`
- `fundamental_action_reconstruction/generated/p427_current_strict_pair1_diagonal_o2_cut_canonical_eigenbasis_audit_probe_summary.json`

## Honest verdict (no false pass)

The audit confirms the reconstruction formulas of `N468` on toy profiles:

1. when `F2(d)=0` (constant/mode-1 profiles), no canonical angle exists (the `pair1` block is scalar),
2. when `F2(d)≠0` (mode-2 and mixed profiles), the computed $\theta_*$ diagonalizes the `pair1` block and the
   eigenvalue formula matches within numerical tolerance.

This still does **not** produce a strict-core `QW-2191` discharge, because the canonical FIN diagonal coefficients
remain symbolic (`P426`) and no strict-derived instantiation deciding `F2(d)≠0` is exported.

