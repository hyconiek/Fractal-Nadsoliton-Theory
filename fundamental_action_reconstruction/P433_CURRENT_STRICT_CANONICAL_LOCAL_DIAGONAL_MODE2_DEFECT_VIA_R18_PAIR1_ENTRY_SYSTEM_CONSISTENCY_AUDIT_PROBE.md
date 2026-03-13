# P433 Current Strict Canonical Local-Diagonal Mode‑2 Defect via `R18` `pair1` Entry System Consistency Audit Probe

Status: `P433_EXECUTED_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_VIA_R18_PAIR1_ENTRY_SYSTEM_CONSISTENCY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

`N473` proves a structural identity:

```text
Re(F2) = 6*(c1c1 - s1s1),   Im(F2) = 12*(c1s1),
and substituting the R18 entry formulas reproduces the N467/P426 six-class defect reduction.
```

`P433` audits this identity numerically on random assignments of the six opposite‑pair sums
`Sigma_psi0_psi6, ..., Sigma_psi5_psi11`, without making any physical admissibility claim.

## Inputs reused

1. `R18`
   - coefficient-class reduction of the canonical local diagonal residual `pair1` block entries.
2. `N467/P426`
   - six-class mode‑2 defect formulas (Re/Im coefficients on the six opposite‑pair sums).
3. `N473`
   - the theorem-level structural identity being audited.

## Probe method

Executed by:

- `fundamental_action_reconstruction/p433_current_strict_canonical_local_diagonal_mode2_defect_via_r18_pair1_entry_system_consistency_audit_probe.py`

For multiple random vectors $(\Sigma_0,\ldots,\Sigma_5)\in\mathbb{R}^6$:

1. compute `pair1` entries $(a_1,b_1,d_1)$ using the `R18` coefficient signatures,
2. compute `(Re(F2), Im(F2))` from the `N467/P426` six-class formulas,
3. compute `(Re(F2), Im(F2))` again as `(6*(a1-d1), 12*b1)`,
4. check equality within numerical tolerance.

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p433_current_strict_canonical_local_diagonal_mode2_defect_via_r18_pair1_entry_system_consistency_audit_probe.json`
- `fundamental_action_reconstruction/generated/p433_current_strict_canonical_local_diagonal_mode2_defect_via_r18_pair1_entry_system_consistency_audit_probe_summary.json`

## Result (no false pass)

On all audited random samples, the identity holds within tolerance (audit uses `tol=1e-10`, robust under the
~12-digit truncation of some `R18` transport coefficients):

```text
Re(F2) (six-class)  ==  Re(F2) (from R18 pair1 entries)
Im(F2) (six-class)  ==  Im(F2) (from R18 pair1 entries)
```

This probe audits only the *mathematical consistency* of the reductions.
It does not provide strict-derived coefficient values and does not decide `F2(d)` for the canonical FIN profile
(`T166` remains open).
