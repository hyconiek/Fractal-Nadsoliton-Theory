# P430 Current Strict `Aut(Z_12)`‑Invariant Diagonal Profile Mode‑2 Defect Orbit‑Reduction Audit Probe

Status: `P430_EXECUTED_CURRENT_STRICT_AUT_Z12_INVARIANT_DIAGONAL_PROFILE_MODE2_DEFECT_ORBIT_REDUCTION_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

`N471` proves a pure structural reduction:

```text
If a diagonal profile on n=12 is Aut(Z_12)-invariant (constant on the 6 quotient orbits from N455),
then F2(d) is real and equals an explicit orbit-linear combination.
```

`P430` provides a small numeric audit of this reduction on explicit Aut-invariant profiles (no false pass).

## Inputs reused

1. `N471` (orbit-reduction formula)
2. `N466` (pair1 diagonal O(2)-cut criterion via F2(d))

## Result

For all audited Aut-invariant example profiles:

1. computed `Im(F2(d))` is below tolerance,
2. computed `Re(F2(d))` matches the orbit-reduced formula from `N471` within tolerance,
3. the `pair1` anisotropy signature matches the `N466` expectation (cuts `O(2)` iff `|F2(d)|>tol`).

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p430_current_strict_aut_z12_invariant_diagonal_profile_mode2_defect_orbit_reduction_audit_probe.json`
- `fundamental_action_reconstruction/generated/p430_current_strict_aut_z12_invariant_diagonal_profile_mode2_defect_orbit_reduction_audit_probe_summary.json`

## Honest verdict (no false pass)

This probe audits only the **mathematical** orbit-reduction statement of `N471`.
It does **not** claim that any physically exported diagonal/local sector is Aut-invariant, and it does not discharge
`QW-2191`.

