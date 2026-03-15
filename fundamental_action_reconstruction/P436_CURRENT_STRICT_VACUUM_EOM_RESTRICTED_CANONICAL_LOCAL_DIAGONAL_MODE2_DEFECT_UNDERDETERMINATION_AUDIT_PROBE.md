# P436 Current Strict Vacuum‑EoM‑Restricted Canonical Local‑Diagonal Mode‑2 Defect Underdetermination Audit Probe

Status: `P436_EXECUTED_CURRENT_STRICT_VACUUM_EOM_RESTRICTED_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_UNDERDETERMINATION_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

`N472/P431` prove that the canonical diagonal/local mode‑2 defect target (`T166`) is underdetermined at the exported
coefficient‑class level (`R15`).

After `N474/N475`, one might still be tempted to think:

```text
maybe imposing the exported constant-vacuum EoM stationarity already decides F2(d) for the canonical diagonal profile
```

`P436` audits (toy‑level) that this is **not** the case:

```text
even after restricting to constant‑vacuum stationarity (and vpsi_i≠0 so N474 applies),
F2(d) can still be 0 or ≠0 depending only on remaining free local coefficient families.
```

No strict-derived value claim is made.

## Probe method

Executed by:

- `fundamental_action_reconstruction/p436_current_strict_vacuum_eom_restricted_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe.py`

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p436_current_strict_vacuum_eom_restricted_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe.json`
- `fundamental_action_reconstruction/generated/p436_current_strict_vacuum_eom_restricted_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe_summary.json`

## Result (no false pass)

The probe constructs **two** toy instantiations which both satisfy constant‑vacuum stationarity by per‑site solving of
`m2_psi{i}`:

1. uniform local couplings ⇒ constant diagonal profile ⇒ `F2(d)=0`,
2. a single-site `g4_psi0` defect (all else equal) ⇒ non‑translation‑invariant diagonal profile ⇒ `F2(d)≠0`.

Therefore, constant‑vacuum stationarity alone does **not** decide `F2(d)` for the canonical diagonal sector, and cannot
discharge `T166` without additional strict‑derived structure.

## Relation to theorems/targets

- Yukawa elimination premise: `N474`
- opposite‑pair sum rewrite: `N475`
- vacuum‑restricted underdetermination theorem: `N476`
- decision target: `T166`

