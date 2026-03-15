# P443 Current Strict QW‑2122 Scalar → Canonical Per‑Site Mapping Underdetermination Audit Probe

Status: `P443_EXECUTED_CURRENT_STRICT_QW2122_SCALAR_TO_PER_SITE_MAPPING_UNDERDETERMINATION_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`T168` remains open because the repo exports strict scalar vacuum closure values (`QW-2122/2123/2124`) but does **not**
export a strict-derived **lift** into the canonical per-site families:

```text
vpsi[0..11], g4[0..11], g6[0..11].
```

One recurring temptation is to treat the scalar outputs (e.g. `rho_star_sq`, `lambda_psi_strict`) as if they already
determined those per-site arrays “up to a trivial choice”.

`P443` audits (toy-level) that this is **not** the case:

```text
even if rho_star_sq is fixed to the strict QW-2122 value and g4/g6 are fixed uniformly by an explicit premise,
different per-site vacuum vectors with the same rho_star_sq can yield different diagonal sigma six-sums and different F2(d).
```

So scalar closure alone cannot discharge `T168` / `T167` / `T166` without exporting an additional strict mapping/selector
ingredient (or a theorem forcing symmetry ⇒ `F2(d)=0`).

No strict-derived value claim is made.

## Probe method

Executed by:

- `fundamental_action_reconstruction/p443_current_strict_qw2122_scalar_to_per_site_mapping_underdetermination_audit_probe.py`

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p443_current_strict_qw2122_scalar_to_per_site_mapping_underdetermination_audit_probe.json`
- `fundamental_action_reconstruction/generated/p443_current_strict_qw2122_scalar_to_per_site_mapping_underdetermination_audit_probe_summary.json`

## Result (no false pass)

The probe constructs **two** toy instantiations that share:

1. the strict scalar norm `rho_star_sq` from `QW-2122` (broken branch), and
2. the same uniform local self-coupling arrays (explicit premise),

but differ only in the per-site vacuum vector `vpsi` while keeping `Σ_i vpsi_i^2 = rho_star_sq`.

It then computes:

- the Yukawa-free diagonal residual profile and six opposite-pair sums using the `N477` `K_total`-numeric formula, and
- the induced `F2(d)` via the `N467` six-sum reduction.

One case yields `F2(d)=0` (uniform `vpsi`), and the other yields `F2(d)≠0` (a two-site perturbation),
demonstrating that the scalar outputs do not canonically fix the per-site arrays.

## Relation to theorems/targets

- strict scalar vacuum closure (inputs): `QW-2122/2123/2124`
- kernel channel + floor (inputs): `R14/R15`
- Yukawa-free `K_total`-numeric diagonal reduction (used): `N477`
- opposite-pair six-sum reduction (used): `N467`
- missing strict provider target: `T168`
- diagonal decision target: `T166`

