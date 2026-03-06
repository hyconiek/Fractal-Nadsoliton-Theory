# RAPORT QW-2186: KTOTAL SPECTRAL STABILITY MARGIN GATE

- Date UTC: 2026-03-05T01:55:08.090990+00:00
- Verdict: **KTOTAL_SPECTRAL_STABILITY_MARGIN_GATE_PASS_STRICT_BRANCH_SCOPE**
- pass_count: `10/11`

## Core certificate
- `lambda_min(A)=0.331677209978` for `A=K_total + m0^2 I` with broken-floor `m0^2=1.013551972358`.
- Certified perturbation radius (Weyl): `||Delta||_2 < 0.331677209978` preserves PSD.

## Deterministic checks
- MC min lambda at safe radius: `0.139516159470`
- MC min lambda near boundary: `0.071742111429`
- witness above radius min lambda: `-0.016583860499`

## Scope
- Closed in branch scope for bounded symmetric operator-norm perturbations.
- Outside-scope classes remain explicit and unclaimed.

## Artifact
- JSON: `report_qw2186_ktotal_spectral_stability_margin_gate.json`
