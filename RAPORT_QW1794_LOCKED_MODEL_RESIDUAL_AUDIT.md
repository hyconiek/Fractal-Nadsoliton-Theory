# RAPORT QW-1794: LOCKED MODEL RESIDUAL AUDIT

- Data UTC: 2026-03-02T20:31:30.623354+00:00
- Cohort size: 91
- Best-fit params: A=-1.9767, q=1.4626, C=0.1646, sigma=0.9678
- rho(theta,resid)=0.1767, p_perm=0.0933
- rho(quality,resid)=-0.0136, p_perm=0.9021
- max |bin mean residual|=0.2801
- z diagnostics: mean=0.0006, std=0.9885, q95|z|=1.4942
- Verdict: **LOCKED_MODEL_RESIDUALS_STRUCTURED**

## Pass Flags
- theta_independence: False
- quality_independence: True
- angle_bin_flatness: False
- z_residual_diagnostics: True

## Artifacts
- JSON: `report_qw1794_locked_model_residual_audit.json`
