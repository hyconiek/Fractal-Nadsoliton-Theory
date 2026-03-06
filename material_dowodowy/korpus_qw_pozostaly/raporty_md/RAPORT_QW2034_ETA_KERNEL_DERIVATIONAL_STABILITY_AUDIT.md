# RAPORT QW-2034: ETA KERNEL DERIVATIONAL STABILITY AUDIT

- Data UTC: 2026-03-03T20:22:10.589690+00:00
- Readiness: **DERIVATIONAL_STABILITY_ACCEPTABLE**
- Verdict: **ETA_KERNEL_DERIVATIONAL_STABILITY_PASS**
- pass_count: 5/6

## Target Kernel (QW-2030)
- omega=0.233750, phi=-0.137500, beta=1.170000, eta=1.800000

## Bootstrap CI95
- omega q02.5/q50/q97.5: 0.020000 / 0.224116 / 0.247708
- phi q02.5/q50/q97.5: -0.257614 / 1.350505 / 2.929739
- beta q02.5/q50/q97.5 (std): 0.567478 / 0.761376 / 0.925256 (0.109906)
- eta q02.5/q50/q97.5 (std): 1.800000 / 1.800000 / 1.800000 (0.000000)

## Flags
- target_omega_in_bootstrap_ci95: True
- target_phi_in_bootstrap_ci95: True
- target_beta_in_bootstrap_ci95: False
- target_eta_in_bootstrap_ci95: True
- beta_std_le_0p60: True
- eta_std_le_0p50: True

## Required Next Step
- LOCK_DERIVATIONAL_APPENDIX_FOR_FROZEN_KERNEL

## Artifacts
- JSON: `report_qw2034_eta_kernel_derivational_stability_audit.json`
