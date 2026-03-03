# RAPORT QW-1968: REFINED KERNEL ROBUSTNESS BOOTSTRAP GATE

- Data UTC: 2026-03-03T08:28:40.038681+00:00
- Verdict: **FRAGILE_PASS_NOT_YET_LOCKABLE**

## Baseline (from QW-1967 best point)
- flavor CKM/PMNS mean rel%: 12.852/11.734
- GW auc/adv/sep/gap: 0.8233/0.3221/0.002699/0.002282
- all_pass: True

## Local Robustness (fixed q_nu, no sector retune)
- radius=0.0025: pass_rate=100.000% (12000/12000)
- radius=0.0050: pass_rate=100.000% (12000/12000)
- radius=0.0100: pass_rate=100.000% (12000/12000)
- radius=0.0200: pass_rate=95.433% (11452/12000)
- radius=0.0500: pass_rate=70.992% (8519/12000)

## GW Bootstrap (n=5000)
- GW pass rate: 71.92%
- Triad pass rate (mass+flavor fixed + GW bootstrap): 71.92%

## Sensitivity Top 5 (one-at-a-time)
- chi_im: sensitivity_abs=0.000, step=0.002000
- diag_iso: sensitivity_abs=0.000, step=0.000300
- diag_q_coeff: sensitivity_abs=0.000, step=0.000400
- diag_sector: sensitivity_abs=0.000, step=0.000400
- lambda_mix: sensitivity_abs=0.000, step=0.004000

## Required Next Step
- INCREASE_PASS_VOLUME_WITHOUT_KERNEL_RETUNE_THEN_REPEAT_QW1968

## Artifacts
- JSON: `report_qw1968_refined_kernel_robustness_bootstrap_gate.json`
