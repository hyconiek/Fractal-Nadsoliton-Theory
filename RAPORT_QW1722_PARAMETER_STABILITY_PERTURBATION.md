# RAPORT QW-1722: PARAMETER STABILITY PERTURBATION

- Data UTC: 2026-03-02T16:01:03.607994+00:00
- Seed: 1722
- Werdykt: **PARAMETERS_UNSTABLE_UNDER_PERTURBATION**

## 1) Mass block
- pass: False
- test_mean median/p90: 8.99% / 9.30%
- test_max p90: 13.55%
- param stats (normalized IQR, sign flip):
  - l1: norm_iqr=0.060, sign_flip=0.271
  - l2: norm_iqr=0.001, sign_flip=0.000
  - l3: norm_iqr=0.001, sign_flip=0.000

## 2) Flavor block
- pass: False
- canonical median errors CKM/PMNS: 23.08% / 13.08%
- perturbed fit score median: 45.73
- param stats (normalized IQR, sign flip):
  - a1: norm_iqr=0.041, sign_flip=0.000
  - a2: norm_iqr=0.331, sign_flip=0.083
  - a3: norm_iqr=0.000, sign_flip=0.000
  - a4: norm_iqr=0.436, sign_flip=0.242
  - delta: norm_iqr=0.031, sign_flip=0.000

## Artefakty
- JSON szczegolowy: `report_qw1722_parameter_stability_perturbation.json`
