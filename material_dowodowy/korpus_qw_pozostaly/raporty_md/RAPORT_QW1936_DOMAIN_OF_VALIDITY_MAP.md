# RAPORT QW-1936: DOMAIN OF VALIDITY MAP

- Data UTC: 2026-03-03T05:44:12.030834+00:00
- Reparam triad: omega=0.373414, phi=-1.310234, beta=0.615938, eta=2.80
- Verdict: **DOMAIN_OF_VALIDITY_MAP_COMPLETE**

## Global Delta MSE (HD - Reparam)
- Primary median [q10, q90]: 0.000056 [0.000015, 0.000111]
- Primary global win_rate(reparam): 0.9167
- Stress median [q10, q90]: -0.005178 [-0.006425, -0.004421]
- Stress global win_rate(reparam): 0.0000

## Domain Statement: Primary
- reparam_dominant_bins: ['30-60']
- hd_dominant_bins: []
- uncertain_bins: ['0-30', '60-90', '90-120', '120-150', '150-180']

## Domain Statement: Stress
- reparam_dominant_bins: []
- hd_dominant_bins: ['0-30', '30-60', '60-90', '90-120', '120-150', '150-180']
- uncertain_bins: []

## Required Next Step
- FREEZE_DOMAIN_CONDITIONAL_PREDICTIONS_AND_PREPARE_EXTERNAL_CONFIRMATORY_PROTOCOL

## Artifacts
- JSON: `report_qw1936_domain_of_validity_map.json`
