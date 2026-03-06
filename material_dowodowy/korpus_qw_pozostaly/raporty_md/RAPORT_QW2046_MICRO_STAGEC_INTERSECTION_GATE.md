# RAPORT QW-2046: MICRO STAGE-C INTERSECTION GATE

- Data UTC: 2026-03-03T21:31:10.017655+00:00
- Verdict: **MICRO_STAGEC_INTERSECTION_GATE_PARTIAL**
- Readiness: **MICRO_TO_STAGEC_BRIDGE_PARTIAL**
- pass_count: 5/7

## Micro Support (from QW-2045)
- beta CI95: [0.056034, 37.175785]
- eta CI95: [1.450000, 3.000000]

## Stage-C Intersection
- pass_examples_count: 10
- intersection_count: 10
- selected_from_intersection: True
- selected kernel omega/phi/beta/eta: 0.185750 / 0.162500 / 1.000000 / 1.800000

## External Primary
- corr: 0.054705 (q95=0.035421, p=0.007399)
- gain: 0.000571 (q95=-0.001273, p=0.007399)

## External Stress
- corr: 0.330346 (q95=0.051744, p=0.000200)
- gain: 0.055985 (q95=-0.031792, p=0.000200)

## Flags
- stagec_candidate_all_pass: True
- micro_beta_overlap: True
- micro_eta_overlap: True
- external_primary_all_pass: True
- external_stress_soft_pass: True
- micro_pointwise_bins_ge_6: False
- micro_phase_condition_strength_ge_0p75: False

## Required Next Step
- ADD_SIGNED_PHASE_TORSION_OBSERVABLE_AND_RE-RUN_POINTWISE_DERIVATION

## Artifacts
- JSON: `report_qw2046_micro_stagec_intersection_gate.json`
