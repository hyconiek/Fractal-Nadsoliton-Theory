# RAPORT QW-2049: SPECTRAL MICRO STAGE-C INTERSECTION GATE

- Data UTC: 2026-03-04T07:19:37.630814+00:00
- Verdict: **SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS**
- Readiness: **TOE_INTERNAL_BRIDGE_STRICTLY_CLOSED_PENDING_EXTERNAL_MULTITEAM_AUDIT**
- pass_count: 7/7

## Micro Support (from QW-2048)
- beta CI95: [0.000054, 106.682605]
- eta CI95: [1.000000, 3.000000]

## Stage-C Intersection
- pass_examples_count: 10
- intersection_count: 10
- selected kernel omega/phi/beta/eta: 0.185750 / 0.162500 / 1.000000 / 1.800000

## External Primary
- corr: 0.054705 (q95=0.036683, p=0.006999)
- gain: 0.000571 (q95=-0.001152, p=0.006999)

## External Stress
- corr: 0.330346 (q95=0.050699, p=0.000200)
- gain: 0.055985 (q95=-0.032107, p=0.000200)

## Flags
- stagec_candidate_all_pass: True
- micro_beta_overlap: True
- micro_eta_overlap: True
- external_primary_all_pass: True
- external_stress_soft_pass: True
- micro_pointwise_bins_ge_6: True
- micro_phase_condition_strength_ge_0p75: True

## Required Next Step
- FREEZE_SPECTRAL_MICRO_BRIDGE_AND_PREPARE_EXTERNAL_MULTITEAM_PACKAGE

## Artifacts
- JSON: `report_qw2049_spectral_micro_stagec_intersection_gate.json`
