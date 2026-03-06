# RAPORT QW-1793: MODEL LOCK-IN GATE

- Data UTC: 2026-03-02T20:29:26.824238+00:00
- Global score: 0.991
- Hard gate: **PASS**
- Readiness: **MODEL_LOCKIN_CONFIRMED**

## Checks
- Campaign robustness (1786): PASS | score=0.966 | note=EMPIRICAL_CAMPAIGN_ROBUSTNESS_CONFIRMED
- Hyperparameter lock (1788+1789): PASS | score=1.000 | note=q_width=0.2, criteria=K0_base
- Negative-control extension rejected (1790): PASS | score=1.000 | note=delta_M3-M2=-0.8776
- Likelihood alternatives vs baseline (1791+1792): PASS | score=1.000 | note=delta_t-gauss=-0.0948, delta_hetero-homo=-0.1012

## Locked Operational Protocol
- signal_model: M2_reparam_HDq_plus_constant
- noise_model: homoscedastic_gaussian
- fraction: 0.95
- q_width: 0.2
- cohort: K0_base

## Artifacts
- JSON: `report_qw1793_model_lockin_gate.json`
