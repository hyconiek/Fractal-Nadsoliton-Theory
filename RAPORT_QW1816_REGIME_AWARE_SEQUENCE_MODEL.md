# RAPORT QW-1816: REGIME-AWARE SEQUENCE MODEL

- Data UTC: 2026-03-02T21:42:24.807928+00:00
- Regime-aware cohort size: 90
- Mean windows per pair: 15.59
- Features: f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch
- Regime counts: [30, 30, 30]
- Full logB M2/M2E/M2ER: 0.1061 / 15.8678 / 14.6801
- Full delta M2ER-M2E: -1.1877
- P(M2ER>M2E): 0.500
- P(M2ER>flat): 1.000
- Std delta (M2E-M2): 1.560
- Std delta (M2ER-M2): 1.666
- Std reduction vs M2E: -0.106
- Verdict: **REGIME_AWARE_SEQUENCE_MODEL_WEAK**

## Pass Flags
- full_gain_vs_m2e: False
- replication_gain_vs_m2e: False
- dispersion_control: False
- stability_improvement_vs_m2e: False

## Artifacts
- JSON: `report_qw1816_regime_aware_sequence_model.json`
