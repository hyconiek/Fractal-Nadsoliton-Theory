# RAPORT QW-1792: HETEROSCEDASTIC REPARAM

- Data UTC: 2026-03-02T20:28:20.994782+00:00
- Operational pairs: 91
- Full homo logB(reparam vs flat): 0.2300
- Full hetero logB(reparam vs flat): 0.1288
- Full gain hetero-homo: -0.1012
- P(hetero>0): 1.000
- P(hetero>homo): 0.214
- P(hetero delta vs legacy > 0): 1.000
- Std gain hetero-homo: 0.051
- Verdict: **HETEROSCEDASTIC_REPARAM_PARTIAL**

## Pass Flags
- full_positive: True
- full_gain_over_homo: False
- rep_positive: True
- rep_gain_over_homo: False
- rep_delta_vs_legacy_positive: True
- dispersion_control: True

## Artifacts
- JSON: `report_qw1792_heteroscedastic_reparam.json`
