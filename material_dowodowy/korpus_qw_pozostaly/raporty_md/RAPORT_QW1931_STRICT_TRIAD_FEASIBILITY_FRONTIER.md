# RAPORT QW-1931: STRICT TRIAD FEASIBILITY FRONTIER

- Data UTC: 2026-03-03T05:32:34.361645+00:00
- Verdict: **STRICT_TRIAD_FEASIBILITY_FAIL**
- strict_pass_count: 0

## Selected Candidate
- lambda_beta: 2.0
- omega: 0.059969
- phi: 0.890000
- beta: 1.062500
- rel_loss_vs_unconstrained: 0.7570
- primary corr/gain ratios: 0.9839 / 0.6714
- stress corr/gain ratios: 0.9650 / 0.9293

## Strict Flags (selected)
- beta_le_1p20: True
- rel_loss_le_0p35: False
- primary_corr_ratio_ge_0p90: True
- primary_gain_ratio_ge_0p90: False
- stress_corr_ratio_ge_0p80: True
- stress_gain_ratio_ge_0p80: True
- omega_interior: True

## Required Next Step
- REFORMULATE_MICROMODEL_OR_PARAMETERIZATION_TO_REMOVE_STRUCTURAL_TRADEOFF

## Artifacts
- JSON: `report_qw1931_strict_triad_feasibility_frontier.json`
