# RAPORT QW-1926: BETA-CHANNEL POWER REQUIREMENTS MAP

- Data UTC: 2026-03-03T04:48:09.091130+00:00
- Verdict: **BETA_CHANNEL_POWER_REQUIREMENTS_MAP_READY**

## Baseline Inputs
- base_effect_min_from_qw1922: 0.9472
- base_sigma_conservative: 0.1000
- alpha_one_sided: 0.025

## Scenario Table (n_holdout)
- optimistic: eff=0.8051, sigma=0.0900, n80_eff=1, n90_eff=1, n95_eff=1, n80_act=400, n90_act=500, n95_act=500
- reference: eff=0.5683, sigma=0.1000, n80_eff=1, n90_eff=1, n95_eff=1, n80_act=400, n90_act=500, n95_act=500
- conservative: eff=0.4262, sigma=0.1250, n80_eff=1, n90_eff=1, n95_eff=2, n80_act=400, n90_act=600, n95_act=600
- stress: eff=0.2842, sigma=0.1600, n80_eff=3, n90_eff=4, n95_eff=5, n80_act=400, n90_act=800, n95_act=800

## Recommended Targets
- n_holdout_reference_power_0p90: 500
- n_holdout_conservative_power_0p90: 600
- n_total_pairs_reference: 1200
- n_total_pairs_conservative: 1600

## Artifacts
- JSON: `report_qw1926_beta_channel_power_requirements_map.json`
