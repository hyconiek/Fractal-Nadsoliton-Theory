# RAPORT QW-1926: BETA-CHANNEL POWER REQUIREMENTS MAP

- Data UTC: 2026-03-03T19:06:54.097440+00:00
- Verdict: **BETA_CHANNEL_POWER_REQUIREMENTS_MAP_READY**

## Baseline Inputs
- base_effect_min_from_qw1922: 1.3174
- base_sigma_conservative: 0.1000
- alpha_one_sided: 0.025

## Scenario Table (n_holdout)
- optimistic: eff=1.1198, sigma=0.0900, n80_eff=1, n90_eff=1, n95_eff=1, n80_act=400, n90_act=500, n95_act=500
- reference: eff=0.7904, sigma=0.1000, n80_eff=1, n90_eff=1, n95_eff=1, n80_act=400, n90_act=500, n95_act=500
- conservative: eff=0.5928, sigma=0.1250, n80_eff=1, n90_eff=1, n95_eff=1, n80_act=400, n90_act=600, n95_act=600
- stress: eff=0.3952, sigma=0.1600, n80_eff=2, n90_eff=2, n95_eff=3, n80_act=400, n90_act=800, n95_act=800

## Recommended Targets
- n_holdout_reference_power_0p90: 500
- n_holdout_conservative_power_0p90: 600
- n_total_pairs_reference: 1200
- n_total_pairs_conservative: 1600

## Artifacts
- JSON: `report_qw1926_beta_channel_power_requirements_map.json`
