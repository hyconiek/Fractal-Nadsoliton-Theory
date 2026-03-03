# RAPORT QW-1935: HEAD-TO-HEAD NADSOLITON VS HD PROXY

- Data UTC: 2026-03-03T05:42:35.710199+00:00
- Reparam triad: omega=0.373414, phi=-1.310234, beta=0.615938, eta=2.80
- Verdict: **HEAD2HEAD_MIXED_PRIMARY_ONLY_REPARAM_BETTER**

## Primary Summary
- reparam corr/gain/rmse medians: 0.0706/0.0024/0.1836
- hd corr/gain/rmse medians: 0.0569/0.0016/0.1839
- reparam_vs_hd delta_rmse median: 0.000153
- reparam_vs_hd win_rate/sign_p: 0.9167/1.794e-05

## Stress Summary
- reparam corr/gain/rmse medians: 0.3177/0.0508/0.2695
- hd corr/gain/rmse medians: 0.4112/0.0879/0.2595
- reparam_vs_hd delta_rmse median: -0.009758
- reparam_vs_hd win_rate/sign_p: 0.0000/1

## Flags
- primary_reparam_beats_hd_win_rate_ge_0p80: True
- primary_reparam_beats_hd_sign_p_le_0p05: True
- stress_reparam_beats_hd_win_rate_ge_0p80: False
- stress_reparam_beats_hd_sign_p_le_0p05: False

## Caveat
- HD here is a pragmatic SM+GR PTA proxy with affine nuisance, not a full detector-noise Bayesian pipeline.

## Required Next Step
- CHARACTERIZE_REGIME_SPLIT_AND_DEFINE_DOMAIN_OF_VALIDITY_BEFORE_STRONG_TOE_CLAIM

## Artifacts
- JSON: `report_qw1935_head2head_nadsoliton_vs_hd_proxy.json`
