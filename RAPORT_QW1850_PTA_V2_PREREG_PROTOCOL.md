# RAPORT QW-1850: PTA V2 PREREG PROTOCOL

- Data UTC: 2026-03-02T23:06:54.272552+00:00
- Verdict: **PTA_V2_PREREG_FROZEN_EXTERNAL_CONFIRMATORY_PENDING**
- Protocol SHA256: `e5bcdc803f5587f790d9c1a70418463ed416760c4fcec72f6cd06b46a92b2f50`

## PTA V2 Thresholds (Pair-Level)
- mean_pair_mean_gain >= 0.04
- bootstrap_lower95_mean_pair_mean_gain >= 0.0
- prob_pair_mean_gain_positive >= 0.667
- one_sided_lower95_prob_pair_mean_gain_positive >= 0.6

## Reference (Design Dataset, Non-Confirmatory)
- mean_pair_mean_gain: 0.07115258005448141
- prob_pair_mean_gain_positive: 0.8202247191011236
- one_sided_lower95_prob_pair_mean_gain_positive: 0.7398456947151129

## Constraints
- External confirmatory dataset required.
- Reusing design dataset for final claim is forbidden.
- GW protocol inherited from QW-1839 (unchanged).

## Artifacts
- JSON: `report_qw1850_pta_v2_prereg_protocol.json`
