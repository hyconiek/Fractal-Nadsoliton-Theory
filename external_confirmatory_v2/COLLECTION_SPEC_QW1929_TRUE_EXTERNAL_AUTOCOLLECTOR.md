# COLLECTION SPEC QW-1929: TRUE EXTERNAL AUTOCOLLECTOR

## Scope
Build `external_confirmatory_v2/beta_channel_true_external` from:
1. `NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz` raw timing archive,
2. GWOSC GWTC event catalog (`https://www.gw-openscience.org/eventapi/json/GWTC/`).

## Non-Negotiable Rule
Do NOT read any internal proxy/rebuild pair tables.

## Required Outputs
- `manifest_beta_channel.json`
- `beta_channel_pairs.csv`
- `intervention_events.csv`
- `protocol_freeze.json`

## Pair Schema
- pair_id, theta_deg, hxy, f_std, f_autoc1, f_switch, f_slope, intervention_id, regime

## Externality Guardrail
- provider must not contain INTERNAL / INTERNAL_PROXY,
- statement must include `independent` and `external`,
- statement must not include `not independent` or `internal proxy`.

## Gate
- run `python3 QW_1927_TRUE_EXTERNAL_BETA_CHANNEL_READINESS_GATE.py`
- expected verdict: `TRUE_EXTERNAL_BETA_CHANNEL_READY`
