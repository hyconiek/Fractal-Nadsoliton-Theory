# COLLECTION SPEC QW-1925: TRUE EXTERNAL BETA CHANNEL

This is a strict operational checklist for collecting true external intervention data
to resolve Stage-B beta identifiability in FIN-ToE.

## Hard Requirements
- External independent provider only (no INTERNAL_PROXY / INTERNAL tokens).
- Manifest must include roles: beta_pairs, intervention_events, protocol_freeze.
- At least two exogenous intervention cohorts.
- Freeze protocol before metric execution.

## File Package
- manifest_beta_channel.json
- beta_channel_pairs.csv
- intervention_events.csv
- protocol_freeze.json

## Minimal Signal Targets
- holdout_effect_beta >= 0.6
- holdout_effect_omega <= 0.25
- holdout_contrast >= 0.35
- holdout_contrast_boot_q05 >= 0.2

## Power Targets
- n_holdout min (80%): 400
- n_holdout min (90%): 500
- n_total_pairs recommended: 1200

## Execution After Collection
- Run QW-1927 readiness gate (schema/externality/protocol lock).
- Run blind intervention evaluation with frozen thresholds.
- Re-run Stage-B gate with updated evidence.
