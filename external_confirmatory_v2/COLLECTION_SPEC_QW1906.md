# COLLECTION SPEC QW-1906

This is an operational checklist for acquiring a truly external frozen PTA/GW package.

## Hard Requirements
- PTA pairs >= 1200 (preferred 2000)
- Keep thresholds fixed from QW-1850 (no retuning)
- Independent provider and explicit externality statement

## PTA Signal Targets
- mean_pair_mean_gain >= 0.04
- prob_pair_mean_gain_positive >= 0.667
- one_sided_lower95_prob_positive >= 0.6

## Proxy Translation (planning)
- alpha_required_for_80pct_power ~= 6.0
- target hxy-feature slope ~= 0.2744
- target hxy-feature corr ~= 0.9724

## Next Execution
- Freeze manifest + file hashes
- Run QW-1852 -> QW-1853 -> QW-1902
