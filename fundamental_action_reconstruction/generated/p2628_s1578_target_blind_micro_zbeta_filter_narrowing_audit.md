# P2628/S1578 target-blind micro Z_beta filter narrowing audit

Status: `P2628_TARGET_BLIND_MICRO_ZBETA_FILTER_NARROWING_NO_STRICT_SOURCE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE`

## Content-first anti-duplication audit

This packet greps research content and quality-filter language, not only ticket numbers or names.

- `micro_pointwise_quality_content`: 4747 hits
- `target_blind_filter_content`: 56 hits
- `zbeta_distribution_content`: 6777 hits
- `tolerance_bridge_content`: 27 hits
- `closure_guard_content`: 13947 hits

## Finite target-blind filter theorem

Within the fixed target-blind rectangular quality-filter class over QW-2048 bins, no filter with at least 3 support bins has both median Z_beta within 1% of 100 and q25-q75 inside the same 1% envelope.

Evaluated filters: `1086`; strict accepting filters: `0`; practical interval accepting filters: `0`.

Best median-only row:

- thresholds: `{'n_min': 8, 'phase_min_median_min': 0.8375652609364894, 'rmse_median_max': 0.11595265662226857}`
- support bins: `13` with `d_bins=[2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]`
- median Z_beta: `109.761471`
- q25/q75: `[28.508074, 364.031509]`
- relative median error to 100: `0.097615`

## Verdict

The target-blind quality-filter class does not export `positive_beta_renormalization_source`.  It also does not repair P2625/P2620, role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure.

## Recommended next honest step

Stop trying to squeeze exact Z_beta=100 from post-hoc quality filters.  The next proof-quality move is to derive a target-independent micro operator or conservation/normalization identity whose value is Z_beta before comparison with K_strict_gate; if that cannot be done, declare only an approximate finite-domain effective bridge with explicit epsilon/domain and keep role-transfer closed.

Fingerprint: `a1fdecacfc4ba68288201ffc0daac0940232a1cafcf923e6da94221e632b2633`
