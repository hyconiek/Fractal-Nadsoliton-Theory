# P2652/S1602 beta=1 gauge unit-map validator and holdout covariance firewall

Status: `P2652_GAUGE_UNIT_MAP_VALIDATOR_READY_NO_REAL_PAYLOAD_NO_FALSE_PASS`

## Content-first anti-duplication audit

This packet greps gauge/unit-map, holdout covariance, fake unit-pass firewall, and source nonclosure content before adding the validator.

- `gauge_unit_map_content`: 117 hits
- `holdout_covariance_content`: 52 hits
- `fake_unit_pass_firewall_content`: 55 hits
- `source_nonclosure_content`: 14445 hits

## Gauge/unit-map validator schema

Required gauge contract keys: `['beta_after_gauge', 'beta_before_gauge', 'beta_gauge_status', 'distance_map_direction', 'eta', 'gauge_scale_a', 'raw_distances_are_not_interpreted_as_beta_one_distances', 'unit_map_source', 'unit_map_source_precommitted_before_holdout']`.
Required unit-aware measurement keys: `['far', 'far_in_beta_one_gauge', 'measured_log_tail_slope', 'measured_tail_ratio', 'near', 'near_in_beta_one_gauge', 'ratio_standard_error', 'raw_far', 'raw_near', 'slope_standard_error']`.
Forbidden unit-map sources: `['holdout_fit', 'per_pair_fit', 'post_unblind_refit', 'strict_kernel_target_fit', 'tail_ratio_fit']`.

## Firewall result

Covariant synthetic fixture admissible? `True`.
Raw-distance beta=1 substitution admissible? `False`.
Target-fit unit map admissible? `False`.
Firewall passes? `True`.

## Verdict

P2652 turns the P2651 beta=1 gauge contract into an executable covariance firewall for P2647/P2648.  A holdout payload may be evaluated in the beta=1 working gauge only if it carries a precommitted distance/unit map whose scaled distances match the locked P2646 pairs.  The covariant synthetic fixture passes that bookkeeping check, while raw-distance beta=1 substitution and target/holdout-fitted unit maps are rejected.  This is still validator readiness only: it supplies no real data and no typed nadsoliton unit-source theorem.

Decision: `GAUGE_UNIT_MAP_VALIDATOR_READY_BUT_NO_REAL_PAYLOAD_NO_UNIT_SOURCE_THEOREM_NO_CLOSURE`.
Real blind holdout payload loaded? `False`.
Unit-map source theorem exported now? `False`.
Beta source exported now? `False`.
Role-bearing L_total now? `False`.
ToE closure now? `False`.

## Professorial closure path

- Require every real P2647/P2648 payload to include beta_gauge_contract and per-measurement raw-to-beta-one distance coordinates.
- Reject any payload whose unit map is fitted from the holdout tail ratios, from the strict target, or per pair after unblinding.
- Run the statistical margin rule only after the gauge/unit covariance validator passes.
- Keep beta=1 labelled as gauge-fixed working normalization unless a typed nadsoliton metric/UV source theorem is proved independently.
- If real data pass P2647/P2648/P2652, rerun only the modified/compressed successor role row; do not revive unchanged legacy inverse-hierarchy or L_total closure.

## Next honest step

Prepare a real blinded payload with a precommitted raw-to-beta=1 unit map, run P2652 first, then P2647/P2648 exactly once; in parallel, attempt the independent typed nadsoliton metric/UV source theorem rather than fitting the unit map from the compression target.

## Negative exports

- `real_blind_holdout_payload_loaded`: `False`
- `empirical_confirmation_claimed`: `False`
- `unit_map_source_theorem_exported`: `False`
- `target_fit_unit_map_allowed`: `False`
- `raw_distance_beta_one_substitution_allowed`: `False`
- `beta_source_exported`: `False`
- `bridge_completion_exported`: `False`
- `role_transfer_revalidated`: `False`
- `q_w_2191_discharged`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `toe_closure_claimed`: `False`
