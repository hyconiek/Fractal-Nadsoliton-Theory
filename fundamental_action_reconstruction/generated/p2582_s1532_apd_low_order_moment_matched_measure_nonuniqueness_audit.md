# P2582/S1532 APD low-order moment-matched measure nonuniqueness audit

Status: `P2582_APD_LOW_ORDER_MOMENT_MATCHED_MEASURE_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Boundary targets audited: `3`.
- Moment-matched positive measures audited: `3`.
- Moment orders matched: `[0, 1, 2]`.
- All targets low-order moment nonunique: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

Matching mass, first moment, and second raw moment does not determine the APD Gram measure.  P2582 constructs a positive moment-nullspace family with identical low-order moments; every audited measure preserves the APD nodes and endpoint slopes, but the selected off-node APD dynamics still changes.

## Recommended next honest step

Do not promote low-order APD measure moments into a strict source. P2582 constructs positive measures with identical mass, first moment, and second raw moment that still induce different Gram metrics and off-node APD dynamics. The next honest step is to derive the full APD measure density or all moment constraints from strict nadsoliton dynamics, rather than fixing only a finite moment prefix.

## Negative controls

No low-order APD moment-law source, moment-matched measure source, positive-measure source, L2-inner-product source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`1cea4d7e7cffb18569d51b524a7e3b30fb53fd867220faee7b3c2756a1f6b36e`
