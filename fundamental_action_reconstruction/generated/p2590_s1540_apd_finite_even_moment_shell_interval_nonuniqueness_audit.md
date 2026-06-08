# P2590/S1540 APD finite even-moment shell interval nonuniqueness audit

Status: `P2590_APD_FINITE_EVEN_MOMENT_SHELL_INTERVAL_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Supports audited: `8`.
- Product-parameter grid: `[300.0, 340.0, 380.0, 420.0, 460.0, 500.0, 540.0, 576.0]`.
- Internal second/fourth/sixth shells: `30.0` / `354.0` / `4890.0`.
- Interval support nonunique: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

P2590 turns the previous finite examples into an interval-style certificate: the Vieta power-sum identity fixes the even-moment shell prefix while leaving the product parameter free.  Across the audited grid, full fixed-support moments still recover positive weights only after support is chosen, and selected off-node APD dynamics changes with the support.

## Recommended next honest step

Do not promote a finite even-moment shell prefix or its free product-parameter interval into an APD support source. P2590 keeps endpoints, mirror symmetry, cardinality, and fixed second/fourth/sixth shells across an audited product-parameter grid, while off-node APD dynamics remains support-dependent. The next honest step is a strict nadsoliton-derived support/density law or a genuine internal APD dynamics theorem.

## Negative controls

No finite even-moment shell source, product-parameter interval source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`83120eceafe6cea28ca7ac321b835e83297e14f9f4c3d451db9c2813093c7ae6`
