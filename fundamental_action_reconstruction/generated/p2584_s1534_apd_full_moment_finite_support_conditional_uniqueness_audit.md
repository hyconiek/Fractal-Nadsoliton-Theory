# P2584/S1534 APD full-moment finite-support conditional uniqueness audit

Status: `P2584_APD_FULL_MOMENT_FINITE_SUPPORT_CONDITIONAL_UNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Supports audited: `3`.
- Full moment orders: `[0, 1, 2, 3]`.
- Fixed-support weights conditionally unique: `True`.
- Targets support-dependent: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

Full moments on a fixed finite support recover the weights by a nonsingular Vandermonde system.  This is a useful conditional lemma, but not a strict source: changing the still-unsourced support changes the induced Gram metric and selected off-node APD dynamics while preserving the same APD node and boundary constraints.

## Recommended next honest step

Use full finite-support moments only as a conditional lemma: once the support is strict-sourced, Vandermonde moments recover the weights. P2584 shows that different strict-unsourced supports each have internally unique weights but still select different APD dynamics. The next honest step is to derive the APD measure support/density from strict nadsoliton dynamics, not merely solve weights on a chosen support.

## Negative controls

No finite-support source, full-moment-law source, support-selection source, positive-measure source, L2-inner-product source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`569234a1aa9a377650bcadf75dfad88847e92cd11b132044fd5153c39273abc1`
