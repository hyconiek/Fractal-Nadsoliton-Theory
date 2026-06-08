# P2581/S1531 APD Gram-measure moment dependence audit

Status: `P2581_APD_GRAM_MEASURE_MOMENT_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Boundary targets audited: `3`.
- Positive measures audited: `5`.
- All targets measure dependent: `True`.
- All Gram metrics positive definite: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

Positive L2/Gram metrics are covariance-compatible candidates, but the positive measure is extra data.  P2581 audits five positive measures and finds that they preserve the same finite APD nodes and boundary slopes while selecting different off-node APD dynamics.  Thus the missing source is pushed from coordinates to the measure/moment law.

## Recommended next honest step

Do not treat an L2/Gram APD inner product as strict unless the positive measure and its moments are derived from strict nadsoliton dynamics. P2581 shows that multiple positive measures give positive-definite, covariance-compatible Gram metrics while selecting different APD off-node dynamics under the same finite values and boundary targets. The next honest step is to derive the APD measure density or kinetic weight, not merely choose a quadrature rule.

## Negative controls

No APD positive-measure source, Gram-moment source, L2-inner-product source, measure-selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`2595b860fc2bcfd118f9ab62b3213a1598255d384a1ae963c0d70e34c98f36f9`
