# P2585/S1535 APD support-geometry selector nonuniqueness audit

Status: `P2585_APD_SUPPORT_GEOMETRY_SELECTOR_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Supports audited: `3`.
- Shared support constraints: cardinality `4`, endpoints `(0.25, 10.75)`, centroid `5.375`.
- Fixed-support weights conditionally unique: `True`.
- Support geometry nonunique: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

Cardinality, endpoints, centroid, and fixed-support full moments still do not determine the APD support.  P2585 keeps those support-geometry constraints fixed, recovers weights conditionally on each support, and still obtains support-dependent off-node APD dynamics.

## Recommended next honest step

Do not promote finite-support geometry constraints into an APD support source. P2585 fixes cardinality, endpoints, and centroid, and still obtains distinct support choices with conditionally unique Vandermonde weights that select different APD off-node dynamics. The next honest step is to derive the actual APD support/density law from strict nadsoliton dynamics.

## Negative controls

No APD support-geometry source, endpoint/centroid source, support-cardinality source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`6775393fd64dfc7ec696edd75d0883122b2ea79010f2f983e19b121017e8dbfb`
