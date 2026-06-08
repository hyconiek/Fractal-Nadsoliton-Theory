# P2586/S1536 APD mirror-symmetric support selector nonuniqueness audit

Status: `P2586_APD_MIRROR_SYMMETRIC_SUPPORT_SELECTOR_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Supports audited: `3`.
- Mirror center: `5.5`.
- All supports pair-symmetric: `True`.
- Mirror support nonunique: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

Mirror/pair symmetry is stronger than the centroid guard in P2585, but it still does not determine the APD support.  P2586 keeps reflection symmetry around the same center, recovers weights conditionally on each fixed support, and still obtains support-dependent off-node APD dynamics.

## Recommended next honest step

Do not promote mirror/reflection symmetry of finite support into an APD support source. P2586 keeps pair symmetry about the same center and still obtains distinct support choices with conditionally unique weights and support-dependent APD dynamics. The next honest step is to derive the actual support/density law from strict nadsoliton dynamics, not another geometric support preference.

## Negative controls

No APD mirror-support source, reflection-symmetry support source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`19f624b3a24aa76a5df24cdd57c738bd02114942d796e6d06e2954cdc0a34675`
