# P2592/S1542 APD Newton-Girard next even-moment sensitivity certificate

Status: `P2592_APD_NEWTON_GIRARD_NEXT_EVEN_MOMENT_SENSITIVITY_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Inherited product-parameter interval: `[300, 576]`.
- Internal eighth-shell formula: `74658.0 - 4*e4`.
- Central eighth-moment formula: `1303578.06546021 - 8*e4`.
- Lower shells constant but next shell varies: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

P2592 identifies the exact next missing coordinate after the P2591 Sturm interval.  Newton-Girard gives `p4 = 74658 - 4*e4` for the internal squared-offset eighth shell, so the same fixed second/fourth/sixth shell prefix leaves the next even shell linearly free.  Supplying that eighth shell would recover the product parameter, but the certificate does not derive that shell from strict nadsoliton dynamics.

## Recommended next honest step

Do not treat the eighth shell as a strict APD source merely because Newton-Girard recovers the free product parameter from it. P2592 proves the next even moment is exactly the missing selector coordinate for the P2591 interval; the honest next step is to derive that next shell/support law from strict nadsoliton dynamics, or else keep the product-parameter freedom explicit.

## Negative controls

No eighth-moment shell source, next-even-moment selector source, Newton-Girard selector source, product-parameter source, support-selection source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`68bb268ce311499f6d0be111c8951ffab041ec7b6ea281e11afb7a6a41a4c9fa`
