# P2572/S1522 APD boundary-penalty selector continuum audit

Status: `P2572_APD_BOUNDARY_PENALTY_SELECTOR_CONTINUUM_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Candidate principle: `fix k=2 and uniform bulk measure, then vary one endpoint slope penalty t >= 0`.
- Stationarity formula: `lambda(t)=-(B_bulk+t*B_endpoint)/(A_bulk+t*A_endpoint)`.
- Left penalty distinct minimizers: `9`.
- Right penalty distinct minimizers: `9`.
- All penalty rows preserve APD nodes: `True`.
- Finite APD values select boundary penalty: `False`.

## Interpretation

With `k=2` and uniform bulk measure fixed, adding a nonnegative endpoint slope penalty gives the explicit one-parameter stationarity law `lambda(t)=-(B_bulk+t*B_endpoint)/(A_bulk+t*A_endpoint)`.  The audited grid already yields multiple left- and right-endpoint minimizers, all preserving the same finite APD nodes.  Boundary penalty strength is therefore a continuous unsourced selector, not strict A/P/D dynamics.

## Recommended next honest step

Do not choose an endpoint penalty by convenience or by post-hoc fit. The next honest step is to derive boundary conditions or boundary penalties from the strict nadsoliton APD action; otherwise even fixed k=2 and uniform bulk measure leave a continuous family of conditional APD selectors.

## Negative controls

No APD boundary penalty source, boundary selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`85a4157be2e8a213415ea61f2c4473f45006bb2ee609a28f0ddb1e0f01c5e4d7`
