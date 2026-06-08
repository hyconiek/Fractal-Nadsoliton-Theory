# P2573/S1523 APD boundary-penalty inverse-target tunability audit

Status: `P2573_APD_BOUNDARY_PENALTY_INVERSE_TARGET_TUNABILITY_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Candidate principle: `given a desired APD lambda in the endpoint-penalty interval, solve the boundary penalty t that realizes it`.
- Inverse formula: `t(lambda_target)=-(lambda_target*A_bulk+B_bulk)/(lambda_target*A_endpoint+B_endpoint)`.
- Inverse rows: `10`.
- All targets recovered: `True`.
- All solved penalties nonnegative: `True`.
- All inverse rows preserve APD nodes: `True`.
- Finite APD values select target/penalty law: `False`.

## Interpretation

The P2572 stationarity law can be inverted: for any audited target `lambda` between the uniform-bulk minimizer and the endpoint-penalty limit, the boundary penalty `t` is solved explicitly and recovers the target while preserving all finite APD nodes.  This makes endpoint penalties post-hoc tunable unless their law is derived before selecting the APD dynamics.

## Recommended next honest step

Do not fit boundary penalties inversely to obtain a desired APD interpolation. The next honest step is to derive admissible boundary conditions and penalty strengths before choosing the APD member; otherwise the boundary term is a post-hoc selector and not strict A/P/D dynamics.

## Negative controls

No APD boundary penalty law source, inverse boundary selector source, target-lambda source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`03cce5b09d814f43ca4cd165a441f0aa1caa75031b0132935120a7074ddbf851`
