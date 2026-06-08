# P2570/S1520 APD Sobolev roughness selector order-dependence audit

Status: `P2570_APD_SOBOLEV_ROUGHNESS_SELECTOR_ORDER_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Candidate principle: `choose q_lambda minimizing a Sobolev roughness J_k(lambda)=int_0^11 |d^k q_lambda/dx^k|^2 dx on the P2569 interpolation family`.
- P2569 interpolation family inherited: `True`.
- Derivative orders audited: `[0, 1, 2, 3, 4, 12]`.
- Orders selecting nonzero lambda: `[0, 1, 2, 3, 4]`.
- Orders selecting base lambda zero: `[12]`.
- Roughness order changes selector: `True`.
- Finite APD values select Sobolev order: `False`.

## Interpretation

On the P2569 family `q_lambda=q_interp+lambda*prod_{d=0}^{11}(x-d)`, the audited Sobolev roughness objective has a closed-form quadratic in `lambda`.  Different derivative orders choose different minimizers: low-order roughness selects nonzero lambda values while the twelfth-derivative roughness selects the base interpolant.  Because every selected member still preserves the finite APD nodes, the derivative order/measure is an extra source obligation.

## Recommended next honest step

Do not promote a minimum-roughness APD interpolation to strict dynamics unless the derivative order, measure, and boundary class are derived from nadsoliton dynamics. The next honest step is to derive a sourced APD variational functional; if no such functional is available, the A/P/D residual atom remains open even when finite APD values are exact.

## Negative controls

No APD roughness selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`8544032263c6e6fe6eec2ef39663e6356d22f09f675622f893de7808f4664af2`
